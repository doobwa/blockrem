#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>
using namespace std;

int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

// Update each s(t,i,j) vector with event (a,b).
// s: List of NxN matrices with named elements.  Each matrix represents the current value for that statistic.

vector< vector< vector<int> > > initializeStatistics(int N, int P) {
  vector< vector< vector<int> > > s;
  for (int i = 0; i < N; i++) {
    vector< vector<int> > s_i;
    for (int j = 0; j < N; j++) {
      //      int arr[P] = {0};
      //      vector<int> s_ij (arr, arr + sizeof(arr) / sizeof(arr[0]) );
      vector<int> s_ij(P);
      s_i.push_back(s_ij);
    }
    s.push_back(s_i);
  }
  return s;
}


void updateStatistics(  vector< vector< vector<int> > > &s, int a, int b, int i, int j) {
  // (a,b) the event that occurred.  (i,j) the event that we are computing statistics for
  if (i != j) { 
    // P-shifts
    s[i][j][1] = (i!=a & i==b & j==a & j!=b); // abba
    s[i][j][2] = (i!=a & i==b & j!=a & j!=b); // abby
    s[i][j][3] = (i!=a & i!=b & j==a & j!=b); // abxa
    s[i][j][4] = (i!=a & i!=b & j!=a & j==b); // abxb
    s[i][j][5] = (i==a & i!=b & j!=a & j!=b); // abay
    s[i][j][6] = (i==a & i!=b & j!=a & j==b); // abab

    // Degree effects
    if (a==i) {// && b!=j) {
      s[i][j][7] += 1;  // sender out degree
    }
    if (a==j) {// && a!=j) {
      s[i][j][8] += 1;  // receiver out degree
    }
    if (b==i) {// && b!=i) { 
      s[i][j][9] += 1;  // sender in degree
    }
    if (b==j) {// && a!=i) {
      s[i][j][10] += 1; // receiver in degree
    }
  }
  // s[i][j][7] += 1; // sender out degree
  // s[i][j][8] += 1; // 
  // s[i][j][9] += 1;
  // s[i][j][10] += 1;
}

double computeLambda(int i, int j, int zi, int zj, vector<int> s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)]; // intercept
  for (int p = 1; p < P; p++) {
    lam += s[p] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}
  
// Compute a data structure for finding m_{last changepoint of ij}.
// Returns tau, where tau[i][j] is a vector of event indices m where lambda_{ij} changed (due to an event involving either i or j).  All vectors begin with 0 (since all intensities are assumed to change at time 0).  Element v of stats[i][j] are the statistics that were applicable up to event tau[i][j][v].  TODO: Should also have M-1?

vector< vector< vector<int> > > precomputeTauDyad(Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int N, int M) {
  vector< vector< vector<int> > > tau;
  for (int i = 0; i < N; i++) {
    vector< vector<int> > tau_i;
    for (int j = 0; j < N; j++) {
      int arr[] = {0};
      vector<int> tau_ij (arr, arr + sizeof(arr) / sizeof(arr[0]) );
      tau_i.push_back(tau_ij);
    }
    tau.push_back(tau_i);
  }
  for (int m = 0; m < M; m++) {
    int i = sen[m];
    int j = rec[m];
    for (int r = 0; r < N; r++) {
      if (r!=j && r!=i) {
        tau[i][r].push_back(m);
        tau[r][i].push_back(m);
        tau[j][r].push_back(m);
        tau[r][j].push_back(m);
      }
    }
    tau[i][j].push_back(m);
    tau[j][i].push_back(m);
  }
  return tau;
}

// element (i,j,v) is a vector of sufficient statistics for dyad (i,j) at its v'th changepoint.  These sufficient statistics apply to the time period leading up to the v'th timepoint.

vector< vector< vector< vector<int> > > > precomputeStats(Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int N, int M, int P) {
  vector< vector< vector<int> > > s = initializeStatistics(N,P);
  vector< vector< vector< vector<int> > > > x;
  for (int i = 0; i < N; i++) {
    vector< vector< vector<int> > > x_i;
    for (int j = 0; j < N; j++) {
      vector< vector<int> > x_ij;
      vector<int> x_ijp(P);
      x_ij.push_back(x_ijp);
      x_i.push_back(x_ij);
    }
    x.push_back(x_i);
  }
  for (int m = 0; m < M; m++) {
    int i = sen[m];
    int j = rec[m];
    for (int r = 0; r < N; r++) {
      if (r!=j && r!=i) {
        updateStatistics(s,i,j,i,r);
        updateStatistics(s,i,j,r,i);
        updateStatistics(s,i,j,j,r);
        updateStatistics(s,i,j,r,j);
        x[i][r].push_back(s[i][r]);
        x[r][i].push_back(s[r][i]);
        x[j][r].push_back(s[j][r]);
        x[r][j].push_back(s[r][j]);
      }
    }
    updateStatistics(s,i,j,i,j);
    updateStatistics(s,i,j,j,i);
    x[i][j].push_back(s[i][j]);
    x[j][i].push_back(s[j][i]);
  }
  return x;
}

// Get previous event index for a given vector of indices.  
// e.g. indx = [0 1 3 4].  getTau(indx,2)=1 and getTau(indx,1)=0.

int getLastChangepoint(vector<int> indx, int m) {
  if (m==0) {
    return 0;
  } else {
    vector<int>::iterator x = std::lower_bound(indx.begin(), indx.end(), m);
    int ix = int(x- indx.begin()) - 1;
    return indx[ix];
  }
}

int getLastChangepointIndex(vector<int> indx, int m) {
  if (m==0) {
    return 0;
  } else {
    vector<int>::iterator x = std::lower_bound(indx.begin(), indx.end(), m);
    return int(x - indx.begin()) - 1;
  }  
}

// (i,j,v) represents the statistic vector that

// Compute the loglikelihood corresponding to a single actor, a.

Rcpp::NumericVector llki(int a, Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P, vector<int> ma, vector< vector< vector<int> > > tau) {

  int i,j,zi,zj,m,t,v;
  double llk,lam;
  llk=0;
  Rcpp::NumericVector llks(M);//ma.size());
  vector< vector< vector< vector<int> > > > stats = precomputeStats(times,sen,rec,N,M,P);
  for (int m = 0; m < M; m++) {
    i = sen[m];
    j = rec[m];
    llk = 0;
    zi = z[i];
    zj = z[j];
    if (i==a | j==a | m==(M-1)) {
      //      llk += computeLambda(i,j,zi,zj,stats[i][j][t],beta,N,K,P);
      for (int r = 0; r < N; r++) {
        int zr = z[r];
        if (r != i) {
          v    = getLastChangepointIndex(tau[i][r],m);
          t    = getLastChangepoint(tau[i][r],m);
          lam  = computeLambda(i,r,zi,zr,stats[i][r][v],beta,N,K,P);
          llk -= (times[m] - times[t]) * exp(lam);
          //          Rprintf("%i (%i,%i) %f %i\n",m,i,r,lam,t);
          v    = getLastChangepointIndex(tau[r][i],m);
          t    = getLastChangepoint(tau[r][i],m);
          lam  = computeLambda(r,i,zr,zi,stats[r][i][v],beta,N,K,P);
          llk -= (times[m] - times[t]) * exp(lam);
          //          Rprintf("%i (%i,%i) %f %i\n",m,r,i,lam,t);
        }
        if (r != j) {
          v    = getLastChangepointIndex(tau[j][r],m);
          t    = getLastChangepoint(tau[j][r],m);
          lam  = computeLambda(j,r,zj,zr,stats[j][r][v],beta,N,K,P);
          llk -= (times[m] - times[t]) * exp(lam);
          v    = getLastChangepointIndex(tau[r][j],m);
          t    = getLastChangepoint(tau[r][j],m);
          lam  = computeLambda(r,j,zr,zj,stats[r][j][v],beta,N,K,P);
          llk -= (times[m] - times[t]) * exp(lam);
          //          Rprintf("%i (%i,%i) %f %i\n",m,r,j,lam,t);
        }
      }
    }
    llks(m) = llk;
  }
  return llks;
}

Rcpp::List gibbs(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P, Rcpp::List mas) {
  Rcpp::NumericVector x;
  Rcpp::List llks;
  vector< vector< vector<int> > > tau = precomputeTauDyad(times,sen,rec,N,M);
  for (int a=0; a < N; a++) {
    Rcpp::List llk_a;
    for (int k=0; k < K; k++) {
      z[a] = k;
      vector<int> ma = mas[a];
      x = llki(a,beta,times,sen,rec,z,N,M,K,P,ma,tau);
      llk_a.push_back(x);
    }
    llks.push_back(llk_a);
  }
  return Rcpp::List::create(Rcpp::Named("llks") = llks,
                            Rcpp::Named("z")    = z);
}
RCPP_MODULE(bremf){
  function( "gibbs", &gibbs ) ;
  function( "initializeStatistics", &initializeStatistics);
  function( "precomputeStats", &precomputeStats);
  //  function( "computeLambda", &computeLambda);
  function( "precomputeTauDyad", &precomputeTauDyad);
  //  function( "getTau", &getTau);
  function( "getLastChangepoint", &getLastChangepoint);
  function( "getLastChangepointIndex", &getLastChangepointIndex);
}
