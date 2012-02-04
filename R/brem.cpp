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

Rcpp::NumericVector initializeStatistics(int N, int P) {
  Rcpp::NumericVector s     = Rcpp::NumericVector(Dimension(P,N,N));

  //  Intercept statistic: all (0,i,j) are equal to 1
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      s[threeDIndex(0,i,j,P,N,N)] = 1;
    }
  }

  return s;
}


Rcpp::NumericVector updateStatistics(Rcpp::NumericVector s, int a, int b, int N, int P) {
  // Create vector of indicators for each pshift.
  // i.e. If I[(i,j) is ab-ba from last event (a,b)]
    // Iterate through dyads having either a or b as a sender or receiver
   for (int r = 0; r < N; r++) {
      Rcpp::IntegerVector sen = Rcpp::IntegerVector::create(a,r,b,r);
      Rcpp::IntegerVector rec = Rcpp::IntegerVector::create(r,a,r,b);
      for (int k = 0; k < sen.size(); k++) {
        int i = sen[k];
        int j = rec[k];
        if (i != j) {
          // P-shifts
          s[threeDIndex(1,i,j,P,N,N)] = (i!=a & i==b & j==a & j!=b);
          s[threeDIndex(2,i,j,P,N,N)] = (i!=a & i==b & j!=a & j!=b);
          s[threeDIndex(3,i,j,P,N,N)] = (i!=a & i!=b & j==a & j!=b);
          s[threeDIndex(4,i,j,P,N,N)] = (i!=a & i!=b & j!=a & j==b);
          s[threeDIndex(5,i,j,P,N,N)] = (i==a & i!=b & j!=a & j!=b);
          s[threeDIndex(6,i,j,P,N,N)] = (i==a & i!=b & j!=a & j==b);

          // Degree effects
          if (a==i & b!=j) {
            s[threeDIndex(7,i,j,P,N,N)] += 1;  // sender out degree
          }
          if (b==i & a!=j) {
            s[threeDIndex(8,i,j,P,N,N)] += 1;  // sender in degree
          }
          if (a==j & b!=i) { 
            s[threeDIndex(9,i,j,P,N,N)] += 1;  // receiver out degree
          }
          if (b==j & a!=i) {
            s[threeDIndex(10,i,j,P,N,N)] += 1; // receiver in degree
          }
        }
      }
    }
   s[threeDIndex(7,a,b,P,N,N)] += 1; // sender out degree
   s[threeDIndex(8,b,a,P,N,N)] += 1; // 
   s[threeDIndex(9,b,a,P,N,N)] += 1;
   s[threeDIndex(10,a,b,P,N,N)] += 1;
  return s;
}

void updateStatistics2(Rcpp::NumericVector &s, int a, int b, int i, int j, int N, int P) {
  // (a,b) the event that occurred.  (i,j) the event that we are computing statistics for
        if (i != j) {
          // P-shifts
          s[threeDIndex(1,i,j,P,N,N)] = (i!=a & i==b & j==a & j!=b);
          s[threeDIndex(2,i,j,P,N,N)] = (i!=a & i==b & j!=a & j!=b);
          s[threeDIndex(3,i,j,P,N,N)] = (i!=a & i!=b & j==a & j!=b);
          s[threeDIndex(4,i,j,P,N,N)] = (i!=a & i!=b & j!=a & j==b);
          s[threeDIndex(5,i,j,P,N,N)] = (i==a & i!=b & j!=a & j!=b);
          s[threeDIndex(6,i,j,P,N,N)] = (i==a & i!=b & j!=a & j==b);

          // Degree effects
          if (a==i & b!=j) {
            s[threeDIndex(7,i,j,P,N,N)] += 1;  // sender out degree
          }
          if (b==i & a!=j) {
            s[threeDIndex(8,i,j,P,N,N)] += 1;  // sender in degree
          }
          if (a==j & b!=i) { 
            s[threeDIndex(9,i,j,P,N,N)] += 1;  // receiver out degree
          }
          if (b==j & a!=i) {
            s[threeDIndex(10,i,j,P,N,N)] += 1; // receiver in degree
          }
        }
   s[threeDIndex(7,a,b,P,N,N)] += 1; // sender out degree
   s[threeDIndex(8,b,a,P,N,N)] += 1; // 
   s[threeDIndex(9,b,a,P,N,N)] += 1;
   s[threeDIndex(10,a,b,P,N,N)] += 1;
  return;
}

//
double computeLambda(int i, int j, int zi, int zj, Rcpp::NumericVector s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = 0;
  for (int p = 0; p < P; p++) {
    lam += s[threeDIndex(p,i,j,P,N,N)] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}
double computeLambda2(int i, int j, int zi, int zj, Rcpp::NumericVector &s, Rcpp::NumericVector &beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)]; // intercept
  for (int p = 1; p < P; p++) {
    lam += s[threeDIndex(p,i,j,P,N,N)] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}

Rcpp::NumericVector llkp(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M, int K, int P) {
  omp_set_num_threads(16);
  double lam = 0;
  int i,j,r,zi,zj;
  double llktotal = 0.0;
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N);

  Rcpp::NumericVector s  = initializeStatistics(N,P);
  s = updateStatistics(s,sen[0],rec[0],N,P);
  updateStatistics2(s,sen[0],rec[0],0,0,N,P);
  Rcpp::NumericVector llk(M);
  Rcpp::NumericVector llks(N);
  
  for (int m = 1; m < (M-1); m++) {
    i = sen[m];
    j = rec[m];
    zi = z[i];
    zj = z[j];
    //llk(m) += computeLambda2(i,j,zi,zj,s,beta,N,K,P);
    
    double llkm = computeLambda2(i,j,zi,zj,s,beta,N,K,P);
    #pragma omp parallel reduction(-:llkm)
    {
    #pragma omp for
      for (r = 0; r < N; r++) {
        int zr = z[r];
        if (r != i) {
          lam  = computeLambda2(i,r,zi,zr,s,beta,N,K,P);
          llkm -= (times[m] - times[mp(i,r)]) * exp(lam);
          lam  = computeLambda2(r,i,zr,zi,s,beta,N,K,P);
          llkm -= (times[m] - times[mp(r,i)]) * exp(lam);
          mp(i,r) = m;
          mp(r,i) = m;
          updateStatistics2(s,sen[m],rec[m],i,r,N,P);
          updateStatistics2(s,sen[m],rec[m],r,i,N,P);
        }
        if (r != j) {
          lam  = computeLambda2(j,r,zj,zr,s,beta,N,K,P);
          llkm -= (times[m] - times[mp(j,r)]) * exp(lam);
          lam  = computeLambda2(r,j,zr,zj,s,beta,N,K,P);
          llkm -= (times[m] - times[mp(r,j)]) * exp(lam);
          mp(j,r) = m;  // update mp
          mp(r,j) = m;
          updateStatistics2(s,sen[m],rec[m],j,r,N,P);
          updateStatistics2(s,sen[m],rec[m],r,j,N,P);
        }
      }
    } // openmp
    //llk(m) = std::accumulate(llks.begin(),llks.end(), 0.0);
  }

  //llktotal = std::accumulate(llk.begin(),llk.end(), 0.0);
  return llk;
}

// Rcpp::IntegerVector precomputeTau(Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec) {
//   Rcpp::IntegerMatrix tau = Rcpp::IntegerMatrix(Dimension(M,N,N));
//   for (int m = 0; m < times.size(); m++) {
//     tau[threeDIndex(m,sen[m],rec[m],M,N,N)] = m;
//   }
//   return tau;
// }
  
// Compute a data structure for finding tau_{ijm}.
// Returns tau, where tau[i][j] is a vector of event indices m where lambda_{ij} changed (due to an event involving either i or j.

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
      tau[i][r].push_back(m);
      tau[r][i].push_back(m);
      tau[j][r].push_back(m);
      tau[r][j].push_back(m);
    }
  }
  // TODO: Add M-1 to everybody?  Assumption: all intensities change on last event....
  return tau;
}

// Get last event for a given vector of times
int getTau(vector<int> indx, int m) {
  if (m==0) {
    return 0;
  } else {
    vector<int>::iterator x = std::lower_bound(indx.begin(), indx.end(), m);
    int ix = int(x- indx.begin()) - 1;
    return indx[ix];
  }
}

Rcpp::NumericVector llki(int a, Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P, Rcpp::IntegerVector ma, vector< vector< vector<int> > > tau) {

  int i,j,zi,zj,m,t;
  double llk,lam;
  llk=0;
  Rcpp::NumericVector llks(M);//ma.size());
  Rcpp::NumericVector s  = Rcpp::NumericVector(Dimension(P,N,N));  
  for (int m = 0; m < M; m++) {
    i = sen[m];
    j = rec[m];
    llk = 0;
    zi = z[i];
    zj = z[j];
    if (i==a | j==a | m==(M-1)) {
      llk += computeLambda2(i,j,zi,zj,s,beta,N,K,P);
      for (int r = 0; r < N; r++) {
        int zr = z[r];
        if (r != i) {
          lam  = computeLambda2(i,r,zi,zr,s,beta,N,K,P);
          t    = getTau(tau[i][r],m);
          //          Rprintf("%i (%i,%i) %f %i\n",m,i,r,lam,t);
          llk -= (times[m] - times[t]) * exp(lam);
          lam  = computeLambda2(r,i,zr,zi,s,beta,N,K,P);
          t    = getTau(tau[r][i],m);
          //          Rprintf("%i (%i,%i) %f %i\n",m,r,i,lam,t);
          llk -= (times[m] - times[t]) * exp(lam);
        }
        if (r != j) {
          lam  = computeLambda2(j,r,zj,zr,s,beta,N,K,P);
          t    = getTau(tau[j][r],m);
          //          Rprintf("%i (%i,%i) %f %i\n",m,j,r,lam,t);
          llk -= (times[m] - times[t]) * exp(lam);
          lam  = computeLambda2(r,j,zr,zj,s,beta,N,K,P);
          t    = getTau(tau[r][j],m);
          //          Rprintf("%i (%i,%i) %f %i\n",m,r,j,lam,t);
          llk -= (times[m] - times[t]) * exp(lam);
        }
      }
    }
    s = updateStatistics(s,i,j,N,P);  // incorrect for dyads not involving a
    llks(m) = llk;
  }
  return llks;
}



Rcpp::List gibbs(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P, Rcpp::List mas) {
  Rcpp::NumericVector x;
  Rcpp::List llks;
  vector< vector< vector<int> > > tau = precomputeTauDyad(times,sen,rec,N,M);
  for (int a=0; a < 1; a++) {
    Rcpp::List llk_a;
    for (int k=0; k < 1; k++) {
      z[a] = k;
      Rcpp::IntegerVector ma = mas[a];
      x = llki(a,beta,times,sen,rec,z,N,M,K,P,ma,tau);
      llk_a.push_back(x);
    }
    llks.push_back(llk_a);
  }
  return Rcpp::List::create(Rcpp::Named("llks") = llks,
                            Rcpp::Named("z")    = z);
}

Rcpp::NumericVector llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P) {


  double llk = 0.0; 
  double lam = 0;
  Rcpp::NumericVector llks(M);
  llks[0] = llk;

  int i,j,r;

  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N);
  Rcpp::NumericVector s  = initializeStatistics(N,P);
  s = updateStatistics(s,sen[0],rec[0],N,P);

  for (int m = 1; m < (M-1); m++) {
    i = sen[m];
    j = rec[m];
    int zi = z[i];
    int zj = z[j];
    llk += computeLambda2(i,j,zi,zj,s,beta,N,K,P);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int r = 0; r < N; r++) {
      int zr = z[r];
      if (r != i) {
        lam  = computeLambda2(i,r,zi,zr,s,beta,N,K,P);
        llk -= (times[m] - times[mp(i,r)]) * exp(lam);
        lam  = computeLambda2(r,i,zr,zi,s,beta,N,K,P);
        llk -= (times[m] - times[mp(r,i)]) * exp(lam);
        mp(i,r) = m;
        mp(r,i) = m;
      }
      if (r != j) {
        lam  = computeLambda2(j,r,zj,zr,s,beta,N,K,P);
        llk -= (times[m] - times[mp(j,r)]) * exp(lam);
        lam  = computeLambda2(r,j,zr,zj,s,beta,N,K,P);
        llk -= (times[m] - times[mp(r,j)]) * exp(lam);
        mp(j,r) = m;  // update mp
        mp(r,j) = m;
      }
    }
    s = updateStatistics(s,sen[m],rec[m],N,P);
    llks[m] = llk;
  }
  // All intensities assumed to change at the last event
  //for (int i = 0; i < N; i++) {
  //  for (int j = 0; j < N; j++) {
  //    int zi = z[i];
  //    int zj = z[j];
  //    if (i != j) {
  //      lam  = computeLambda2(i,j,zi,zj,s,beta,N,K,P);
  //      llk -= (times[M-1] - times[mp(i,j)]) * exp(lam);
  //    }
  //  }
  //}
  llks[M-1] = llk;
  return llks;
}



// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).

Rcpp::NumericVector lrm(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec,Rcpp::IntegerVector z, int N, int M,int K, int P){

  Rcpp::NumericVector lrmat = Rcpp::NumericVector(Dimension(M,N,N));
  Rcpp::NumericVector s = initializeStatistics(N,P);

  //int a,b,u,v,i,j;
  for (int m = 1; m < M; m++) {
    s = updateStatistics(s,sen[m-1],rec[m-1],N,P);
    int i = 1;
    int j = 0;
    int zi = z[i];
    int zj = z[j];
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int zi = z[i];
        int zj = z[j];
        lrmat[threeDIndex(m,i,j,M,N,N)] = computeLambda(i,j,zi,zj,s,beta,N,K,P);
      }
    }
  }
  return lrmat;
}
// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).
double llk2(Rcpp::NumericVector lrm,
            Rcpp::NumericVector times,
            Rcpp::IntegerVector sen, 
            Rcpp::IntegerVector rec, 
            int N, int M){
  double llk = 0;
  double den;
  double delta = 0;
  int i,j;
  llk = lrm[threeDIndex(0,sen[0],rec[0],M,N,N)];
  for (int m = 1; m < M; m++) {
    den = 0.0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (i != j) {
          den += exp(lrm[threeDIndex(m,i,j,M,N,N)]);
        }
      }
    }
    llk += lrm[threeDIndex(m,sen[m],rec[m],M,N,N)];
    delta = times[m] - times[m-1];
    llk -= delta * den;
  }
  return llk;
}
RCPP_MODULE(brem){
  function( "llk", &llk ) ;
  function( "llkp", &llkp ) ;
  //  function( "llki", &llki ) ;
  function( "gibbs", &gibbs ) ;
  function( "llk2", &llk2 ) ;
  function( "lrm", &lrm ) ;
  function( "updateStatistics", &updateStatistics);
  function( "initializeStatistics", &initializeStatistics);
  function( "computeLambda", &computeLambda);
  function( "precomputeTauDyad", &precomputeTauDyad);
  function( "getTau", &getTau);
}
