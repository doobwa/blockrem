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


class Stat {
public:
  Stat(Rcpp::NumericVector times_, Rcpp::IntegerVector sen_, 
        Rcpp::IntegerVector rec_, int N_, int M_, int P_) : 
    times(times_),sen(sen_),rec(rec_),N(N_),M(M_),P(P_) {

    // Current vectors of statistics
    for (int i = 0; i < N; i++) {
      vector< vector<int> > s_i;
      for (int j = 0; j < N; j++) {
        vector<int> s_ij(P);
        s_i.push_back(s_ij);
      }
      s.push_back(s_i);
    }

    // Actor changepoitn indices
    for (int i = 0; i < N; i++) {
      vector<int> u_i;
      u.push_back(u_i);
    }

    // Dyad changepoint indices
    for (int i = 0; i < N; i++) {
      vector< vector<int> > v_i;
      for (int j = 0; j < N; j++) {
        int arr[] = {0};
        vector<int> v_ij (arr, arr + sizeof(arr) / sizeof(arr[0]) );
        v_i.push_back(v_ij);
      }
      v.push_back(v_i);
    }
    
    // Dyadic level statistics at each changepoing
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
  }


  void update(int a, int b, int i, int j) {
    // (a,b) the event that occurred.  
    // (i,j) the event that we are computing statistics for
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
  }


// element (i,j,v) is a vector of sufficient statistics for dyad (i,j) at its v'th changepoint.  These sufficient statistics apply to the time period leading up to the v'th timepoint.
// Compute a data structure for finding m_{last changepoint of ij}.
// Returns v, where v[i][j] is a vector of event indices m where lambda_{ij} changed (due to an event involving either i or j).  All vectors begin with 0 (since all intensities are assumed to change at time 0).  Element v of stats[i][j] are the statistics that were applicable up to event v[i][j][v].  TODO: Should also have M-1?

  void precompute() {
    for (int m = 0; m < M; m++) {
      Rprintf(".");
      int i = sen[m];
      int j = rec[m];
      for (int r = 0; r < N; r++) {
        if (r!=j && r!=i) {
          update(i,j,i,r);
          update(i,j,r,i);
          update(i,j,j,r);
          update(i,j,r,j);
          x[i][r].push_back(s[i][r]);
          x[r][i].push_back(s[r][i]);
          x[j][r].push_back(s[j][r]);
          x[r][j].push_back(s[r][j]);
          v[i][r].push_back(m);
          v[r][i].push_back(m);
          v[j][r].push_back(m);
          v[r][j].push_back(m);
        }
      }
      update(i,j,i,j);
      update(i,j,j,i);
      x[i][j].push_back(s[i][j]);
      x[j][i].push_back(s[j][i]);
      v[i][j].push_back(m);
      v[j][i].push_back(m);
      u[i].push_back(m);
      u[j].push_back(m);
    }
    for (int i = 0; i < N; i++) {
      u[i].push_back(M-1);
    }
    Rprintf("\n");
  }

  vector<int> get_s(int m, int i, int j) {
    int r = get_prev_index(v[i][j],m);
    return x[i][j][r];
  }

  double get_tau(int m, int i, int j) {
    int r = get_prev(v[i][j],m);
    return times[r];
  }

  vector< vector< vector<int> > >  get_all_v() {
    return v;
  }

  vector< vector< vector< vector<int> > > > get_all_s() {
    return x;
  }

  vector< vector<int> > get_all_u() {
    return u;
  }

  vector<int> get_u(int a) {
    return u[a];
  }

// Get previous event index for a given vector of indices.  
// e.g. indx = [0 1 3 4].  get_prev(indx,2)=1 and get_prev(indx,1)=0.

  int get_prev(vector<int> indx, int m) {
    if (m==0) {
      return 0;
    } else {
      return indx[get_prev_index(indx,m)];
    }
  }

  int get_prev_index(vector<int> ms, int m) {
    if (m==0) {
      return 0;
    } else {
      vector<int>::iterator low;
      sort (ms.begin(), ms.end());                // 10 10 10 20 20 20 30 30
      low = lower_bound (ms.begin(), ms.end(), m); //          &
      int ans = low - ms.begin();
      // Rprintf("%i\n",ans);
      // Rprintf("%i\n",ans-1);
      return ans - 1;
    }
  }

  Rcpp::IntegerVector get_prev_index2(vector<int> ms, int m) {
    Rcpp::IntegerVector a(get_prev_index(ms,m));
    return a;
  }
  SEXP ptr() {
    //Rprintf("%i",this);
    return wrap(XPtr<Stat>(this, true));
  }

  // Number of nodes, events, parameters, and clusters.
  int N, M, P;
  Rcpp::NumericVector times;
  Rcpp::IntegerVector sen;
  Rcpp::IntegerVector rec;

private:

  // Temporary stats built up incrementally
  vector< vector< vector<int> > > s;  

  // Data structure of stats vectors for each dyad at each of its changepoints
  vector< vector< vector< vector<int> > > > x;

  // Data structure of event indices for dyad's changepoint.
  vector< vector< vector<int> > > v;

  // Data structure of event indices for dyad's changepoint.
  vector< vector<int> > u;
};




double computeLambda(int i, int j, int zi, int zj, vector<int> s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)]; // intercept
  for (int p = 1; p < P; p++) {
    lam += s[p] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}


// (i,j,v) represents the statistic vector that

// Compute the loglikelihood corresponding to a single actor, a.

Rcpp::NumericVector llki(int a, Rcpp::NumericVector beta, Rcpp::IntegerVector z, Stat *s, int K) {
  int N = s->N;
  int M = s->M;
  int P = s->P;
  int i,j,zi,zj,m,t,v;
  double llk,lam;
  vector<int> smij;
  vector<int> ma = s->get_u(a);
  Rcpp::NumericVector llks(ma.size());
  for (int ix = 0; ix < ma.size(); ix++) {
    m = ma[ix];
    i = s->sen[m];
    j = s->rec[m];
    llk = 0;
    zi = z[i];
    zj = z[j];
    if (i==a | j==a | m==(M-1)) {
      smij = s->get_s(m,i,j);
      llk += computeLambda(i,j,zi,zj,smij,beta,N,K,P);
      for (int r = 0; r < N; r++) {
        int zr = z[r];
        if (r != i) {
          lam  = computeLambda(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
          llk -= (s->times[m] - s->get_tau(m,i,r)) * exp(lam);
          //          Rprintf("%i (%i,%i) %f %i\n",m,i,r,lam,t);
          lam  = computeLambda(r,i,zr,zi,s->get_s(m,r,i),beta,N,K,P);
          llk -= (s->times[m] - s->get_tau(m,r,i)) * exp(lam);
          //          Rprintf("%i (%i,%i) %f %i\n",m,r,i,lam,t);
        }
        if (r != j) {
          lam  = computeLambda(j,r,zj,zr,s->get_s(m,j,r),beta,N,K,P);
          llk -= (s->times[m] - s->get_tau(m,j,r)) * exp(lam);
          lam  = computeLambda(r,j,zr,zj,s->get_s(m,r,j),beta,N,K,P);
          llk -= (s->times[m] - s->get_tau(m,r,j)) * exp(lam);
          // //          Rprintf("%i (%i,%i) %f %i\n",m,r,j,lam,t);
        }
      }
    }
    llks(ix) = llk;
  }
  return llks;
}

Rcpp::List gibbs(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K) {
  Rcpp::NumericVector x;
  Rcpp::List llks;
  Stat *s = XPtr<Stat>(statptr_);
  for (int a = 0; a < s->N; a++) {
    Rcpp::List llk_a;
    // Compute p(z[a] = k | all else) for each k
    for (int k = 0; k < K; k++) {
      z[a] = k;
      x = llki(a,beta,z,s,K);
      llk_a.push_back(x);
    }
    llks.push_back(llk_a);
    Rprintf(".");
    // Sample z[a]
  }
    Rprintf("\n");
  return Rcpp::List::create(Rcpp::Named("llks") = llks,Rcpp::Named("z") = z);
}
RCPP_MODULE(bremf){
  //  function( "computeLambda", &computeLambda);
  //  function( "getTau", &getTau);
  class_<Stat>( "Stat" )
    .constructor<Rcpp::NumericVector,Rcpp::IntegerVector,Rcpp::IntegerVector,
                 int,int,int>()
    //    .property( "x", &Stat::get_x, &Stat::set_x,
    //               "Docstring for x" )
    .method( "precompute", &Stat::precompute,
             "Precompute the data structure of REM statistics")
    .method( "get_s", &Stat::get_s,
             "")
    .method( "get_tau", &Stat::get_tau,
             "")
    .method( "get_tau", &Stat::get_u,
             "")
    .method( "get_all_s", &Stat::get_all_s,
             "")
    .method( "get_all_v", &Stat::get_all_v,
             "")
    .method( "get_all_u", &Stat::get_all_u,
             "")
    .method( "get_prev", &Stat::get_prev,
             "")
    .method( "get_prev_index", &Stat::get_prev,
             "")
    .method( "get_prev_index2", &Stat::get_prev,
             "")
    .method( "ptr", &Stat::ptr,
             "")
    ;
  function( "gibbs", &gibbs ) ;
}
