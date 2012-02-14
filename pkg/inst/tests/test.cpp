/* 
Christopher DuBois

C++ code for block model relational events.  There are currently two ways of fitting such models with the following code.  The first computes the likelihood of the model by updating statistic vectors after each event.  These functions include:
llk: compute the likelihood
llk2: compute the likelihood using a log rate array
lrm:  compute a log rate array (T x N x N)
initializeStatistics: initialize an P x N x N
updateStatistics: update the current most recent set of statistics with an observed event.
computeLambda: compute the log intentisty function, lambda, using a vector of statistics and parameters beta.

The second method precomputes the statistic vector for every dyad at every observed event.  Though this requires a large amount of memory, searching for the last observed time will not be required.  The RemStat class encapsulates this data structure and allows one to precompute all the necessary statistics given the observed data, as well as query for s(t,i,j) and tau(t,i,j).  The folowing functions use a RemStat object.

llkfast: 
lrmfast
gibbs
computeLambdaFast

rcategorical is a helper function for gibbs.  It draws from a categorical distribution when the log probabilities are given.
 */


#include <iostream>
#include <algorithm>
#include <vector>
//#include <omp.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;
RNGScope scope;

// Return the (j,k,l) element of a (J,K,L) array represented as a vector.
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

// Update each s(t,i,j) vector with event (a,b).
// s: List of NxN matrices with named elements.  Each matrix represents the current value for that statistic.


class RemStat {
public:
  RemStat(Rcpp::NumericVector times_, Rcpp::IntegerVector sen_, 
          Rcpp::IntegerVector rec_, int N_, int M_, int P_) : 
    times(times_),sen(sen_),rec(rec_),N(N_),M(M_),P(P_) {

    // Current vectors of statistics
    for (int i = 0; i < N; i++) {
      vector< vector<int> > w_i;
      for (int j = 0; j < N; j++) {
        vector<int> w_ij(M);
        w_i.push_back(w_ij);
      }
      w.push_back(w_i);
    }

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
        //        int arr[] = {0};
        vector<int> v_ij(1);// (arr, arr + sizeof(arr) / sizeof(arr[0]) );
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
        //        x_ijp[0] = 1;
        x_ij.push_back(x_ijp);
        x_i.push_back(x_ij);
      }
      x.push_back(x_i);
    }
  }

  // Update the statistics for dyad (i,j) with the information that dyad (a,b) just occurred.  s_1{(i,j)} (the statistic for an abba effect)  will now be 1 if b==i and a==j.

  void update(int a, int b, int i, int j) {
    s[i][j][0] = 1;
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
    // Save nothing to v or x (as they were anlready initalized with 0)

    for (int m = 0; m < (M-1); m++) {
      //      Rprintf(".");

      int i = sen[m];
      int j = rec[m];
     // Update statistics for those affected by previous event
      for (int r = 0; r < N; r++) {
        if (r!=j && r!=i) {
          update(i,j,i,r);
          update(i,j,r,i);
          update(i,j,j,r);
          update(i,j,r,j);
        }
      }
      update(i,j,i,j);
      update(i,j,j,i);


      // Update statistics and changepoints for all dyads affected by (i,j)@m
      for (int r = 0; r < N; r++) {
        if (r!=i && r!=j) {
          x[i][r].push_back(s[i][r]); // pushing to m element
          x[r][i].push_back(s[r][i]);
          x[j][r].push_back(s[j][r]);
          x[r][j].push_back(s[r][j]);
          v[i][r].push_back(m);
          v[r][i].push_back(m);
          v[j][r].push_back(m);
          v[r][j].push_back(m);
        }
      }
      x[i][j].push_back(s[i][j]);
      x[j][i].push_back(s[j][i]);
      v[i][j].push_back(m);
      v[j][i].push_back(m);

      // Update index of changepoint for this dyad
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          w[i][j][m] = v[i][j].size() - 1;
        }
      }

      u[i].push_back(m);
      u[j].push_back(m);
 
    }

    // Handle last event separately: all v[i][j] get M-1 at the end.
    // Update w[i][j][M-1] to be the final size of v[i][j]
    int m = M-1;
    int a = sen[m];
    int b = rec[m];

    for (int i = 0; i < N; i++) {
      u[i].push_back(m);
      for (int j = 0; j < N; j++) {
        if (i != j) {
          update(a,b,i,j);
          x[i][j].push_back(s[i][j]);  // should go unused.
          v[i][j].push_back(m);
          w[i][j][m] = v[i][j].size() - 1;
        }
      }
    }
      
    Rprintf("\n");
  }

  // Get the statistic vector s_(i,j) that is used before event m.
  vector<int> get_s(int m, int i, int j) {
    int r = w[i][j][m - 1];
    return x[i][j][r];
  }

  // Get the time of the event that last affected dyad (i,j)
  double get_tau(int m, int i, int j) {
    // If asking for an event past the observed events, return the last time
    if (m > M-1) {
      m = M-1;
    }
    // If asking for the first (or negative) event, return the first time
    if (m < 1) {
      m = 1;
    }
    int r = w[i][j][m - 1];
    int y = v[i][j][r];
    return times[y];
  }

  // For each (i,j) a vector of event indicies for lambda_ij changepoints.
  // All vectors start with 0 (since all are assumed to have a changepoint at 0)
  vector< vector< vector<int> > >  get_all_v() {
    return v;
  }

  // For each (i,j,v_ij) a vector of P statistics.  First vector for each (i,j) is all 0 (the vector of statistics for that dyad before the first event occurs)
  vector< vector< vector< vector<int> > > > get_all_s() {
    return x;
  }

  vector<int> get_v(int i, int j) {
    return v[i][j];
  }
  vector<int> get_w(int i, int j) {
    return w[i][j];
  }

  vector< vector<int> > get_all_u() {
    return u;
  }

  // Vector of changepoints for a given actor
  vector<int> get_u(int a) {
    return u[a];
  }


  SEXP ptr() {
    //Rprintf("%i",this);
    //return wrap(XPtr<RemStat>(this, true));
    return Rcpp::XPtr<RemStat>(this, false);
  }

  // Number of nodes, events, parameters, and clusters.
  int N, M, P;
  Rcpp::NumericVector times;
  Rcpp::IntegerVector sen;
  Rcpp::IntegerVector rec;

private:

  // Temporary stats built up incrementally
  vector< vector< vector<int> > > s;  

  // Data structure of current u_{ijm} 
  vector< vector< vector<int> > > w;  

  // Data structure of stats vectors for each dyad at each of its changepoints
  vector< vector< vector< vector<int> > > > x;

  // Data structure of event indices for dyad's changepoint.
  vector< vector< vector<int> > > v;

  // Data structure of event indices for dyad's changepoint.
  vector< vector<int> > u;
};


// Use the precomputed set of statistics to copmute lambda

double computeLambdaFast(int i, int j, int zi, int zj, vector<int> s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)]; // intercept
  for (int p = 1; p < 7; p++) {
    lam += s[p] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  for (int p = 7; p < P; p++) {
    lam += log(s[p]+1) * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}

// Compute the loglikelihood corresponding to a single actor, a.

double llki(int a, Rcpp::NumericVector beta, Rcpp::IntegerVector z, RemStat *s, int K) {
  int N = s->N;
  int M = s->M;
  int P = s->P;
  int i,j,zi,zj,m,t,v;
  double llk;
  double total = 0;
  vector<int> smij;
  vector<int> ma = s->get_u(a);
  Rcpp::NumericVector llks(ma.size());
  // m = 1259;
  // i = 2;
  // v = 6;
  // double tmp = computeLambdaFast(i,v,z[i],z[v],s->get_s(m,i,v),beta,N,K,P);
  // Rprintf("%f\n",tmp);
  // tmp = (s->times[m] - s->get_tau(m,i,v)) * exp(tmp);
  // Rprintf("%f\n",tmp);
  for (int ix = 0; ix < ma.size(); ix++) {
    m = ma[ix];
    i = s->sen[m];
    j = s->rec[m];
    llk = 0;
    zi = z[i];
    zj = z[j];
    llk += computeLambdaFast(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);

    // #pragma omp parallel reduction(+:llk)
    // {
    // #pragma omp for

    for (int r = 0; r < N; r++) {
      double lam;
      int zr = z[r];
      if (r != i && r !=j) {
        Rprintf("%i %i %i %i %i %i %i %i:",a,m,i,j,r,zi,zj,zr);
        //        if (m==1259 && i==2 && j==3 && r==6 && a==2) {
        if (a==7 && m==1205 && i==7 && j==6 && r==5) {
          // Rprintf("%f ",s->get_tau(m,i,r));
          // Rprintf("%f ",s->get_tau(m,r,i));
          // Rprintf("%f ",s->get_tau(m,j,r));
          // Rprintf("%f ",s->get_tau(m,r,j));
          // Rprintf("%i ",(s->get_s(m,i,r))[0]);
          // Rprintf("%i ",(s->get_s(m,r,i))[0]);
          // Rprintf("%i ",(s->get_s(m,j,r))[0]);
          // Rprintf("%i\n",(s->get_s(m,r,j))[0]);
          lam  = computeLambdaFast(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
          Rprintf("%f ",lam); 
          llk -= (s->times[m] - s->get_tau(m,i,r)) * exp(lam);
          Rprintf("%f \n",llk);
          lam  = computeLambdaFast(r,i,zr,zi,s->get_s(m,r,i),beta,N,K,P);
          Rprintf("%f ",lam); 
          llk -= (s->times[m] - s->get_tau(m,r,i)) * exp(lam);
          Rprintf("%f \n",llk); 
          lam  = computeLambdaFast(j,r,zj,zr,s->get_s(m,j,r),beta,N,K,P);
          Rprintf("%f ",lam); 
          llk -= (s->times[m] - s->get_tau(m,j,r)) * exp(lam);
          Rprintf("%f \n",llk); 
          lam  = computeLambdaFast(r,j,zr,zj,s->get_s(m,r,j),beta,N,K,P);
          Rprintf("%f ",lam); 
          llk -= (s->times[m] - s->get_tau(m,r,j)) * exp(lam);
          Rprintf("%f \n",llk);
          }
        lam  = computeLambdaFast(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
        Rprintf("%f \n",lam); 
      }
    }

    //} // openmp

    llks(ix) = llk;
    total += llk;
  }
  return total;
}


 bool IsFiniteNumber(double x) 
    {
        return (x <= DBL_MAX && x >= -DBL_MAX); 
    } 
int rcategorical (Rcpp::NumericVector lp) {
  int max = 0;
  int min = 0;
  int K = lp.size();
  for (int i = 0; i < K; i++) {
    if (lp[i] > lp[max]) {
      max = i;
    }
    if (lp[i] < lp[min]) {
      min = i;
    }
  }
  
  // If all small/large values, subtract max/min value for numeric stability.
  if (lp[max] < 0) {
    double lmax = lp[max];
    for (int i = 0; i < K; i++) {
      lp[i] -= lmax;
    }
  }
  if (lp[min] > 0) {
    double lmin = lp[min];
    for (int i = 0; i < K; i++) {
      lp[i] -= lmin;
    }
  }

  // Inverse CDF method
  Rcpp::NumericVector p = exp(lp);

  // Get sum.  If any Infinite values, return that index. 
  // TODO: This probably introduces a bias towards smaller k.
  int k;
  double cuml = 0;
  double sum = 0;
  for (k=0;k<K;k++) {
    if (!IsFiniteNumber(p[k])) {
      return k;
    }
    sum += p[k];
  }

  double r = as<double>(Rcpp::runif(1));
  for (k=0;k<K;k++) {
    cuml += p[k]/sum;
    if (r < cuml) {
      break;
    } 
  }
  return k;
}


// TODO: sample z
Rcpp::List gibbs(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K) {
  Rcpp::NumericVector x;
  Rcpp::List llks;
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  Rcpp::IntegerVector counts;
  double alpha = 1.0;

  // Get counts for the number assigned each class
  for (int i = 0; i < N; i++) {
    counts[z[i]] += 1;
  }

  for (int a = 0; a < N; a++) {

    counts[z[a]] -= 1;
    Rcpp::NumericVector y(K);
    // Compute p(z[a] = k | all else) for each k
    for (int k = 0; k < K; k++) {
      z[a] = k;
      y[k] = llki(a,beta,z,s,K);// + log(counts[k] + alpha) - log(N + alpha);
    }
    llks.push_back(y);
    
    // Sample z[a]
    z[a] = rcategorical(y);
    counts[z[a]] += 1;
  }
  return Rcpp::List::create(Rcpp::Named("llks") = llks,Rcpp::Named("z") = z);
}


RCPP_MODULE(brem){
  
  class_<RemStat>( "RemStat" )
    .constructor<Rcpp::NumericVector,Rcpp::IntegerVector,Rcpp::IntegerVector,
                 int,int,int>()
    .method( "precompute", &RemStat::precompute,
             "Precompute the data structure of REM statistics")
    .method( "get_s", &RemStat::get_s,
             "retrieve the statistics vector prior to event m for dyad (i,j)")
    .method( "get_tau", &RemStat::get_tau,
             "retrieve the last time that lambda_(i,j) changed")
    // .method( "get_u", &Stat::get_u,
    //          "")
    .method( "get_all_s", &RemStat::get_all_s,
             "")
    .method( "get_all_v", &RemStat::get_all_v,
             "")
    .method( "get_all_u", &RemStat::get_all_u,
             "")
    .method( "get_v", &RemStat::get_v,
             "vector where element k is the event index of the k'th changepoint for (i,j).  i.e. if w(i,j)[m] = k then v[i,j,k] = m")
    .method( "get_w", &RemStat::get_w,
             "vector where element m is the index of the previous changepoint for (i,j).  i.e. if w(i,j)[m] = k then v[i,j,k] = m")
    .method( "ptr", &RemStat::ptr,
             "")
    ;
  function( "gibbs", &gibbs ) ;
  function( "loglikelihood_fast", &llkfast ) ;
  function( "compute_lambda_fast", &computeLambdaFast);
  function( "rcategorical", &rcategorical);
}
