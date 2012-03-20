/* 
Christopher DuBois

C++ code for block model relational events.  There are currently two ways of fitting such models with the following code.  The first computes the likelihood of the model by updating statistic vectors after each event.  These functions include:
Loglikelihood: compute the likelihood
LoglikelihoodFromLogRateArray: compute the likelihood using a log rate array
Lrm:  compute a log rate array (T x N x N)
InitializeStatistics: initialize an P x N x N
UpdateStatistics: update the current most recent set of statistics with an observed event.
ComputeLambda: compute the log intentisty function, lambda, using a vector of statistics and parameters beta.

The second method precomputes the statistic vector for every dyad at every observed event.  Though this requires a large amount of memory, searching for the last observed time will not be required.  The RemStat class encapsulates this data structure and allows one to precompute all the necessary statistics given the observed data, as well as query for s(t,i,j) and tau(t,i,j).  The folowing functions use a RemStat object.

llkfast: 
lrmfast
gibbs
computeLambdaFast

Currently implemented statistics for a particular (t,i,j):
- Pshifts: ab-ba, ab-by, ab-ay, ab-xa, ab-xb, ab-ab
- Degree effects: 
  - sender out degree
  - sender in degree
  - receiver out degree
  - receiver in degree
  - number of times this dyad has occurred
  - number of changepoints fo this dyad

old notes:
// element (i,j,v) is a vector of sufficient statistics for dyad (i,j) at its v'th changepoint.  These sufficient statistics apply to the time period leading up to the v'th timepoint.
// Compute a data structure for finding m_{last changepoint of ij}.
// Returns v, where v[i][j] is a vector of event indices m where lambda_{ij} changed (due to an event involving either i or j).  All vectors begin with 0 (since all intensities are assumed to change at time 0).  Element v of stats[i][j] are the statistics that were applicable up to event v[i][j][v].  TODO: Should also have M-1?


rcategorical is a helper function for gibbs.  It draws from a categorical distribution when the log probabilities are given.
 */


#include <iostream>
#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include <utils.h>
using namespace std;
using namespace Rcpp;
RNGScope scope;

/*
Data structure for precomputing relational event statistics
*/

class RemStat {
public:
  RemStat(Rcpp::NumericVector times_, Rcpp::IntegerVector sen_,  Rcpp::IntegerVector rec_, int N_, int M_, int P_) : 
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

  void update(int m, int a, int b, int i, int j) {
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
      if (a==i && b==j) {
        s[i][j][11] += 1; // dyad count
      }
      s[i][j][12] = m;   // changepoint count
    }
  }


  void precompute() {
    // Save nothing to v or x (as they were anlready initalized with 0)

    for (int m = 0; m < (M-1); m++) {
      //      Rprintf(".");

      int i = sen[m];
      int j = rec[m];
     // Update statistics for those affected by previous event
      for (int r = 0; r < N; r++) {
        if (r!=j && r!=i) {
          update(m,i,j,i,r);
          update(m,i,j,r,i);
          update(m,i,j,j,r);
          update(m,i,j,r,j);
        }
      }
      update(m,i,j,i,j);
      update(m,i,j,j,i);


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
          update(m,a,b,i,j);
          x[i][j].push_back(s[i][j]);  // should go unused.
          v[i][j].push_back(m);
          w[i][j][m] = v[i][j].size() - 1;
        }
      }
    }
  }

  // Get the statistic vector s_(i,j) that is used before event m.
  vector<int> get_s(int m, int i, int j) {
    int r = w[i][j][m - 1];
    return x[i][j][r];
  }

  // Get the time of the event previous to m that affected dyad (i,j)
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

  // Return all changepoints.  Each element corresponds to 
  vector< vector<int> > get_all_u() {
    return u;
  }

  // Return vector of changepoints for a given actor
  vector<int> get_u(int a) {
    return u[a];
  }

  // Return 
  SEXP ptr() {
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

double LogLambdaPc(int i, int j, int zi, int zj, vector<int> s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)]; // intercept
  for (int p = 1; p < 7; p++) {
    lam += s[p] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  double numEvents = double(s[12]);
  for (int p = 7; p < 12; p++) {
    //    lam += s[p]/(numEvents + 1) * beta[threeDIndex(p,zi,zj,P,K,K)];
    lam += log((s[p]+1)/(numEvents + N*(N-1))) * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}

// Compute the loglikelihood corresponding to a single actor, a.

double RemLogLikelihoodActorPc(int a, Rcpp::NumericVector beta, Rcpp::IntegerVector z, RemStat *s, int K) {
  int N = s->N;
  int M = s->M;
  int P = s->P;
  int i,j,zi,zj,m,t,v;
  double llk;
  double total = 0;
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
    llk += LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);

    // #pragma omp parallel reduction(+:llk)
    // {
    // #pragma omp for

    for (int r = 0; r < N; r++) {
      double lam;
      int zr = z[r];
      if (r != i && r !=j) {
        lam  = LogLambdaPc(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
        llk -= (s->times[m] - s->get_tau(m,i,r)) * exp(lam);
        lam  = LogLambdaPc(r,i,zr,zi,s->get_s(m,r,i),beta,N,K,P);
        llk -= (s->times[m] - s->get_tau(m,r,i)) * exp(lam);
        lam  = LogLambdaPc(j,r,zj,zr,s->get_s(m,j,r),beta,N,K,P);
        llk -= (s->times[m] - s->get_tau(m,j,r)) * exp(lam);
        lam  = LogLambdaPc(r,j,zr,zj,s->get_s(m,r,j),beta,N,K,P);
        llk -= (s->times[m] - s->get_tau(m,r,j)) * exp(lam);
      }
    }

    //    } // openmp

    llks(ix) = llk;
    total += llk;
  }
  return total;
}

// Sample class assignments for each actor using RemLogLikelihoodActorPc

Rcpp::List RemGibbsPc(Rcpp::IntegerVector ix, Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K) {
  Rcpp::List llks;
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  int a;
  Rcpp::IntegerVector counts(N);
  double alpha = 1.0;

  // Get counts for the number assigned each class
  for (int i = 0; i < N; i++) {
    counts[z[i]] += 1;
  }

  for (int i = 0; i < ix.size(); i++) {
    a = ix[i];

    counts[z[a]] -= 1;
    Rcpp::NumericVector y(K);
    // Compute p(z[a] = k | all else) for each k
    for (int k = 0; k < K; k++) {
      z[a] = k;
      y[k] = RemLogLikelihoodActorPc(a,beta,z,s,K);// + log(counts[k] + alpha) - log(N + alpha);
    }
    llks.push_back(y);
    
    // Sample z[a]
    z[a] = rcategorical(y);
    counts[z[a]] += 1;
  }
  return Rcpp::List::create(Rcpp::Named("llks") = llks,Rcpp::Named("z") = z);
}



/*
 Use the precomputed statistics s to compute the likelihood
 mx: vector of events to compute by default be a vector of 0 to M-1 (inclusive)
 */

Rcpp::NumericVector RemLogLikelihoodPcSubset(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K, Rcpp::IntegerVector mx) {
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  int P = s->P;
  int M = s->M;

  double llk = 0.0; 
  double lam;
  int i,j,r,zi,zj;
  Rcpp::NumericVector llks(M);
  Rcpp::IntegerVector sen = s->sen;
  Rcpp::IntegerVector rec = s->rec;

  for (int v = 0; v < mx.size(); v++) {
    int m = mx[v];
    i = sen[m];
    j = rec[m];
    zi = z[i];
    zj = z[j];
    llk = 0;

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int r = 0; r < N; r++) {
      double lambda;
      int zr = z[r];
      if (r != i && r != j) {
        lambda  = LogLambdaPc(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
        llk  = llk - (s->times[m] - s->get_tau(m,i,r)) * exp(lambda);
        lambda  = LogLambdaPc(r,i,zr,zi,s->get_s(m,r,i),beta,N,K,P);
        llk  = llk - (s->times[m] - s->get_tau(m,r,i)) * exp(lambda);
        lambda  = LogLambdaPc(j,r,zj,zr,s->get_s(m,j,r),beta,N,K,P);
        llk  = llk - (s->times[m] - s->get_tau(m,j,r)) * exp(lambda);
        lambda  = LogLambdaPc(r,j,zr,zj,s->get_s(m,r,j),beta,N,K,P);
        llk  = llk - (s->times[m] - s->get_tau(m,r,j)) * exp(lambda);
      }
    }

    lam  = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
    llk -= (s->times[m] - s->get_tau(m,i,j)) * exp(lam);
    lam  = LogLambdaPc(j,i,zj,zi,s->get_s(m,j,i),beta,N,K,P);
    llk -= (s->times[m] - s->get_tau(m,j,i)) * exp(lam);
    
    // observed event
    llk += LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
    llks[m] = llk;
  }

  // All intensities assumed to change at the last event
  int m = M-1;
  i = sen[m];
  j = rec[m];
  zi = z[i];
  zj = z[j];
  llk = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
  for (int i = 0; i < N; i++) {
   for (int j = 0; j < N; j++) {
     zi = z[i];
     zj = z[j];
     if (i != j) {
       lam  = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
       llk -= (s->times[m] - s->get_tau(m,i,j)) * exp(lam);
     }
   }
  }
  llks[M-1] = llk;
  return llks;
}

Rcpp::NumericVector RemLogLikelihoodPc(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K) {
  RemStat *s = XPtr<RemStat>(statptr_);
  int M = s->M;
  Rcpp::IntegerVector mx = Rcpp::seq(0,M-1);
  return RemLogLikelihoodPcSubset(beta,z,statptr_,K,mx);
}

// Gradient for beta_.,k,l
// mx: should index the events where either z_i or z_j is k or l.

Rcpp::NumericVector RemGradientPcSubset(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K, Rcpp::IntegerVector mx) {
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  int P = s->P;
  int M = s->M;

  double llk = 0.0; 
  double lam;
  int i,j,r,zi,zj;
  Rcpp::NumericVector grad(P);
  Rcpp::NumericVector llks(M);
  Rcpp::IntegerVector sen = s->sen;
  Rcpp::IntegerVector rec = s->rec;

  for (int v = 0; v < mx.size(); v++) {
    int m = mx[v];
    i = sen[m];
    j = rec[m];
    zi = z[i];
    zj = z[j];

    for (int p = 0; p < 7; p++) {
      grad[p] += s->get_s(m,i,j)[p];
    }
    for (int p = 7; p < 12; p++) {
      grad[p] -= log((s->get_s(m,i,j)[p] + 1.0) / (s->get_s(m,i,j)[12] + N*(N-1)));
    }
    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    double a,b,c,d;
    for (int r = 0; r < N; r++) {
      double lambda;
      int zr = z[r];
      if (r != i && r != j) {
        lambda  = LogLambdaPc(i,r,zi,zr,s->get_s(m,i,r),beta,N,K,P);
        a = (s->times[m] - s->get_tau(m,i,r)) * exp(lambda);
        lambda  = LogLambdaPc(r,i,zr,zi,s->get_s(m,r,i),beta,N,K,P);
        b  = (s->times[m] - s->get_tau(m,r,i)) * exp(lambda);
        lambda  = LogLambdaPc(j,r,zj,zr,s->get_s(m,j,r),beta,N,K,P);
        c  = (s->times[m] - s->get_tau(m,j,r)) * exp(lambda);
        lambda  = LogLambdaPc(r,j,zr,zj,s->get_s(m,r,j),beta,N,K,P);
        d  = (s->times[m] - s->get_tau(m,r,j)) * exp(lambda);
      }
      for (int p = 0; p < 7; p++) {
        grad[p] -= a * s->get_s(m,i,r)[p];
        grad[p] -= b * s->get_s(m,r,i)[p];
        grad[p] -= c * s->get_s(m,j,r)[p];
        grad[p] -= d * s->get_s(m,r,j)[p];
      }
      for (int p = 7; p < 12; p++) {
        grad[p] -= a * log((s->get_s(m,i,r)[p] + 1.0) / (s->get_s(m,i,r)[12] + N*(N-1)));
        grad[p] -= b * log((s->get_s(m,r,i)[p] + 1.0) / (s->get_s(m,r,i)[12] + N*(N-1)));
        grad[p] -= c * log((s->get_s(m,j,r)[p] + 1.0) / (s->get_s(m,j,r)[12] + N*(N-1)));
        grad[p] -= d * log((s->get_s(m,r,j)[p] + 1.0) / (s->get_s(m,r,j)[12] + N*(N-1)));
      }
    }

    lam  = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
    a    = (s->times[m] - s->get_tau(m,i,j)) * exp(lam);
    lam  = LogLambdaPc(j,i,zj,zi,s->get_s(m,j,i),beta,N,K,P);
    b    = (s->times[m] - s->get_tau(m,j,i)) * exp(lam);
    for (int p = 0; p < 7; p++) {
      grad[p] -= a * s->get_s(m,i,j)[p];
      grad[p] -= b * s->get_s(m,j,i)[p];
    }
    for (int p = 7; p < 12; p++) {
      //      Rprintf("%f %f %i\n",smij,smji,s->get_s(m,i,j)[p]);
      grad[p] -= a * log((s->get_s(m,i,j)[p] + 1.0) / (s->get_s(m,i,j)[12] + N*(N-1)));
      grad[p] -= b * log((s->get_s(m,j,i)[p] + 1.0) / (s->get_s(m,j,i)[12] + N*(N-1)));
    }
  }

  return grad;
}


// Create numeric array with dimensions P x N x N
Rcpp::NumericVector InitializeStatisticsArray(int N, int P) {
  Rcpp::NumericVector s = Rcpp::NumericVector(Dimension(P,N,N));
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     s[0,i,j] = 1;
  //   }
  // }
  return s;
}

// Update statistics array s with an event (m,a,b).
Rcpp::NumericVector UpdateStatisticsArray(Rcpp::NumericVector s, int m, int a, int b, int N, int P) {
  // Create vector of indicators for each pshift.
  // i.e. If I[(i,j) is ab-ba from last event (a,b)]
    // Iterate through dyads having either a or b as a sender or receiver
   for (int r = 0; r < N; r++) {
     Rcpp::IntegerVector sen = Rcpp::IntegerVector::create(a,r,b,r);
     Rcpp::IntegerVector rec = Rcpp::IntegerVector::create(r,a,r,b);
      for (int k = 0; k < sen.size(); k++) {
        int i = sen[k];
        int j = rec[k];
        if (!(i==a && j==b) & i!=j) {  // deal with (a,b) at the end
          // intercept
          s[threeDIndex(0,i,j,P,N,N)] = 1;

          // P-shifts
          s[threeDIndex(1,i,j,P,N,N)] = (i!=a & i==b & j==a & j!=b);
          s[threeDIndex(2,i,j,P,N,N)] = (i!=a & i==b & j!=a & j!=b);
          s[threeDIndex(3,i,j,P,N,N)] = (i!=a & i!=b & j==a & j!=b);
          s[threeDIndex(4,i,j,P,N,N)] = (i!=a & i!=b & j!=a & j==b);
          s[threeDIndex(5,i,j,P,N,N)] = (i==a & i!=b & j!=a & j!=b);
          s[threeDIndex(6,i,j,P,N,N)] = (i==a & i!=b & j!=a & j==b);

          // Degree effects
          if (a==i && b!=j) {
            s[threeDIndex(7,i,j,P,N,N)] += 1;  // sender out degree
          }
          if (b==i && a!=j) {
            s[threeDIndex(9,i,j,P,N,N)] += 1;  // sender in degree
          }
          if (a==j && b!=i) { 
            s[threeDIndex(8,i,j,P,N,N)] += 1;  // receiver out degree
          }
          if (b==j && a!=i) {
            s[threeDIndex(10,i,j,P,N,N)] += 1; // receiver in degree
          }
          if (a==i && b==j) {
            s[threeDIndex(11,i,j,P,N,N)] += 1; // dyad count
          }
          s[threeDIndex(12,i,j,P,N,N)] = m; // changepoint count
        }
      }
    }

   // Fixed bug that didn't update stats 0 through 6
   s[threeDIndex(0,a,b,P,N,N)] = 1;

   s[threeDIndex(1,a,b,P,N,N)] = 0;
   s[threeDIndex(2,a,b,P,N,N)] = 0;
   s[threeDIndex(3,a,b,P,N,N)] = 0;
   s[threeDIndex(4,a,b,P,N,N)] = 0;
   s[threeDIndex(5,a,b,P,N,N)] = 0;
   s[threeDIndex(6,a,b,P,N,N)] = 1;

   s[threeDIndex(7,a,b,P,N,N)] += 1; // sender out degree
   s[threeDIndex(8,b,a,P,N,N)] += 1; // 
   s[threeDIndex(9,b,a,P,N,N)] += 1;
   s[threeDIndex(10,a,b,P,N,N)] += 1;
   s[threeDIndex(11,a,b,P,N,N)] += 1;
   s[threeDIndex(12,a,b,P,N,N)] = m;
  return s;
}


// Return lambda for a particular dyad (i,j) given a covariate array s and parameter array beta.
// Compute log lambda
double LogLambda(int i, int j, int zi, int zj, Rcpp::NumericVector s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = beta[threeDIndex(0,zi,zj,P,K,K)];//s[threeDIndex(0,i,j,P,N,N)];
  for (int p = 1; p < 7; p++) {
    lam += s[threeDIndex(p,i,j,P,N,N)] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  double numEvents = double(s[threeDIndex(12,i,j,P,N,N)]);//s[12];
  for (int p = 7; p < 12; p++) {
    //    lam += s[threeDIndex(p,i,j,P,N,N)]/(numEvents+1) * beta[threeDIndex(p,zi,zj,P,K,K)];
    lam += log((s[threeDIndex(p,i,j,P,N,N)] + 1)/(numEvents + N*(N-1))) * beta[threeDIndex(p,zi,zj,P,K,K)];
  }

  return lam;
}


Rcpp::NumericVector RemLogLikelihood(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector z, int N, int M,int K, int P) {


  double llk = 0.0; 
  double lam = 0;
  Rcpp::NumericVector llks(M);
  llks[0] = llk;

  int i,j,r;

  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N);
  Rcpp::NumericVector s  = InitializeStatisticsArray(N,P);

  for (int m = 0; m < M; m++) {
    i = sen[m];
    j = rec[m];
    int zi = z[i];
    int zj = z[j];
    llk = LogLambda(i,j,zi,zj,s,beta,N,K,P);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int r = 0; r < N; r++) {
      int zr = z[r];
      if (r != i) {
        lam  = LogLambda(i,r,zi,zr,s,beta,N,K,P);
        llk -= (times[m] - times[mp(i,r)]) * exp(lam);
        lam  = LogLambda(r,i,zr,zi,s,beta,N,K,P);
        llk -= (times[m] - times[mp(r,i)]) * exp(lam);
        mp(i,r) = m;
        mp(r,i) = m;
      }
      if (r != j) {
        lam  = LogLambda(j,r,zj,zr,s,beta,N,K,P);
        llk -= (times[m] - times[mp(j,r)]) * exp(lam);
        lam  = LogLambda(r,j,zr,zj,s,beta,N,K,P);
        llk -= (times[m] - times[mp(r,j)]) * exp(lam);
        mp(j,r) = m;  // update mp
        mp(r,j) = m;
      }
    }
    s = UpdateStatisticsArray(s,m,sen[m],rec[m],N,P);
    llks[m] = llk;
  }
  return llks;
}

// Compute (M,N,N) array of log rates (see lrm()).  This one uses LogLambdaPc.
Rcpp::NumericVector LogIntensityArrayPc(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K) {
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  int P = s->P;
  int M = s->M;
  Rcpp::NumericVector lrmat = Rcpp::NumericVector(Dimension(M,N,N));
  for (int m = 0; m < M; m++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int zi = z[i];
        int zj = z[j];
       double lam  = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
       lrmat[threeDIndex(m,i,j,M,N,N)] = lam;
      }
    }
  }
  return lrmat;
}

Rcpp::NumericVector LogIntensityArrayPcSubset(Rcpp::NumericVector beta, Rcpp::IntegerVector z, SEXP statptr_, int K, Rcpp::IntegerVector ix) {
  RemStat *s = XPtr<RemStat>(statptr_);
  int N = s->N;
  int P = s->P;
  int M = s->M;
  Rcpp::NumericVector lrmat = Rcpp::NumericVector(Dimension(ix.size(),N,N));
  for (int k = 0; k < ix.size(); k++) {
    int m = ix[k];
    //    Rprintf("%i %i\n",k,m);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int zi = z[i];
        int zj = z[j];
        double lam  = LogLambdaPc(i,j,zi,zj,s->get_s(m,i,j),beta,N,K,P);
        lrmat[threeDIndex(k,i,j,ix.size(),N,N)] = lam;
      }
    }
  }
  return lrmat;
}


// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).

Rcpp::NumericVector LogIntensityArray(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec,Rcpp::IntegerVector z, int N, int M,int K, int P){

  Rcpp::NumericVector lrmat = Rcpp::NumericVector(Dimension(M,N,N));
  Rcpp::NumericVector s = InitializeStatisticsArray(N,P);
  //int a,b,u,v,i,j;
  for (int m = 0; m < M; m++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int zi = z[i];
        int zj = z[j];
        lrmat[threeDIndex(m,i,j,M,N,N)] = LogLambda(i,j,zi,zj,s,beta,N,K,P);
      }
    }
    s = UpdateStatisticsArray(s,m,sen[m],rec[m],N,P);
  }
  return lrmat;
}

// Use lrm to compute the likelihood in the naive way
double RemLogLikelihoodFromArray(Rcpp::NumericVector lrm, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int N, int M){
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

Rcpp::NumericVector RemLogLikelihoodVecFromArray(Rcpp::NumericVector lrm, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int N, int M){
  double llk = 0;
  double den;
  double delta = 0;
  int i,j;
  llk = lrm[threeDIndex(0,sen[0],rec[0],M,N,N)];
  Rcpp::NumericVector llks(M);
  llks[0] = llk;
  for (int m = 1; m < M; m++) {
    den = 0.0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (i != j) {
          den += exp(lrm[threeDIndex(m,i,j,M,N,N)]);
        }
      }
    }
    llk = lrm[threeDIndex(m,sen[m],rec[m],M,N,N)];
    delta = times[m] - times[m-1];
    llk -= delta * den;
    llks[m] = llk;
  }
  return llks;
}



RCPP_MODULE(brem){
  
  class_<RemStat>( "RemStat" )
    .constructor<Rcpp::NumericVector,Rcpp::IntegerVector,Rcpp::IntegerVector,
                 int,int,int>()
    .method( "precompute", &RemStat::precompute,
             "Precompute the data structure of REM statistics")
    .method( "get_s", &RemStat::get_s,
             "Retrieve the statistics vector prior to event m for dyad (i,j)")
    .method( "get_tau", &RemStat::get_tau,
             "Retrieve the last time that lambda_(i,j) changed")
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

  // API using precomputed statistics via a RemStats object
  function( "LogLambdaPc", &LogLambdaPc);
  function( "RemLogLikelihoodPc", &RemLogLikelihoodPc ) ;
  function( "RemLogLikelihoodPcSubset", &RemLogLikelihoodPcSubset ) ;
  function( "RemGradientPcSubset", &RemGradientPcSubset ) ;
  //  function( "RemLogLikelihoodActorPc", &RemLogLikelihoodActorPc ) ;
   function( "RemGibbsPc", &RemGibbsPc ) ;
  // API using full array for statistics
  function( "InitializeStatisticsArray", &InitializeStatisticsArray ) ;
  function( "UpdateStatisticsArray", &UpdateStatisticsArray ) ;
  function( "LogLambda", &LogLambda);
  function( "RemLogLikelihood", &RemLogLikelihood ) ;
  //  function( "MultLogLikelihood", &MultLogLikelihood ) ;
  // // Create or use full log intensity array
  function( "LogIntensityArray", &LogIntensityArray ) ;
  function( "LogIntensityArrayPc", &LogIntensityArrayPc ) ;
  function( "LogIntensityArrayPcSubset", &LogIntensityArrayPcSubset ) ;
  function( "RemLogLikelihoodFromArray", &RemLogLikelihoodFromArray ) ;
  function( "RemLogLikelihoodVecFromArray", &RemLogLikelihoodVecFromArray ) ;
  //  function( "MultLogLikelihoodFromArray", &RemLogLikelihoodFromArray ) ;
  //  function( "MultLogLikelihoodVecFromArray", &RemLogLikelihoodVecFromArray ) ;
  // Utilities
  function( "rcategorical", &rcategorical);
}
