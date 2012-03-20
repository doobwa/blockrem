#include <iostream>
#include <algorithm>
#include <vector>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

class RemDP {
public:

  RemDP(Rcpp::List edgelist_,  int N_) : 
    edgelist(edgelist_), N(N_) {
    times = edgelist(0);
    sen = edgelist(1);
    rec = edgelist(2);
    M = sen.size();
    P = 13;
    initialize();
    //precompute();
  }

  // DEPRECATED
  RemDP(Rcpp::NumericVector times_, 
        Rcpp::IntegerVector sen_,  
        Rcpp::IntegerVector rec_, 
        int N_, int M_, int P_) : 
    times(times_),sen(sen_),rec(rec_),
    N(N_),M(M_),P(P_) {
    initialize();
    //precompute();
  }

  void initialize() {

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

  // Return pointer
  SEXP get_ptr() {
    return Rcpp::XPtr<RemDP>(this, false);
  }

  double LogLambda(int m, int i, int j) {
    vector< int > s = this->get_s(m,i,j);
    double lam = eta(z[i],z[j]); // intercept
    for (int p = 1; p < 7; p++) {
      lam += s[p] * (beta(z[i],p) + gamma(z[j],p));
    }
    double numEvents = double(s[12]);
    for (int p = 7; p < 12; p++) {
      lam += log((s[p]+1)/(numEvents + N*(N-1))) *
        (beta(z[i],p) + gamma(z[j],p));
    }
    return lam;
  }

  double LogLikelihoodByEvent(Rcpp::IntegerVector ix) {
    return 2.0;
  }

  double LogLikelihoodByNode() {

  }

  double LogLikelihood() {
    return 1.0;
  }

  Rcpp::NumericVector LogIntensityArraySubset(Rcpp::IntegerVector ix) {

  }

  // Return gradient with respect to a subset of variables and events ix
  Rcpp::List Gradient(int vx, Rcpp::IntegerVector ix) {

  }


  // Saves probability of each assignment
  void Gibbs() {

  }

  // Make sure several aspects of the model are in order
  bool Check() {

    // beta, gamma, and eta have same number of rows (K)
    int k1 = eta.nrow();
    int k2 = beta.nrow();
    int k3 = gamma.nrow();
    
    // all z values are less than K
    int zmax = 0;
    for (int i = 0; i < z.size(); i++) {
      if (z[i] > zmax) {
        zmax = z[i];
      }
    }

    Rprintf("%i %i %i %i\n",k1,k2,k3,zmax);

    return (k1 == k2 && k1 == k3 && (zmax<k1));
  }

  Rcpp::DataFrame get_edgelist() {
    return edgelist;
  }

  Rcpp::List get_params() {
    Rcpp::List ret; 
    ret["eta"]   = eta; 
    ret["beta"]  = beta;
    ret["gamma"] = gamma;
    ret["z"] = z;
    return ret;
  }

  void set_params(Rcpp::List params) {
    eta   = Rcpp::as<Rcpp::NumericMatrix>(params["eta"]);
    beta  = Rcpp::as<Rcpp::NumericMatrix>(params["beta"]);
    gamma = Rcpp::as<Rcpp::NumericMatrix>(params["gamma"]);
    z     = Rcpp::as<Rcpp::IntegerVector>(params["z"]);
  }

  // Number of nodes, events, parameters, and clusters.
  int N, M, P;
  Rcpp::DataFrame edgelist;
  Rcpp::NumericVector times;
  Rcpp::IntegerVector sen;
  Rcpp::IntegerVector rec;

private:

  // Variables relating to the statistics s:

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

  // Model parameters:

  Rcpp::NumericMatrix eta;
  Rcpp::NumericMatrix beta;
  Rcpp::NumericMatrix gamma;
  Rcpp::IntegerVector z;

};


RCPP_MODULE(remdp){
  
  class_<RemDP>( "RemDP" )
    .constructor<Rcpp::DataFrame,int>()
    .constructor<Rcpp::NumericVector,Rcpp::IntegerVector,Rcpp::IntegerVector,
                 int,int,int>()
    .method( "precompute", &RemDP::precompute,
             "Precompute the data structure of REM statistics.")
    .method( "get_s", &RemDP::get_s,
             "Retrieve statistics vector prior to event m for dyad (i,j)")
    .method( "get_tau", &RemDP::get_tau,
             "Retrieve the last time that lambda_(i,j) changed")
    .method( "get_v", &RemDP::get_v,
             "vector where element k is the event index of the k'th changepoint for (i,j).  i.e. if w(i,j)[m] = k then v[i,j,k] = m")
    .method( "get_w", &RemDP::get_w,
             "vector where element m is the index of the previous changepoint for (i,j).  i.e. if w(i,j)[m] = k then v[i,j,k] = m")
    .method( "get_ptr", &RemDP::get_ptr,
             "Get pointer to this object.")
    .method( "get_edgelist", &RemDP::get_edgelist,
             "Get edgelist.")
    .method( "get_params", &RemDP::get_params,
             "Get a list of  eta, beta, gamma, and z.")
    .method( "set_params", &RemDP::set_params,
             "Set parameter values with a list of  eta, beta, gamma, and z.")
    .method( "Check", &RemDP::Check,
             "Make sure parameter matrices and latent assignments have proper dimensions")
    .method( "LogLambda", &RemDP::LogLambda,
             "Log(lambda_ij(t_m)) for a particular dyad at event m.")
    .method( "LogLikelihood", &RemDP::LogLikelihood,
             "Log likelihood on all observed events.")
    .method( "LogLikelihoodByEvent", &RemDP::LogLikelihoodByEvent,
             "Log likelihood for only a subset of events.")
    .method( "LogLikelihoodByNode", &RemDP::LogLikelihoodByNode,
             "Log likelihood for a single node.")
    .method( "LogIntensityArraySubset", &RemDP::LogIntensityArraySubset,
             "Log intensity array for only a subset of events.")
    ;

}
