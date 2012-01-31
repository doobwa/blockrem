library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

// Update each s(t,i,j) vector with event (a,b).
// s: List of NxN matrices with named elements.  Each matrix represents the current value for that statistic.

Rcpp::NumericVector initializeStatistics(int N, int P) {
  Rcpp::NumericVector s     = Rcpp::NumericVector(Dimension(P,N,N));

  // Intercept statistic: all (0,i,j) are equal to 1
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

//
double computeLambda(int i, int j, int zi, int zj, Rcpp::NumericVector s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = 0;
  for (int p = 0; p < P; p++) {
    lam += s[threeDIndex(p,i,j,P,N,N)] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}


double llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M){
  // last event id that lam_ij changed
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N); 

  int a,b,i,j,r;
  double lam;

  double llk = 0.0;  // initial event assumed uniform across risk set
  Rcpp::NumericVector llks(M);
  llks[0] = llk;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    //llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,px);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int v = 0; v < jx.size(); v++) {
      r = jx[v];
      if (r != i) {
        // Sender/receiver of last event involving i or r
        a = sen[mp(i,r)];
        b = rec[mp(i,r)];
        //lam = computeLambda(i,r,a,b,beta,px);
        llk -= (times[m] - times[mp(i,r)]) * exp(lam);
        //lam = computeLambda(r,i,a,b,beta,px);
        llk -= (times[m] - times[mp(r,i)]) * exp(lam);
        mp(i,r) = m;  // update mp
        mp(r,i) = m;
      }
    }
    for (int v = 0; v < ix.size(); v++) {
      r = ix[v];
      if (r != j) {
        a = sen[mp(j,r)];
        b = rec[mp(j,r)];
        //lam = computeLambda(j,r,a,b,beta,px);
        llk -= (times[m] - times[mp(j,r)]) * exp(lam);
        //lam = computeLambda(r,j,a,b,beta,px);
        llk -= (times[m] - times[mp(r,j)]) * exp(lam);
        mp(j,r) = m;  // update mp
        mp(r,j) = m;
      }
    }
    llks[m] = llk;
  }
  return llk;
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
  function( "llk2", &llk2 ) ;
  function( "lrm", &lrm ) ;
  function( "updateStatistics", &updateStatistics);
  function( "initializeStatistics", &initializeStatistics);
  function( "computeLambda", &computeLambda);
}
', plugin="Rcpp")

brem <- Module("brem",getDynLib(fx))
