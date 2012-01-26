library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}
Rcpp::IntegerVector testThreeDIndex() {
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(Dimension(3,4,5));
  mp[threeDIndex(1,1,1,3,4,5)] = 1;
  return mp;
}

int testIntConv(Rcpp::IntegerVector ps) {
  int r;
  if (ps[0] == 1) {
    r = 1;
  } else {
    r = 2;
  }
  int k = ps[ps[0]];
  return k;
}
double computeLambdaOld(int i, int j, int a, int b, Rcpp::NumericVector beta) {
  double lam = 0;
    lam += beta[0];
    if (i!=a & i==b & j==a & j!=b) { // ab-ba
      lam += beta[1];
    }
    if (i!=a & i==b & j!=a & j!=b) { // ab-by
      lam += beta[2];
    }
    if (i!=a & i!=b & j==a & j!=b) { // ab-xa
      lam += beta[3];
    }
    if (i!=a & i!=b & j!=a & j==b) { // ab-xb
      lam += beta[4];
    }
    if (i==a & i!=b & j!=a & j!=b) { // ab-ay
      lam += beta[5];
    }
  return lam;
}
double computeLambda(int i, int j, int a, int b, Rcpp::NumericVector beta,Rcpp::IntegerVector ps) {
  double lam = 0;
  //for (int p = 0; p < ps.size(); p++) {
    //ps[p] = 1;
    if (ps[0] == 1) {
      lam += beta[0];
    }
    if (ps[1] == 1 & i!=a & i==b & j==a & j!=b) { // ab-ba
      lam += beta[1];
    }
    if (ps[2] == 1 & i!=a & i==b & j!=a & j!=b) { // ab-by
      lam += beta[2];
    }
    if (ps[3] == 1 & i!=a & i!=b & j==a & j!=b) { // ab-xa
      lam += beta[3];
    }
    if (ps[4] == 1 & i!=a & i!=b & j!=a & j==b) { // ab-xb
      lam += beta[4];
    }
    if (ps[5] == 1 & i==a & i!=b & j!=a & j!=b) { // ab-ay
      lam += beta[5];
    }
  // if (ps[6] == 6 & i==a & i!=b & j!=a & j==b) { // ab-ab
  //    lam += beta[6];
  //  }
  //}
  return lam;
}
// All senders, receivers must be 0-indexed.
// Current "weirdness": Assumes all events "occur" at time 0.
Rcpp::List llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M, int P){
  // last event id that lam_ij changed
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N); 

  double llk = 0;
  int a,b,i,j,r;
  double lam;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,px);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int v = 0; v < jx.size(); v++) {
      r = jx[v];
      // Sender/receiver of last event involving i or r
      a = sen[mp(i,r)];
      b = rec[mp(i,r)];
      lam = computeLambda(i,r,a,b,beta,px);
      llk -= (times[m] - times[mp(i,r)]) * exp(lam);
      lam = computeLambda(r,i,a,b,beta,px);
      llk -= (times[m] - times[mp(r,i)]) * exp(lam);
      mp(i,r) = m;  // update mp
      mp(r,i) = m;
    }
    for (int v = 0; v < ix.size(); v++) {
      r = ix[v];
      a = sen[mp(j,r)];
      b = rec[mp(j,r)];
      lam = computeLambda(j,r,a,b,beta,px);
      llk -= (times[m] - times[mp(j,r)]) * exp(lam);
      lam = computeLambda(r,j,a,b,beta,px);
      llk -= (times[m] - times[mp(r,j)]) * exp(lam);
      mp(j,r) = m;  // update mp
      mp(r,j) = m; 
    }
  }
  return Rcpp::List::create(Rcpp::Named( "llk" ) = llk);
}

// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).
Rcpp::NumericVector lrm(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M, int P){

  Rcpp::NumericMatrix mp = Rcpp::NumericMatrix(N,N); 
  Rcpp::NumericVector lrm = Rcpp::NumericVector(Dimension(M,N,N));

  int a,b;
  for (int m = 1; m<M; m++) {
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        // get dyad that occurred last time either i or j was involved
        a = sen[mp(i,j)];
        b = rec[mp(i,j)];
        // compute rate for lambda_ij(t_m)
        lrm[threeDIndex(m,i,j,M,N,N)] = computeLambda(i,j,a,b,beta,px);
      }
    }
    // save which lambdas changed because of (i,j) occurring at time m
    for (int v = 0; v < N; v++) {
      mp(sen[m],v) = m;
      mp(v,sen[m]) = m;
      mp(rec[m],v) = m;
      mp(v,rec[m]) = m;
    }
  }
  return lrm;
}
RCPP_MODULE(drem){
  function( "llk", &llk ) ;
  function( "lrm", &lrm ) ;
  function( "testIntConv", &testIntConv );
  function( "computeLambda", &computeLambda);
  function( "computeLambdaOld", &computeLambdaOld ) ;
  function( "testThreeDIndex", &testThreeDIndex ) ;
}
', plugin="Rcpp")

drem <- Module("drem",getDynLib(fx))
