library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return j*K*L + k*L + l;
}
double computeLambda(int i, int j, int a, int b, Rcpp::NumericVector beta,int x) {
  double lam = 0;
  if (x == 0 || x > 1) {
    lam += beta[0];
  }
  if (x > 0) {
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
  // if (i==a & i!=b & j!=a & j==b) { // ab-ab
  //    lam += beta[6];
  //  }
  }
  return lam;
}
// All senders, receivers must be 0-indexed.
// Current "weirdness": Assumes all events "occur" at time 0.
Rcpp::List llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx, int N, int M, int P, int q){
  // last event id that lam_ij changed
  Rcpp::NumericMatrix mp = Rcpp::NumericMatrix(N+1,N+1); // +1 hack 

  double llk = 0;
  int a,b,i,j,r;
  double lam;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,q);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int v = 0; v < jx.size(); v++) {
      r = jx[v];
      a = sen[mp(i,r)];
      b = rec[mp(i,r)];
      lam = computeLambda(i,r,a,b,beta,q);
      llk -= (times[m] - times[mp(i,r)]) * exp(lam);
      mp(i,r) = m;  // update mp
      a = sen[mp(r,i)];
      b = rec[mp(r,i)];
      lam = computeLambda(r,i,a,b,beta,q);
      llk -= (times[m] - times[mp(r,i)]) * exp(lam);
      mp(r,i) = m;  // update mp
    }
    for (int v = 0; v < ix.size(); v++) {
      r = ix[v];
      a = sen[mp(j,r)];
      b = rec[mp(j,r)];
      lam = computeLambda(j,r,a,b,beta,q);
      llk -= (times[m] - times[mp(j,r)]) * exp(lam);
      mp(j,r) = m;  // update mp
      a = sen[mp(r,j)];
      b = rec[mp(r,j)];
      lam = computeLambda(r,j,a,b,beta,q);
      llk -= (times[m] - times[mp(r,j)],q) * exp(lam);
      mp(r,j) = m;  // update mp
    }
  }
  return Rcpp::List::create(Rcpp::Named( "llk" ) = llk);
}
RCPP_MODULE(drem){
  function( "llk", &llk ) ;
}
', plugin="Rcpp")

drem <- Module("drem",getDynLib(fx))