library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return j*K*L + k*L + l;
}
double computeLambda(int i, int j, int a, int b, Rcpp::NumericVector beta) {
  double lam = beta[0];
  if (i==b & j==a) { // ab-ba
    lam += beta[1];  
  }
  if (i==b & j!=a) { // ab-by
    lam += beta[2]; 
  }
  if (i!=b & j==a) { // ab-xa
    lam += beta[3];
  }
  if (i!=a & j==b) { // ab-xb
    lam += beta[4];
  }
  if (i==b & j!=b) { // ab-ay
    lam += beta[5];
  }
  return exp(lam);
}
// All senders, receivers must be 0-indexed.
// Current "weirdness": Assumes all events "occur" at time 0.
Rcpp::List llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int M, int N1, int N2, int P){

  Rcpp::NumericMatrix mp = Rcpp::NumericMatrix(N1,N2);  // last event id that lam_ij changed

  double llk = 0;
  int a,b,i,j;
  double lam_ij, lam_ir, lam_rj;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int r = 0; r<N1; r++) {
      a = sen[mp(i,r)];
      b = rec[mp(i,r)];
      lam_ir = computeLambda(i,r,a,b,beta);
      llk -= (times[m] - times[mp(i,r)]) * lam_ir;
      mp(i,r) = m;  // update mp
    }
    for (int r = 0; r<N2; r++) {
      a = sen[mp(r,j)];
      b = rec[mp(r,j)];
      lam_rj = computeLambda(r,j,a,b,beta);
      llk -= (times[m] - times[mp(r,j)]) * lam_rj;
      mp(r,j) = m;
    }
  }
  return Rcpp::List::create(Rcpp::Named( "llk" ) = llk,
                            Rcpp::Named( "mp" )  = mp);
}
RCPP_MODULE(foo){
  function( "llk", &llk ) ;
}
', plugin="Rcpp")

foo <- Module("foo",getDynLib(fx))
M <- 10000
N <- 100
P <- 6
times <- sort(runif(M,0,100))
sen <- sample(1:N,M,replace=TRUE) - 1
rec <- sample(1:N,M,replace=TRUE) - 1
beta <- rnorm(P)
system.time(llk <- foo$llk(beta,times,sen,rec,M,N,N,P))

