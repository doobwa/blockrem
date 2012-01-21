
fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return j*K*L + k*L + l;
}
// beta currently: AB-BA,AB-BY,AB-XA,
Rcpp::List llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int M, int N1, int N2, int P){
  Rcpp::NumericVector lrm = Rcpp::NumericVector( Dimension(M,N1,N1) );
  double llk = 0;
  int a,b;
  for (int m = 1; m<M; m++) {
    a = sen[m-1];
    b = rec[m-1];
    for (int i = 0; i<N1; i++) {
      for (int j = 0; j<N2; j++) {
        // TODO: Make interface for using different effects cleaner
        lrm[threeDIndex(m,i,j,M,N1,N2)] = beta[0];

        // PSHIFTS
        if (i==b & j==a) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[1]; //ab-ba
        }
        if (i==b & j!=a) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[2]; //ab-by
        }
        if (i!=b & j==a) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[3]; //ab-xa
        }
        if (i!=a & j==b) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[4]; //ab-xb
        }
        if (i!=a & i!=b & j!=a & j!=b) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[5]; //ab-xy
        }
        if (i==b & j!=b) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[6]; //ab-ay
        }

        // Add to loglikelihood
        double lambda = lrm[threeDIndex(m,i,j,M,N1,N2)];
        llk += exp( - (times[m] - times[m-1]) * lambda);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named( "lrm" ) = lrm,
                            Rcpp::Named( "llk" ) = llk);
}
RCPP_MODULE(foo){
  function( "llk", &llk ) ;
}
', plugin="Rcpp")


foo <- Module("foo",getDynLib(fx))
M <- 10000
N <- 20
P <- 7
times <- sort(runif(M,0,100))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
beta <- rnorm(P)
system.time(llk <- foo$llk(beta,times,sen,rec,M,N,N,P))


fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return j*K*L + k*L + l;
}
// beta currently: AB-BA,AB-BY,AB-XA,
Rcpp::List llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int M, int N1, int N2, int P){
  Rcpp::NumericVector lrm = Rcpp::NumericVector( Dimension(M,N1,N1) );
  double llk = 0;
  int a,b;
  for (int m = 1; m<50; m++) {
    a = sen[m-1];
    b = rec[m-1];
    for (int i = 0; i<N1; i++) {
      for (int j=0; j<N2; j++) {
        // TODO: Make interface for using different effects cleaner
        //lrm[threeDIndex(m,i,j,M,N1,N2)] = beta[0];

        // PSHIFTS
        if (i==b & j==a) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[1]; //ab-ba
        }
        if (i==b & j!=a) {
          lrm[threeDIndex(m,i,j,M,N1,N2)] += beta[2]; //ab-by
        }

        // Add to loglikelihood
        double lambda = lrm[threeDIndex(m,i,j,M,N1,N2)];
      llk += exp( - (times[m] - times[m-1]) * lambda);
    }
    }
  }
  return Rcpp::List::create(Rcpp::Named( "lrm" ) = lrm,
                            Rcpp::Named( "llk" ) = llk);
}
RCPP_MODULE(foo){
  function( "llk", &llk ) ;
}
', plugin="Rcpp")

foo <- Module("foo",getDynLib(fx))
M <- 10000
N <- 20
P <- 7
times <- sort(runif(M,0,100))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
beta <- rnorm(P)
system.time(llk <- foo$llk(beta,times,sen,rec,M,N,N,P))

