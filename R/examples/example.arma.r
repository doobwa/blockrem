library(inline)
library(Rcpp)
library(RcppArmadillo)

# settings <- getPlugin("Rcpp")
# settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
# settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
code <- '
  arma::mat coeff = Rcpp::as<arma::mat>(a);
  arma::mat errors = Rcpp::as<arma::mat>(e);
  int m = errors.n_rows; int n = errors.n_cols;
  arma::mat simdata(m,n);
  simdata.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    simdata.row(row) = simdata.row(row-1)*trans(coeff)+errors.row(row);
  }
  return Rcpp::wrap(simdata);
'
rcppSim <- cxxfunction(signature(a="numeric",e="numeric"),code,plugin="RcppArmadillo")

a <- matrix(c(0.5,0.1,0.1,0.5),nrow=2)
e <- matrix(rnorm(10000),ncol=2)
rcppSim(a,e)

code <- '
  int a = as<int>(a_);
  int b = as<int>(b_);
  Rcpp::IntegerVector x = Rcpp::IntegerVector(Dimension(a,b,b)); 
  return Rcpp::wrap(x);
'
makeBigArray <- cxxfunction(signature(a_="numeric",b_="numeric"),code,plugin="Rcpp")

code <- '
  int a = as<int>(a_);
  int b = as<int>(b_);
  arma::cube x(a,b,b);
  return Rcpp::wrap(x);
'
makeBigCube <- cxxfunction(signature(a_="numeric",b_="numeric"),code,plugin="RcppArmadillo")

code <- '
  int a = as<int>(a_);
  int b = as<int>(b_);
  arma::cube x(a,b,b);
  return Rcpp::wrap(x);
'
makeBigCube <- cxxfunction(signature(a_="numeric",b_="numeric"),code,plugin="RcppArmadillo")
 
fx <- cxxfunction(,"",includes=
'
using namespace arma;
double inner(int i, int j, int zi, int zj, arma::cube s, arma::cube beta) {
  arma::cube x = s.subcube(span::all,span(i,i),span(j,j));
  arma::cube y = s.subcube(span::all,span(zi,zi),span(zj,zj));
  return as_scalar(accu(x%y));
}
double innerprod(int i, int j, int zi, int zj, Rcpp::NumericVector stats, Rcpp::NumericVector params,int N, int K, int P) {
  arma::cube s(stats.begin(),P,N,N,false);
  arma::cube beta(params.begin(),P,N,N,false);
  return inner(i,j,zi,zj,s,beta);
}
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

double innerprod2(int i, int j, int zi, int zj, Rcpp::NumericVector s, Rcpp::NumericVector beta, int N, int K, int P) {
  double lam = 0;
  for (int p = 0; p < P; p++) {
    lam += s[threeDIndex(p,i,j,P,N,N)] * beta[threeDIndex(p,zi,zj,P,K,K)];
  }
  return lam;
}
RCPP_MODULE(example){
  function( "innerprod", &innerprod ) ;
  function( "innerprod2", &innerprod2 ) ;
}
', plugin="RcppArmadillo")

example <- Module("example",getDynLib(fx))
a <- array(1:27,c(3,3,3))
b <- array(1:27,c(3,3,3))
system.time(for (i in 1:10000) example$innerprod(1,2,1,2,a,b,3,3,3))
system.time(for (i in 1:10000) example$innerprod2(1,2,1,2,a,b,3,3,3))
system.time(for (i in 1:10000) sum(a[,2,3] * b[,2,3]))
# example$fn1(10,10)
# a <- example$fn1(100,100)
# #example$fn1(1000,1000)
# 
# example$fn2(10,10)
# b <- example$fn2(100,100)  
# #example$fn2(1000,1000) 
# a
# b

system.time(example$fn2(1,100000000))
system.time(example$fn3(1,100000000))

# seg fault
#example$fn2(1000,10000) # core dump
