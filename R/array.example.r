library(inline)
library(Rcpp)
fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}
Rcpp::IntegerVector putval(int v, int i, int j, int k, int I, int J, int K){
  Rcpp::IntegerVector lrm = Rcpp::IntegerVector(Dimension(I,J,K));
  lrm[threeDIndex(i-1,j-1,k-1,I,J,K)] = v;
  return lrm;
}
RCPP_MODULE(foo){
  function( "putval", &putval ) ;
}
', plugin="Rcpp")

foo <- Module("foo",getDynLib(fx))

foo$putval(5,1,2,3,2,3,4)[1,2,3] # 5


m <- matrix(0,3,4)
mind <- function(j,k,J,K) (k-1)*J + j
m[mind(1,2,3,4)] <- 100
m
a <- array(0,c(2,3,4))
aind <- function(i,j,k,I,J,K) k*I*J + j*I + i
a[aind(1-1,2-1,3-1,2,3,4) + 1] <- 1
a[1,2,3]

tmp <- foo$fn(100,1,2,3,100,10,10)