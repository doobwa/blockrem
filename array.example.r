library(inline)
library(Rcpp)
fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return j*K*L + k*L + l;
}
Rcpp::List fn(int M, int N1, int N2){
  Rcpp::NumericVector lrm = Rcpp::NumericVector( Dimension(M,N1,N1) );
  lrm[threeDIndex(1,2,3,M,N1,N2)] = 3;
  return Rcpp::List::create(Rcpp::Named( "lrm" ) = lrm);
}
RCPP_MODULE(foo){
  function( "fn", &fn ) ;
}
', plugin="Rcpp")

foo <- Module("foo",getDynLib(fx))
foo$fn(3,2,1)
