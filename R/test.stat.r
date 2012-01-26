require(inline)
require(Rcpp)
require(RcppEigen)
require(testthat)

fx <- cxxfunction(,"",includes=
'
int llk(Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, int N1, int N2, int P ){
  int llk = 0;
  for (int m = 0; m<M; m++) {
    for (int i = 0; i<N1; i++) {
      for (int j = 0; j<N2; j++) {
        for (int p = 0; p<P; p++) {
          llk += 1;
        }
      }
    }
  }
  return llk;
}
RCPP_MODULE(foo){
  function( "llk", &llk ) ;
}
', plugin="Rcpp")

foo <- Module("foo",getDynLib(fx))
# norm <- function(M,N1,N2,P) M*N1*N2*P

system.time(foo$norm(10000,200,200,4))

transCpp <- '
using Eigen::Map;
using Eigen::MatrixXi;
// Map the integer matrix AA from R
const Map<MatrixXi> A(as<Map<MatrixXi> >(AA));
// evaluate and return the transpose of A
const MatrixXi At(A.transpose());
return wrap(At);
'
ftrans <- cxxfunction(signature(AA="matrix"), transCpp, plugin="RcppEigen") 
a <- matrix(1:6,nr=2)   
all.equal(ftrans(a),t(a))


