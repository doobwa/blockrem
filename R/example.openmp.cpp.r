library(inline)
library(Rcpp)

settings <- getPlugin("Rcpp")
settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)

fx <- cxxfunction(,"",includes=
  '
#include <iostream>
using namespace std;
#include <omp.h>

double mysqrt1(double x, std::vector<double> s) {
  return sqrt(x);
}
double mysqrt2(double x, Rcpp::NumericVector &s) {
  return sqrt(x);
}
double fn1(int K, int N) {
  std::vector<double> x(N);
  std::vector<double> y(N);
  for (int k = 0; k < K; k++) {
   // #pragma omp parallel
   // {
   // #pragma omp for
      for (int i=0; i<N; i++) {
        y[i] = mysqrt1(x[i],x);
      }
    //}
  }
  return 0.0;
}
double fn2(int K, int N) {
  omp_set_num_threads(16);
  Rcpp::NumericVector x(N);
  Rcpp::NumericVector y(N);
  for (int k = 0; k < K; k++) {
    #pragma omp parallel
    {
    #pragma omp for
      for (int i=0; i<N; i++) {
        y(i) = mysqrt2(x(i),x);
      }
    }
  }
  return 0.0;
}
double fn3(int K, int N) {
  Rcpp::NumericVector x(N);
  Rcpp::NumericVector y(N);
  for (int k = 0; k < K; k++) {
      for (int i=0; i<N; i++) {
        y(i) = mysqrt2(x(i),x);
      }
  }
  return 0.0;
}

RCPP_MODULE(example){
  function( "fn1", &fn1 ) ;
  function( "fn2", &fn2 ) ;
  function( "fn3", &fn3 ) ;
}
', plugin="Rcpp",settings=settings)


example <- Module("example",getDynLib(fx))

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