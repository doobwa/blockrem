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

double mysqrt1(double x, Rcpp::NumericVector s) {
  return sqrt(x);
}
double mysqrt2(double x, Rcpp::NumericVector &s) {
  return sqrt(x);
}
double fn1(int n, Rcpp::NumericVector y) {
int   i;
float a[n], b[n], sum; 

/* Some initializations */
for (i=0; i < n; i++)
  a[i] = b[i] = i * 1.0;
sum = 0.0;

  for (i=0; i < n; i++)
    sum = sum + a[i];
  return sum;
}
double fn2(int n, Rcpp::NumericVector y) {
  omp_set_num_threads(8);
int   i;
float a[n], b[n], sum; 

/* Some initializations */
for (i=0; i < n; i++)
  a[i] = b[i] = i * 1.0;
sum = 0.0;

#pragma omp parallel for reduction(+:sum)
  for (i=0; i < n; i++)
    sum = sum + a[i];

  return sum;
}

RCPP_MODULE(example){
  function( "fn1", &fn1 ) ;
  function( "fn2", &fn2 ) ;
}
', plugin="Rcpp",settings=settings)


example <- Module("example",getDynLib(fx))
example$fn2(10,rnorm(5))
system.time(example$fn1(1000000,1:1000))
system.time(example$fn2(100,1:1000))

# seg fault
#example$fn2(1000,10000) # core dump
