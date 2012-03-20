#include <iostream>
#include <algorithm>
#include <vector>
#include <Rcpp.h>

/*
 Utilities
*/

// Return the (j,k,l) element of a (J,K,L) array represented as a vector.
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

bool IsFiniteNumber(double x) {
  return (x <= DBL_MAX && x >= -DBL_MAX); 
} 

// For a vector of log probabilities, return a draw from the categorical distribution
int rcategorical (Rcpp::NumericVector lp) {
  int max = 0;
  int min = 0;
  int K = lp.size();
  for (int i = 0; i < K; i++) {
    if (lp[i] > lp[max]) {
      max = i;
    }
    if (lp[i] < lp[min]) {
      min = i;
    }
  }
  
  // If all small/large values, subtract max/min value for numeric stability.
  if (lp[max] < 0) {
    double lmax = lp[max];
    for (int i = 0; i < K; i++) {
      lp[i] -= lmax;
    }
  }
  if (lp[min] > 0) {
    double lmin = lp[min];
    for (int i = 0; i < K; i++) {
      lp[i] -= lmin;
    }
  }

  // Inverse CDF method
  Rcpp::NumericVector p = exp(lp);

  // Get sum.  If any Infinite values, return that index. 
  // TODO: This probably introduces a bias towards smaller k.
  int k;
  double cuml = 0;
  double sum = 0;
  for (k=0;k<K;k++) {
    if (!IsFiniteNumber(p[k])) {
      return k;
    }
    sum += p[k];
  }

  double r = Rcpp::as<double>(Rcpp::runif(1));
  for (k=0;k<K;k++) {
    cuml += p[k]/sum;
    if (r < cuml) {
      break;
    } 
  }
  return k;
}
