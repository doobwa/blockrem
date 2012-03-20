#include <iostream>
#include <algorithm>
#include <vector>
#include <Rcpp.h>

/*
 Utilities
*/

// Return the (j,k,l) element of a (J,K,L) array represented as a vector.
int threeDIndex(int j, int k, int l, int J, int K, int L);

bool IsFiniteNumber(double x);

// For a vector of log probabilities, return a draw from the categorical distribution
int rcategorical (Rcpp::NumericVector lp);
