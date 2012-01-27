library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}
// http://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
Rcpp::IntegerVector sampleWithoutReplacement
(
  Rcpp::IntegerVector y,    // size of set sampling from
  int sampleSize
  )
{
  RNGScope scope;
  // Use Knuths variable names
  int& n = sampleSize;
  int N = y.size();
  Rcpp::IntegerVector samples(n);
  int t = 0; // total input records dealt with
  int m = 0; // number of items selected so far
  while (m < n) {
    Rcpp::NumericVector u = runif(1);
    if ( (N - t)*u[0] >= n - m ) {
      t++;
    } else {
      samples[m] = y[t];
      t++; m++;
    }
  }
  return samples;
}
double computeLambda(int i, int j, int a, int b, Rcpp::NumericVector beta,Rcpp::IntegerVector ps) {
  double lam = 0;
  //for (int p = 0; p < ps.size(); p++) {
    //ps[p] = 1;
    if (ps[0] == 1) {
      lam += beta[0];
    }
    if (ps[1] == 1 & i!=a & i==b & j==a & j!=b) { // ab-ba
      lam += beta[1];
    }
    if (ps[2] == 1 & i!=a & i==b & j!=a & j!=b) { // ab-by
      lam += beta[2];
    }
    if (ps[3] == 1 & i!=a & i!=b & j==a & j!=b) { // ab-xa
      lam += beta[3];
    }
    if (ps[4] == 1 & i!=a & i!=b & j!=a & j==b) { // ab-xb
      lam += beta[4];
    }
    if (ps[5] == 1 & i==a & i!=b & j!=a & j!=b) { // ab-ay
      lam += beta[5];
    }
  // if (ps[6] == 6 & i==a & i!=b & j!=a & j==b) { // ab-ab
  //    lam += beta[6];
  //  }
  //}
  return lam;
}
// All senders, receivers (ix,jx) must be 0-indexed.
// Current "weirdness": Assumes all events "occur" at time 0.
double llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M){
  // last event id that lam_ij changed
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N); 

  double llk = 0;
  int a,b,i,j,r;
  double lam;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,px);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int v = 0; v < jx.size(); v++) {
      r = jx[v];
      // Sender/receiver of last event involving i or r
      a = sen[mp(i,r)];
      b = rec[mp(i,r)];
      lam = computeLambda(i,r,a,b,beta,px);
      llk -= (times[m] - times[mp(i,r)]) * exp(lam);
      lam = computeLambda(r,i,a,b,beta,px);
      llk -= (times[m] - times[mp(r,i)]) * exp(lam);
      mp(i,r) = m;  // update mp
      mp(r,i) = m;
    }
    for (int v = 0; v < ix.size(); v++) {
      r = ix[v];
      a = sen[mp(j,r)];
      b = rec[mp(j,r)];
      lam = computeLambda(j,r,a,b,beta,px);
      llk -= (times[m] - times[mp(j,r)]) * exp(lam);
      lam = computeLambda(r,j,a,b,beta,px);
      llk -= (times[m] - times[mp(r,j)]) * exp(lam);
      mp(j,r) = m;  // update mp
      mp(r,j) = m; 
    }
  }
  return llk;
}

// Approximate likelihood
// S: sample size (where S<N)
// TODO: Consider the case when ix and jx are different sizes and the ramifications of choosing S.
double allk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M, int S){
  // last event id that lam_ij changed
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N); 

  double llk,lam,den = 0;
  int a,b,i,j,r;
  Rcpp::IntegerVector ixt, jxt; // Subsampled indices

  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,px);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    den = 0;
    jxt = sampleWithoutReplacement(jx,S);
    for (int v = 0; v < jxt.size(); v++) {
      r = jxt[v];
      // Sender/receiver of last event involving i or r
      a = sen[mp(i,r)];
      b = rec[mp(i,r)];
      lam = computeLambda(i,r,a,b,beta,px);
      den += (times[m] - times[mp(i,r)]) * exp(lam);
      lam = computeLambda(r,i,a,b,beta,px);
      den += (times[m] - times[mp(r,i)]) * exp(lam);
      mp(i,r) = m;  // update mp
      mp(r,i) = m;
    }
    ixt = sampleWithoutReplacement(ix,S);
    for (int v = 0; v < ixt.size(); v++) {
      r = ixt[v];
      a = sen[mp(j,r)];
      b = rec[mp(j,r)];
      lam = computeLambda(j,r,a,b,beta,px);
      den += (times[m] - times[mp(j,r)]) * exp(lam);
      lam = computeLambda(r,j,a,b,beta,px);
      den += (times[m] - times[mp(r,j)]) * exp(lam);
      mp(j,r) = m;  // update mp
      mp(r,j) = m; 
    }
    llk -= (N/S) * den;
  }
  return llk;
}


// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).
Rcpp::NumericVector lrm(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M){

  Rcpp::NumericMatrix mp = Rcpp::NumericMatrix(N,N); 
  Rcpp::NumericVector lrm = Rcpp::NumericVector(Dimension(M,N,N));

  int a,b,u,v,i,j;
  for (int m = 1; m < M; m++) {
    for (int u = 0; u < ix.size(); u++) {
      i = ix[u];
      for (int v = 0; v < jx.size(); v++) {
        j = jx[v];
        // get dyad that occurred last time either i or j was involved
        a = sen[mp(i,j)];
        b = rec[mp(i,j)];
        // compute rate for lambda_ij(t_m)
        lrm[threeDIndex(m,i,j,M,N,N)] = computeLambda(i,j,a,b,beta,px);
        if (i==j) {
          lrm[threeDIndex(m,i,j,M,N,N)] = -15;
        }
      }
    }
    // save which lambdas changed because of (i,j) occurring at time m
    for (int v = 0; v < N; v++) {
      mp(sen[m],v) = m;
      mp(v,sen[m]) = m;
      mp(rec[m],v) = m;
      mp(v,rec[m]) = m;
    }
  }
  return lrm;
}
RCPP_MODULE(drem){
  function( "llk", &llk ) ;
  function( "allk", &allk ) ;
  function( "lrm", &lrm ) ;
  function( "computeLambda", &computeLambda);
  function( "sampleWithoutReplacement", &sampleWithoutReplacement);
}
', plugin="Rcpp")

drem <- Module("drem",getDynLib(fx))
