library(inline)
library(Rcpp)

fx <- cxxfunction(,"",includes=
  '
int threeDIndex(int j, int k, int l, int J, int K, int L) { 
  return l*J*K + k*J + j;
}

// Update each s(t,i,j) vector with event (a,b).
// s: List of NxN matrices with named elements.  Each matrix represents the current value for that statistic.
Rcpp::List updateStatistics(Rcpp::List s, int a, int b, int N) {
  // Create vector of indicators for each pshift.
  // i.e. If I[(i,j) is ab-ba from last event (a,b)]
    Rcpp::NumericMatrix s_abba = s["abba"];
    Rcpp::NumericMatrix s_abby = s["abby"];
    Rcpp::NumericMatrix s_abxa = s["abxa"];
    Rcpp::NumericMatrix s_abxb = s["abxb"];
    Rcpp::NumericMatrix s_abay = s["abay"];
    Rcpp::NumericMatrix s_abab = s["abab"];
    Rcpp::NumericMatrix s_sod = s["sod"];
    Rcpp::NumericMatrix s_rod = s["rod"];
    Rcpp::NumericMatrix s_sid = s["sid"];
    Rcpp::NumericMatrix s_rid = s["rid"];

    // Iterate through dyads having either a or b as a sender or receiver
    for (int r = 0; r < N; r++) {
      Rcpp::IntegerVector sen = Rcpp::IntegerVector::create(a,r,b,r);
      Rcpp::IntegerVector rec = Rcpp::IntegerVector::create(r,a,r,b);
      for (int k = 0; k < sen.size(); k++) {
        int i = sen[k];
        int j = rec[k];
        if (i != j) {
          // P-shifts
          s_abba(i,j) = (i!=a & i==b & j==a & j!=b);
          s_abby(i,j) = (i!=a & i==b & j!=a & j!=b);
          s_abxa(i,j) = (i!=a & i!=b & j==a & j!=b);
          s_abxb(i,j) = (i!=a & i!=b & j!=a & j==b);
          s_abay(i,j) = (i==a & i!=b & j!=a & j!=b);
          s_abab(i,j) = (i==a & i!=b & j!=a & j==b);

          // Degree effects
          if (a==i & b!=j) {
            s_sod(i,j) += 1;
          }
          if (b==i & a!=j) {
            s_sid(i,j) += 1;
          }
          if (a==j & b!=i) { 
            s_rod(i,j) += 1;
          }
          if (b==j & a!=i) {
            s_rid(i,j) += 1;
          }
        }
      }
    }
   s_sid(b,a) += 1;
   s_rid(a,b) += 1;
   s_sod(a,b) += 1;
   s_rod(b,a) += 1;
    s["abba"] = s_abba;
    s["abby"] = s_abby;
    s["abxa"] = s_abxa;
    s["abxb"] = s_abxb;
    s["abay"] = s_abay;
    s["abab"] = s_abab;
    s["sod"] = s_sod;
    s["sid"] = s_sid;
    s["rod"] = s_rod;
    s["rid"] = s_rid;

  // Increment in/out degrees for i and j
  return s;
}

//
double computeLambda(int i, int j, int zi, int zj, Rcpp::List s, Rcpp::List beta) {
    Rcpp::NumericMatrix s_abba = s["abba"];
    Rcpp::NumericMatrix s_abby = s["abby"];
    Rcpp::NumericMatrix s_abxa = s["abxa"];
    Rcpp::NumericMatrix s_abxb = s["abxb"];
    Rcpp::NumericMatrix s_abay = s["abay"];
    Rcpp::NumericMatrix s_abab = s["abab"];
    Rcpp::NumericMatrix s_sod = s["sod"];
    Rcpp::NumericMatrix s_rod = s["rod"];
    Rcpp::NumericMatrix s_sid = s["sid"];
    Rcpp::NumericMatrix s_rid = s["rid"];

  Rcpp::NumericMatrix beta_intercept = beta["intercept"];
  Rcpp::NumericMatrix beta_abba = beta["abba"];
  Rcpp::NumericMatrix beta_abby = beta["abby"];
  Rcpp::NumericMatrix beta_abxa = beta["abxa"];
  Rcpp::NumericMatrix beta_abxb = beta["abxb"];
  Rcpp::NumericMatrix beta_abay = beta["abay"];
  Rcpp::NumericMatrix beta_abab = beta["abab"];
  Rcpp::NumericMatrix beta_sod  = beta["sod"];
  Rcpp::NumericMatrix beta_rod  = beta["rod"];
  Rcpp::NumericMatrix beta_sid  = beta["sid"];
  Rcpp::NumericMatrix beta_rid  = beta["rid"];

  double lam = beta_intercept(zi,zj) + 
              s_abba(i,j) * beta_abba(zi,zj) + 
              s_abby(i,j) * beta_abby(zi,zj) +
              s_abxa(i,j) * beta_abxa(zi,zj) +
              s_abxb(i,j) * beta_abxb(zi,zj) +
              s_abay(i,j) * beta_abay(zi,zj) +
              s_abab(i,j) * beta_abab(zi,zj) + 
              s_sod(i,j)  * beta_sod(zi,zj)  + 
              s_rod(i,j)  * beta_rod(zi,zj)  + 
              s_sid(i,j)  * beta_sid(zi,zj)  + 
              s_rid(i,j)  * beta_rid(zi,zj);

  return lam;
}
// All senders, receivers (ix,jx) must be 0-indexed.
// Current "weirdness": Assumes all events "occur" at time 0.
double llk(Rcpp::NumericVector beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec, Rcpp::IntegerVector ix, Rcpp::IntegerVector jx,Rcpp::IntegerVector px, int N, int M){
  // last event id that lam_ij changed
  Rcpp::IntegerMatrix mp = Rcpp::IntegerMatrix(N,N); 

  int a,b,i,j,r;
  double lam;

  double llk = 0.0;  // initial event assumed uniform across risk set
  Rcpp::NumericVector llks(M);
  llks[0] = llk;
  for (int m = 1; m<M; m++) {

    i = sen[m];
    j = rec[m];
    //llk += computeLambda(i,j,sen[mp(i,j)],rec[mp(i,j)],beta,px);

    // Loop through dyads (i,r) and (r,j) whose intensities change due to event m
    for (int v = 0; v < jx.size(); v++) {
      r = jx[v];
      if (r != i) {
        // Sender/receiver of last event involving i or r
        a = sen[mp(i,r)];
        b = rec[mp(i,r)];
        //lam = computeLambda(i,r,a,b,beta,px);
        llk -= (times[m] - times[mp(i,r)]) * exp(lam);
        //lam = computeLambda(r,i,a,b,beta,px);
        llk -= (times[m] - times[mp(r,i)]) * exp(lam);
        mp(i,r) = m;  // update mp
        mp(r,i) = m;
      }
    }
    for (int v = 0; v < ix.size(); v++) {
      r = ix[v];
      if (r != j) {
        a = sen[mp(j,r)];
        b = rec[mp(j,r)];
        //lam = computeLambda(j,r,a,b,beta,px);
        llk -= (times[m] - times[mp(j,r)]) * exp(lam);
        //lam = computeLambda(r,j,a,b,beta,px);
        llk -= (times[m] - times[mp(r,j)]) * exp(lam);
        mp(j,r) = m;  // update mp
        mp(r,j) = m;
      }
    }
    llks[m] = llk;
  }
  return llk;
}


// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).

Rcpp::NumericVector lrm(Rcpp::List beta, Rcpp::NumericVector times, Rcpp::IntegerVector sen, Rcpp::IntegerVector rec,Rcpp::IntegerVector z, int N, int M){

  Rcpp::NumericVector lrm = Rcpp::NumericVector(Dimension(M,N,N));
  Rcpp::List s;
  s["abba"] = Rcpp::NumericMatrix(N,N);
  s["abby"] = Rcpp::NumericMatrix(N,N);
  s["abxa"] = Rcpp::NumericMatrix(N,N);
  s["abxb"] = Rcpp::NumericMatrix(N,N);
  s["abay"] = Rcpp::NumericMatrix(N,N);
  s["abab"] = Rcpp::NumericMatrix(N,N);
  s["sod"] = Rcpp::NumericMatrix(N,N);
  s["rod"] = Rcpp::NumericMatrix(N,N);
  s["sid"] = Rcpp::NumericMatrix(N,N);
  s["rid"] = Rcpp::NumericMatrix(N,N);

  //int a,b,u,v,i,j;
  for (int m = 1; m < M; m++) {
    s = updateStatistics(s,sen[m-1],rec[m-1],N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        int zi = z[i];
        int zj = z[j];
   //     if (i==0 && j==5) {
    //      Rprintf("%i %i\\n",zi,zj);
    //    }
        lrm[threeDIndex(m,i,j,M,N,N)] = computeLambda(i,j,zi,zj,s,beta);
      }
    }
  }
  return lrm;
}
// Compute (M,N,N) array of log rates, where the (m,i,j) element is log lambda_{i,j}(t_m) (and is therefore the value of that intensity function since the last time lambda_{i,j} changed).
double llk2(Rcpp::NumericVector lrm,
            Rcpp::NumericVector times,
            Rcpp::IntegerVector sen, 
            Rcpp::IntegerVector rec, 
            int N, int M){
  double llk = 0;
  double den;
  double delta = 0;
  int i,j;
  llk = lrm[threeDIndex(0,sen[0],rec[0],M,N,N)];
  for (int m = 1; m < M; m++) {
    den = 0.0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (i != j) {
          den += exp(lrm[threeDIndex(m,i,j,M,N,N)]);
        }
      }
    }
    llk += lrm[threeDIndex(m,sen[m],rec[m],M,N,N)];
    delta = times[m] - times[m-1];
    llk -= delta * den;
  }
  return llk;
}
RCPP_MODULE(brem){
  function( "llk", &llk ) ;
  function( "llk2", &llk2 ) ;
  function( "lrm", &lrm ) ;
  function( "updateStatistics", &updateStatistics);
  function( "computeLambda", &computeLambda);
}
', plugin="Rcpp")

brem <- Module("brem",getDynLib(fx))
