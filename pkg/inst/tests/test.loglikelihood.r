context("loglikelihood and log intensity check")

set.seed(1)
M <- 100
N <- 10
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
ix <- which(sen==rec)
times <- times[-ix]
sen <- sen[-ix]
rec <- rec[-ix]
M <- length(times)
K <- 2
beta <- list("intercept"=matrix(-1,K,K),
             "abba" = matrix(c(1,2,3,4),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(1,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K),
             "dc"=matrix(1,K,K),
             "cc"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)

s <- new(RemStat,times,sen-1,rec-1,N,M,P)
s$precompute()

test_that("Test that statistics the same when updating vs. grabing from precomputed",{
  failed <- FALSE
  r <- InitializeStatisticsArray(N,P)
  for (m in 0:(M-1)) {
    a <- sen[m+1]-1
    b <- rec[m+1]-1
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i != j) {
          x <- s$get_s(m,i,j)
          y <- r[,i+1,j+1]
          if (!all.equal(x,y)) failed <- TRUE
        }
      }
    }
    r <- UpdateStatisticsArray(r,m,a,b,N,P)
  }
  expect_that(failed,is_false())
})

test_that("Test LogLambda and logLambdaPc the same given same statistics, and intensity arrays match",{

  lam1 <- array(0,c(M,N,N))
  lam2 <- array(0,c(M,N,N))
  for (m in 1:M) {
    for (i in 1:N) {
      for (j in 1:N) {
        s_mij <- s$get_s(m-1,i-1,j-1)
        lam1[m,i,j] <- LogLambdaPc(i-1,j-1,z[i]-1,z[j]-1,s_mij,beta,N,K,P)
        S <- array(0,c(P,N,N))
        S[,i,j] <- s_mij
        lam2[m,i,j] <- LogLambda(i-1,j-1,z[i]-1,z[j]-1,S,beta,N,K,P)
      }
    }
  }
  expect_that(lam1, equals(lam2))
   
  
  lrm1 <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)
  expect_that(lam1, equals(lrm1))
  
  lrm2 <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  expect_that(lam2, equals(lrm2))
})

test_that("Different ways of computing the loglikelihood agree",{
  
  lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)
  
  taus <- test_taus(lrm,times,sen-1,rec-1,M,N)
  taus2 <- test_taus_from_s(times,sen-1,rec-1,N,M,P)
  expect_that(all.equal(taus,taus2),is_true())
  
  llks <- RemLogLikelihood(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  llksa <- RRemLogLikelihoodFromArraySlow(lrm,times,sen-1,rec-1,N,M)
  llk2 <-  RemLogLikelihoodFromArray(lrm,times,sen-1,rec-1,N,M)
  llk2a <-  RemLogLikelihoodVecFromArray(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llksa),equals(llk2))
  expect_that(llksa,equals(llk2a))
  
  for (m in 1:M) diag(lrm[m,,]) <- -Inf
  m <- 20
  a <- RemLogLikelihoodVecFromArray(lrm[1:m,,],times,sen-1,rec-1,N,m)
  b <- RRemLogLikelihoodFromArrayAlt(lrm[1:m,,],times,sen-1,rec-1,N,m)
  expect_that(sum(a), equals(sum(b)))
  
  llk3 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
  llk3a <- RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,0:(M-1))
  llk4 <- RRemLogLikelihoodFromArrayAlt(lrm,times,sen-1,rec-1,N,M)
  eqs <- sapply(1:M,function(m) all.equal(llk3[m],llk4[m]))
  expect_that(llk3,equals(llk4))
  expect_that(llk3a,equals(llk4))
  
  expect_that(sum(llk4),equals(sum(llk2)))
  
})
