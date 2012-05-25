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
             "cc"=matrix(0,K,K),
             "rrs"=matrix(0,K,K),
             "rss"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
A <- cbind(times,sen,rec)

ego <- 0
s <- new(RemStat,times,sen-1,rec-1,N,M,ego)
s$precompute()

test_that("Test that statistics the same when updating vs. grabbing from precomputed",{
  failed <- FALSE
  r <- InitializeStatisticsArray(N,P)
  for (m in 0:(M-1)) {
    a <- sen[m+1]-1
    b <- rec[m+1]-1
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i != j) {
          x <- s$get_s(m,i,j)[1:13]  # TODO: Include recency in this test
          y <- r[,i+1,j+1][1:13]
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



test_that("Check ActorPc agrees with using entire dataset", {
  for (a in 1:N) {
    z[a] <- 1
    o1 <- RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K)
    o2 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
    z[a] <- 2
    c1 <- RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K)
    c2 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
    o1-c1
    o2-c2
    expect_that(o1-c1, equals(o2-c2))
  }
})

test_that("Precomputed likelihoods and normalizing function run",{
  m <- 3
  k <- 1
  l <- 2

  knodes <- which(z==k)
  lnodes <- which(z==l)
  m <- 3
  LogNormalizing(beta,z-1,s$ptr(),K,5-1,sen[m]-1,rec[m]-1,1:N-1,1:N-1)
  RemLogLikelihoodBlockPc(k-1,l-1,knodes-1,lnodes-1,beta,z-1,s$ptr(),K)
  RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
  RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,m-1)
})

test_that("Block version is faster", {
  library(rbenchmark)
  k <- l <- 2
  knodes <- which(z==k)
  lnodes <- which(z==l)
  b <- benchmark(block = RemLogLikelihoodBlockPc(k-1,l-1,knodes-1,lnodes-1,beta,z-1,s$ptr(),K),
                 full  = RemLogLikelihoodPc(beta,z-1,s$ptr(),K),
                 replications = 10)
  expect_that(b$elapsed[1] < b$elapsed[2], is_true())
})


test_that("Check BlockPc agrees with using entire dataset", {
  for (k in 1:K) {
    for (l in 1:K) {
      beta[1:3,k,l] <- c(2,1,0)
      knodes <- which(z==k)
      lnodes <- which(z==l)
      o1 <- RemLogLikelihoodBlockPc(k-1,l-1,knodes-1,lnodes-1,beta,z-1,s$ptr(),K)
      o2 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
      beta[1:3,k,l] <- c(0,0,0)
      c1 <- RemLogLikelihoodBlockPc(k-1,l-1,knodes-1,lnodes-1,beta,z-1,s$ptr(),K)
      c2 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
                                        #    o1-c1
                                        #    o2-c2
      cbind(z[A[,2]],z[A[,3]],o1-c1,o2-c2)[1:20,]
      expect_that(o1-c1, equals(o2-c2))
    }
  }
})
