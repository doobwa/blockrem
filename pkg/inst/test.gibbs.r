context("likelihoods for gibbs sampling")

set.seed(1)
M <- 100
N <- 100
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
             "abba" = matrix(c(1,0,0,2),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(1,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(c(.01,0,.01,0),K,K),
             "dc"=matrix(c(0,.03,0,0),K,K),
             "cc"=matrix(0,K,K),
             "rrs"=matrix(0,K,K),
             "rss"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
px <- rep(0,P)
px[1:6] <- 1

A <- cbind(times,sen,rec)
ego <- 0
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),ego)
s$precompute()

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

test_that("Functions run",{
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


## Work with LogNormalizing function.
## library(plyr)
## k <- 1
## l <- 1
## # Debug just lognormalizing
## a <- ldply(2:20,function(m) {
##   knodes <- which(z==k)
##   lnodes <- which(z==l)
##   beta[1:3,k,l] <- c(2,1,0)
##   a <- LogNormalizing(beta,z-1,s$ptr(),K,m-1,sen[m]-1,rec[m]-1,knodes-1,lnodes-1)
##   beta[1:3,k,l] <- c(0,0,0)
##   b <- LogNormalizing(beta,z-1,s$ptr(),K,m-1,sen[m]-1,rec[m]-1,knodes-1,lnodes-1)
##   return(data.frame(before=a,after=b))
## })

## lrm <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)

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

#gibbs.collapsed(1:N,beta,z,s$ptr(),px,N,K,nextra=1)
#b <- array(rnorm(P*K*K),c(P,K,K))
#gibbs.collapsed(1:N,b,z,s$ptr(),px,N,K,nextra=1)
#fit <- mcmc(A,N,K,px,slice,beta=b,z=z,niter=1)
#b <- slice(b,lposterior,m=20)

## set.seed(1)
## M <- 7
## N <- 5
## P <- 13
## times <- seq(.1,.7,by=.1)
## sen <- c(1,3,3,1,2,5,2)
## rec <- c(3,1,1,3,5,1,4)
## K <- 2

## beta <- list("intercept"=matrix(1,K,K),
##              "abba" = matrix(c(1,2,3,4),K,K),
##              "abby"=matrix(0,K,K),
##              "abxa"=matrix(0,K,K),
##              "abxb"=matrix(0,K,K),
##              "abay"=matrix(0,K,K),
##              "abab"=matrix(0,K,K),
##              "sod"=matrix(0,K,K),
##              "rod"=matrix(0,K,K),
##              "sid"=matrix(0,K,K),
##              "rid"=matrix(0,K,K),
##              "dc" =matrix(c(1,0,0,2),K,K),
##              "cc" =matrix(0,K,K))
## P <- length(beta)
## beta <- abind(beta,rev.along=3)
## z <- c(rep(1,N-1),2)
## s <- new(RemStat,times,sen-1,rec-1,N,M,P)
## s$precompute()

## test_that("gibbs runs on small example",{
##   b <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
##   b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
## })

## test_that("no -Inf in loglikelihoods",{
##   beta <- matrix(rnorm(P*K*K),c(P,K,K))
##   b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
##   expect_that(brem.lpost.fast(A,N,K,z,s,beta) > -Inf, is_true())
## })

