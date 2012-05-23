context("gibbs")

set.seed(2)
M <- 1000
N <- 10
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
             "cc"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)

set.seed(1)
N <- 10
P <- 13
px <- rep(0,P)
px[1:6] <- 1
sen <- c(1,2,7,3,3,2,5,3)
rec <- c(3,4,8,2,1,1,4,1)
M <- length(sen)
times <- seq(.1,1,length.out=M)
A <- cbind(times,sen,rec)
z <- c(1,1,1,2,2,2,2,2,2,2)
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
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

