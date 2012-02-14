context("gibbs")

require(abind)
set.seed(1)
M <- 7
N <- 5
P <- 11
times <- seq(.1,.7,by=.1)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)
K <- 2

beta <- list("intercept"=matrix(1,K,K),
             "abba" = matrix(c(1,2,3,4),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(0,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)
z <- c(rep(1,N-1),2)

test_that("gibbs runs on small example",{
  s <- new(RemStat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  b <- loglikelihood_fast(beta,z-1,s$ptr(),K)
  b <- gibbs(beta,z-1,s$ptr(),K)
})

test_that("gibbs runs with K>2",{
  s <- new(RemStat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  
  K <- 5
  beta <- array(rnorm(K*K*P),c(P,K,K))
  z <- sample(1:K,N,replace=TRUE)
  b <- loglikelihood_fast(beta,z-1,s$ptr(),K)
  b <- gibbs(beta,z-1,s$ptr(),K)
  
  
  K <- 5
  z <- c(4,4,4,4,4,5,5,5,5,5)
  b <- loglikelihood_fast(beta,z-1,s$ptr(),K)
  b <- gibbs(beta,z-1,s$ptr(),K)
})

library(brem)
load("data/synthetic.rdata")
set.seed(2)
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,M,P)
s$precompute()

K <- 5
load("results/synthetic/full.10.rdata")
b <- loglikelihood_fast(res$beta,res$z-1,s$ptr(),K)
b <- gibbs(res$beta,res$z-1,s$ptr(),K)



# 
# # TODO: Test gibbs probabilities
# lrm <- brem$lrm(beta, times, sen-1, rec-1, z-1, N, M, K, P)
# lrm[1,,] <- 1
# for (i in 1:M) diag(lrm[i,,]) <- -Inf
# 
# llks <- rep(0,M)
# for (i in 1:M) {
#   llks[i] <- a[i,sen[i],rec[i]] - (times[i]-times[i-1]) * sum(exp(a[i,,]))
# }
# sum(llks)

# a <- llk_indiv(1,lrm,times,sen-1,rec-1)
# b <- bremf$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx)$llks[[1]][[1]]
# expect_that(a,equals(b))
