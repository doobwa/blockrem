context("gibbs")

set.seed(2)
M <- 1000
N <- 10
K <- 3
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
sim <- generate.brem(M,N,beta,z)

A <- sim$edgelist
px <- rep(1,13)
px[13] <- 0
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
s$precompute()

set.seed(1)
N <- 10
P <- 13
times <- c(.1,1,1.5,2,2.1)#seq(.1,.5,by=.1)
sen <- c(1,2,7,3,3)
rec <- c(3,4,8,2,1)
A <- cbind(times,sen,rec)
M <- length(times)
K <- 3
z <- c(1,1,1,2,2,2,3,3,3,3)
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
s$precompute()

a <- 1
z[1] <- 1
RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
z[1] <- 2
RemLogLikelihoodPc(beta,z-1,s$ptr(),K)

lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)
lrm2 <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)

ix <- which(z[A[,2]] == z[a] | z[A[,3]] == z[a])-1
RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
#RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ix)
z[1] <- 2
ix <- which(z[A[,2]] == z[a] | z[A[,3]] == z[a])-1
RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
#RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ix)


## Check ActorPc agrees with using entire dataset
a <- 1
z[a] <- 1
ma <- which(A[,2] == a | A[,3] == a)-1
olp1 <- RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K)$llk
olp2 <- sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
olp3 <- sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ma))
z[a] <- 2
clp1 <- RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K)$llk
clp2 <- sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
clp3 <- sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ma))

olp1-clp1
olp2-clp2
olp3-clp3

#gibbs.collapsed(1:N,beta,z,s$ptr(),px,N,K,nextra=1)
b <- array(rnorm(P*K*K),c(P,K,K))
#gibbs.collapsed(1:N,b,z,s$ptr(),px,N,K,nextra=1)
#fit <- mcmc(A,N,K,px,slice,beta=b,z=z,niter=1)
#b <- slice(b,lposterior,m=20)

set.seed(1)
M <- 7
N <- 5
P <- 13
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
             "rid"=matrix(0,K,K),
             "dc" =matrix(c(1,0,0,2),K,K),
             "cc" =matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)
z <- c(rep(1,N-1),2)
s <- new(RemStat,times,sen-1,rec-1,N,M,P)
s$precompute()


test_that("gibbs runs on small example",{
  b <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
  b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
})

test_that("no -Inf in loglikelihoods",{
  beta <- matrix(rnorm(P*K*K),c(P,K,K))
  b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
  expect_that(brem.lpost.fast(A,N,K,z,s,beta) > -Inf, is_true())
  
#   load("../../../data/synthetic.rdata")
#   s <- new(RemStat,train[,1],train[,2]-1,train[,3]-1,N,M,P)
#   s$precompute()
#   beta <- array(rnorm(P*K*K),c(P,K,K))
#   z <- sample(1:K,N,rep=TRUE)
#   llks <- loglikelihood_fast(beta,z-1,s$ptr(),K)
#   expect_that(all(llks > -Inf), is_true())
#   
#   load("../../../data/eckmann-small.rdata")
#   s <- new(RemStat,train[,1],train[,2]-1,train[,3]-1,N,M,P)
#   s$precompute()
#   beta <- array(rnorm(P*K*K),c(P,K,K))
#   z <- sample(1:K,N,rep=TRUE)
#   llks <- loglikelihood_fast(beta,z-1,s$ptr(),K)
#   expect_that(all(llks > -Inf), is_true())
})


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
