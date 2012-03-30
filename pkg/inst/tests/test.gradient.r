
context("gradient")
library(brem)
set.seed(2)
M <- 2000
N <- 500
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
ix <- which(sen==rec)
times <- times[-ix]
sen <- sen[-ix]
rec <- rec[-ix]
K <- 1
beta <- list("intercept"=matrix(-1,K,K),
             "abba" = matrix(1,K,K),
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
z <- rep(1,N)#c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
## context("gradient")
## library(brem)
## library(multicore)
## library(testthat)
## source("pkg/R/brem.r")
## source("pkg/R/samplers.r")
## source("pkg/R/hmc.r")


## set.seed(2)
## M <- 5000
## N <- 10
## K <- 1
## beta <- list("intercept"=matrix(-1,K,K),
##              "abba" = matrix(1,K,K),
##              "abby"=matrix(0,K,K),
##              "abxa"=matrix(0,K,K),
##              "abxb"=matrix(0,K,K),
##              "abay"=matrix(1,K,K),
##              "abab"=matrix(0,K,K),
##              "sod"=matrix(0,K,K),
##              "rod"=matrix(0,K,K),
##              "sid"=matrix(0,K,K),
##              "rid"=matrix(0,K,K),
##              "dc"=matrix(0,K,K),
##              "cc"=matrix(0,K,K))
## z <- rep(1,N)#c(rep(1,N/2),rep(2,N/2))
## P <- length(beta)
## beta <- abind(beta,rev.along=3)
## set.seed(1)
## sim <- generate.brem(M,N,beta,z)
## A <- sim$edgelist
## px <- rep(1,13)
## px[13] <- 0
## s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
## s$precompute()
## s$transform()
## k1 <- k2 <- 1
## priors <- list(beta=list(mu=0,sigma=1))
## llk <- sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
## lp <- block(A,s,beta,z,k1,k2,px)


## test_that("Log posterior for single block works",{
##   expect_that(lposterior(beta[which(px==1),k1,k2]), equals(lp))
## })

## test_that("Samplers return new values",{
##   px <- c(rep(1,6),rep(0,7))
##   current <- beta[which(px==1),k1,k2]
##   set.seed(1)
##   q1 <- mh(current,lposterior,sd=.001)
##   q2 <- slice(current,lposterior,m=20)
##   q3 <- hmc(current,lposterior,lgrad,eps,L)
## })

## M <- length(times)
## s <- new(RemStat,times,sen-1,rec-1,N,M,P)
## s$precompute()
## system.time(RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1)))
## system.time(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
## system.time(RemGibbsPc(beta,z-1,s$ptr(),K))
  
  
##   M <- length(times)
##   s <- new(RemStat,times,sen-1,rec-1,N,M,P)
##   s$precompute()
##   RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))  
  
## beta[12]=.5
## RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))



## k1 <- 1
## k2 <- 1
## px <- rep(1,P)
## px[c(7,13)] <- 0
## px[2:13] <- 0
## A <- cbind(times,sen,rec)
## U(beta[which(px==1),1,1])
## gU(beta[which(px==1),1,1])
## source("pkg/R/hmc.r")
## L <- 10
## eps <- .0001
## niter <- 100
## qs <- matrix(0,niter,sum(px))
## #qs[1,] <- beta[which(px==1),k1,k2]
## for (i in 2:niter) {
##   qs[i,] <- HMC(U, gU, eps, L, qs[i-1,])
##   print(qs[i,])
## }

## ## M <- length(times)
## ## s <- new(RemStat,times,sen-1,rec-1,N,M,P)
## ## s$precompute()
## ## system.time(RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1)))
## ## system.time(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
## ## system.time(RemGibbsPc(beta,z-1,s$ptr(),K))
  
  
## ##   M <- length(times)
## ##   s <- new(RemStat,times,sen-1,rec-1,N,M,P)
## ##   s$precompute()
## ##   RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))  
  
## ## beta[12]=.5
## ## RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))



## ## k1 <- 1
## ## k2 <- 1
## ## px <- rep(1,P)
## ## px[c(7,13)] <- 0
## ## px[2:13] <- 0
## ## A <- cbind(times,sen,rec)
## ## U(beta[which(px==1),1,1])
## ## gU(beta[which(px==1),1,1])
## ## source("pkg/R/hmc.r")
## ## L <- 10
## ## eps <- .0001
## ## niter <- 100
## ## qs <- matrix(0,niter,sum(px))
## ## #qs[1,] <- beta[which(px==1),k1,k2]
## ## for (i in 2:niter) {
## ##   qs[i,] <- HMC(U, gU, eps, L, qs[i-1,])
## ##   print(qs[i,])
## ## }

## ## xs <- seq(-3,3,by=.1)
## ## plot(xs,sapply(xs,function(x) U(x)),type="l")

## ## library(brem)
## ## opts=list(dataset="synthetic",numclusters=2,model.type="full",gibbs=TRUE,numiterations=100,slice=TRUE,initialize=FALSE,fixz=FALSE,skip.intercept=FALSE)
## ## priors <- list(beta=list(mu=0,sigma=1))
## ## load("results/synthetic/full.2.rdata")
## ## load(paste("data/",opts$dataset,".rdata",sep=""))
## ## # Check llk functions
## ## M <- nrow(A)
## ## s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
## ## s$precompute()
## ## U(beta[which(px==1),1,1])
## ## gU(beta[which(px==1),1,1])

## ## xs <- seq(-3,5,by=.2)
## ## px <- rep(0,P)
## ## px[2] <- 1
## ## k1 <- 1
## ## k2 <- 1
## ## par(mfrow=c(2,1))
## ## plot(xs,sapply(xs,function(x) U(x)),type="l")
## ## plot(xs,sapply(xs,function(x) gU(x)),type="l")

## ## source("pkg/R/hmc.r")
## ## L <- 10
## ## eps <- .001
## ## niter <- 100
## ## qs <- matrix(0,niter,sum(px))
## ## #qs[1,] <- beta[which(px==1),k1,k2]
## ## for (i in 2:niter) {
## ##   qs[i,] <- HMC(U, gU, eps, L, qs[i-1,])
## ##   print(qs[i,])
## ## }

## ## brem.slice(A,N,2,P,z,s,beta,px)
