context("gradient")
library(brem)
set.seed(2)
M <- 10
N <- 6
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


brem.lpost.block <- function(A,N,K,z,s,beta,k1,k2,px,priors=list(beta=list(mu=0,sigma=1))) { 
  beta[which(px==0),,] <- 0
  M <- nrow(A)
  zs <- z[A[,2]]
  zr <- z[A[,3]]
  ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
  if (length(ix) == 0) return(0)
  sum(dnorm(unlist(beta[which(px==1),k1,k2]),priors$beta$mu,priors$beta$sigma,log=TRUE)) +  sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ix))
}
brem.lpost.block.grad <- function(A,N,K,z,s,beta,k1,k2,px,priors=list(beta=list(mu=0,sigma=1))) {   
  M <- nrow(A)
  zs <- z[A[,2]]
  zr <- z[A[,3]]
  ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
  if (length(ix) == 0) return(0)
  - 2 * (beta[which(px==1),k1,k2] - priors$beta$mu)/priors$beta$sigma + RemGradientPcSubset(beta,z-1,s$ptr(),K,ix)[which(px==1)]
}
U <- function(betak1k2) {
  beta[which(px==1),k1,k2] <- betak1k2
  - brem.lpost.block(A,N,K,z,s,beta,k1,k2,px)
}
gU <- function(betak1k2) {
  beta[which(px==1),k1,k2] <- betak1k2
  - brem.lpost.block.grad(A,N,K,z,s,beta,k1,k2,px)
}
RRemGradient <- function(lrm,times,sen,rec,N,M,P) {
  
  s <- InitializeStatisticsArray(N,P);
  for (m in 2:M) {
    
    i <- A[m-1,2]
    j <- A[m-1,3]
    
    s <- UpdateStatisticsArray(s,m-1,i-1,j-1,N,P)
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  llks[1] <- lrm[1,sen[1]+1,rec[1]+1]
  for (m in 0:(M-2)) {
    i <- sen[m+1]
    j <- rec[m+1]
    llks[m+1] <- lrm[m+1,i+1,j+1]
    for (r in 0:(N-1)) {
      if (r != i & r != j) {
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * exp(lrm[m+1,i+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * exp(lrm[m+1,r+1,i+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * exp(lrm[m+1,j+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * exp(lrm[m+1,r+1,j+1])
      }
    }
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * exp(lrm[m+1,i+1,j+1])
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * exp(lrm[m+1,j+1,i+1])
    
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  llks[m+1] <- lrm[m+1,sen[m+1]+1,rec[m+1]+1]
  other <- 0
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        other <- other - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  llks[m+1] <- llks[m+1] + other
  return(llks)
}


M <- length(times)
s <- new(RemStat,times,sen-1,rec-1,N,M,P)
s$precompute()
RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))
beta[12]=.5
RemGradientPcSubset(beta,z-1,s$ptr(),K,1:(M-1))



k1 <- 1
k2 <- 1
px <- rep(1,P)
px[c(7,13)] <- 0
px[2:13] <- 0
A <- cbind(times,sen,rec)
U(beta[which(px==1),1,1])
gU(beta[which(px==1),1,1])
source("pkg/R/hmc.r")
L <- 10
eps <- .0001
niter <- 100
qs <- matrix(0,niter,sum(px))
#qs[1,] <- beta[which(px==1),k1,k2]
for (i in 2:niter) {
  qs[i,] <- HMC(U, gU, eps, L, qs[i-1,])
  print(qs[i,])
}

xs <- seq(-3,3,by=.1)
plot(xs,sapply(xs,function(x) U(x)),type="l")

library(brem)
opts=list(dataset="synthetic",numclusters=2,model.type="full",gibbs=TRUE,numiterations=100,slice=TRUE,initialize=FALSE,fixz=FALSE,skip.intercept=FALSE)
priors <- list(beta=list(mu=0,sigma=1))
load("results/synthetic/full.2.rdata")
load(paste("data/",opts$dataset,".rdata",sep=""))
# Check llk functions
M <- nrow(A)
s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
s$precompute()
U(beta[which(px==1),1,1])
gU(beta[which(px==1),1,1])

xs <- seq(-3,5,by=.2)
px <- rep(0,P)
px[2] <- 1
k1 <- 1
k2 <- 1
par(mfrow=c(2,1))
plot(xs,sapply(xs,function(x) U(x)),type="l")
plot(xs,sapply(xs,function(x) gU(x)),type="l")

source("pkg/R/hmc.r")
L <- 10
eps <- .001
niter <- 100
qs <- matrix(0,niter,sum(px))
#qs[1,] <- beta[which(px==1),k1,k2]
for (i in 2:niter) {
  qs[i,] <- HMC(U, gU, eps, L, qs[i-1,])
  print(qs[i,])
}

brem.slice(A,N,2,P,z,s,beta,px)