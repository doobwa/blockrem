source("R/brem.cpp.r")
source("R/utils.r")
require(abind)
require(testthat)
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

# Make sure gibbs runs
s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()
b <- brem$llkfast(beta,z-1,s$ptr(),K)

s$get_w(2,0)
for (m in 0:(M-1)) {
  s$get_s(m,2,0)
}
s$get_s(6,1,3)
st <- proc.time()
b <- brem$gibbs(beta,z-1,s$ptr(),K)
proc.time() - st



# Try on twitter dataset
# load("data/twitter.example.rdata")
# z <- c(rep(1,N-1),2)
# s <- new(bremf$Stat,B[,1],B[,2]-1,B[,3]-1,N,M,P)
# s$precompute()
# 
# st <- proc.time()
# b <- bremf$gibbs(beta,z-1,s$ptr(),K)
# proc.time() - st


lrm <- brem$lrm(beta, times, sen-1, rec-1, z-1, N, M, K, P)
lrm[1,,] <- 1
for (i in 1:M) diag(lrm[i,,]) <- -Inf

llks <- rep(0,M)
for (i in 1:M) {
  llks[i] <- a[i,sen[i],rec[i]] - (times[i]-times[i-1]) * sum(exp(a[i,,]))
}
llks <- c(a[1,sen[1],rec[1]],
          a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
          a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
          a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])) )
sum(llks)



a <- llk_indiv(1,lrm,times,sen-1,rec-1)
#b <- bremf$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx)$llks[[1]][[1]]
#expect_that(a,equals(b))

# Timing test
set.seed(1)
M <- 1000
N <- 500
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
z <- sample(1:2,N,replace=TRUE)
ix <- which(sen==rec)
sen <- sen[-ix]
rec <- rec[-ix]
times <- times[-ix]
K <- 2
s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()
b <- brem$gibbs(beta,z-1,s$ptr(),K)
b <- brem$llkfast(beta,z-1,s$ptr(),K)
