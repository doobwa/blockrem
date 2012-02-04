source("R/brem.cpp.r")
source("R/utils.r")
require(abind)
require(testthat)
set.seed(1)
M <- 10
N <- 5
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
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
z <- c(rep(1,4),2)

# Make sure gibbs runs
s <- new(bremf$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()
b <- bremf$gibbs(beta,z-1,s$ptr(),K)

A <- cbind(times,sen,rec)
indx <- get.indices(A,N)

lrm <- brem$lrm(beta, times, sen-1, rec-1, z-1, N, M, K, P)
for (i in 1:M) diag(lrm[i,,]) <- -Inf

llk_indiv <- function(a,lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 0:(M-1)) {
    i <- sen[m+1]
    j <- rec[m+1]
    if (i==a | j==a | m==(M-1)) {
      llks[m+1] <- lrm[m+1,i+1,j+1]
      for (r in 0:(N-1)) {
        if (r != i) {
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * sum(exp(lrm[m+1,i+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * sum(exp(lrm[m+1,r+1,i+1]))
        }
        if (r != j) {
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * sum(exp(lrm[m+1,j+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * sum(exp(lrm[m+1,r+1,j+1]))
        }
      }
    }
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  return(llks)
}


a <- llk_indiv(0,lrm,times,sen-1,rec-1)
#b <- bremf$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx)$llks[[1]][[1]]
#expect_that(a,equals(b))

# Timing test
set.seed(1)
M <- 1000
N <- 100
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
z <- sample(1:2,N,replace=TRUE)
K <- 2
A <- cbind(times,sen,rec)
indx <- get.indices(A,N)
#system.time(bremf$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx))
