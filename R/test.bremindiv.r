
source("R/brem.cpp.r")
source("R/utils.r")
library(testthat)
library(abind)
set.seed(1)
M <- 4
N <- 5
times <- c(1,2,3,4)
sen <- c(1,3,3,1)
rec <- c(3,1,1,3)
A <- cbind(times,sen,rec)
K <- 1
beta <- list("intercept"=matrix(1,1,1),
             "abba" = matrix(1,1,1),
             "abby"=matrix(0,1,1),
             "abxa"=matrix(0,1,1),
             "abxb"=matrix(0,1,1),
             "abay"=matrix(0,1,1),
             "abab"=matrix(0,1,1),
             "sod"=matrix(0,1,1),
             "rod"=matrix(0,1,1),
             "sid"=matrix(0,1,1),
             "rid"=matrix(0,1,1))
P <- length(beta)
beta <- abind(beta,rev.along=3)

a <- 1
b <- 3
s <- array(0,c(P,N,N))
#s <- brem$updateStatistics(s,a-1,b-1,N,P)
#brem$computeLambda(1,0,0,0,s,beta,N,K,P)  # 

# Constract log rate matrix by hand and compare to drem$lrm
a <- array(1,c(M,N,N))
a[1,,] <- matrix(0,N,N)
a[2,3,1] <- 2
a[3,1,3] <- 2
a[4,1,3] <- 2
z <- rep(1,N)
K <- 1
#lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
#expect_that(lrm,equals(a))

# Compute log likelihood by hand.  
diag(a[1,,]) <- diag(a[2,,]) <- diag(a[3,,]) <- diag(a[4,,]) <- -Inf
llks <- c(a[1,sen[1],rec[1]],
          a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
          a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
          a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])) )
sum(llks)


i <- 1
indx <- get.indices(A,N)
tau <- precomputeTau(A,N)
system.time(brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau))
brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau)
#sapply(1:N,function(i) brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx))