source("R/brem.cpp.r")
source("R/utils.r")
require(abind)
set.seed(1)
M <- 3000
N <- 200
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
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
z <- rep(1,N)

A <- cbind(times,sen,rec)
indx <- get.indices(A,N)
brem$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx)

system.time(brem$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx))

