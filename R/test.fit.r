library(ggplot2)
source("rem.cpp.r")
source("fns.r")
N <- 10
K <- 1
z <- rep(1,N)
P <- 7
beta <- array(0,c(K,K,P))
beta[1,1,] <- c(0,3,0,0,0,2,0)

M <- 5000
set.seed(1)
sim <- simulate.brem(M,N,z,beta)

drem.mle(sim$A,beta,1:N,1:N)
