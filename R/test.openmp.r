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

i <- 1
A <- cbind(times,sen,rec)

cat("precomputing indices and tau")
indx <- get.indices(A,N)
tau <- precomputeTau(A,N)

cat("llk")
brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau)
brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau)

brem$gibbs(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau)


system.time(brem$llki(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau))
system.time(brem$gibbs(i-1,beta,times,sen-1,rec-1,z-1,N,M,K,P,indx[[i]],tau))

# lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
# llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)

#system.time(llk1 <- brem$llk(beta,times,sen-1,rec-1,z-1,N,M,K,P))
#system.time(llk2 <- brem$llkp(beta,times,sen-1,rec-1,z-1,N,M,K,P))



#llk <- brem$llk(beta,times,sen-1,rec-1,z-1,N,M,K,P)
#llkp <- brem$llkp(beta,times,sen-1,rec-1,z-1,N,M,K,P)