library(ggplot2)

N <- 10
K <- 2
beta <- list("intercept"=matrix(c(0,-1,-1,0),K,K),
             "abba" = matrix(c(3,0,0,0),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(c(0,0,0,2),K,K),
             "abab"=matrix(c(0,-1,-1,0),K,K),
             "sod"=matrix(c(.25,0,0,.25),K,K),
             "rod"=matrix(c(0,0,0,0),K,K),
             "sid"=matrix(c(0,0,0,0),K,K),
             "rid"=matrix(c(.25,0,0,.25),K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)

M <- 1000
set.seed(1)
z <- c(rep(1,5),rep(2,5))

sim <- simulate.brem(M*5,N,z,beta)
mat <- table(sim$A[,2],sim$A[,3])
mat <- melt(as.matrix(mat))
colnames(mat) <- c("X1","X2","value")
plotmat(mat)
plot(sim$A[,1])

A <- sim$A
beta.hat1 <- brem.mle(A,N,K,P,z,beta=beta)
beta.hat2 <- brem.mle(A,N,K,P,z,beta=array(rnorm(P*K*K), c(P, K, K)))
beta.hat3 <- brem.mle(A,N,K,P,z,beta=array(rnorm(P*K*K,0,.05), c(P, K, K)))
sum(brem.llk(A,N,z,beta))
sum(brem.llk(A,N,z,beta.hat1))
sum(brem.llk(A,N,z,beta.hat2))
sum(brem.llk(A,N,z,beta.hat3))

