context("slice sampling")

require(abind)
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
A <- cbind(times,sen,rec)
priors <- list(beta=list(mu=0,sigma=1))

olp <- brem.lpost.fast(A, N, K, z, s, beta, priors)
k1 <- 1
k2 <- 1
brem.lpost.fast.block(A, N, K, z, s, beta, k1, k2, priors)

nsim <- 100
lps <- rep(0,nsim)
for (i in 1:nsim) {
  px <- 1:12
  b <- brem.slice(A,N,K,P,z,s,beta,px,model.type="baserates",priors)
  beta <- b$current
  lps[i] <- brem.lpost.fast(A, N, K, z, s, beta, priors)
}

load("data/synthetic.rdata")
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,M,P)
s$precompute()
priors <- list(beta=list(mu=0,sigma=1))
olp <- brem.lpost.fast(A, N, K, z, s, beta, priors)
brem.lpost.fast.block(A, N, K, z, s, beta, 1, 1, priors)
nsim <- 50
lps <- rep(0,nsim)
current <- beta * 0
params <- matrix(0,nsim,P*K*K)
for (i in 1:nsim) {
  px <- 1:12
  b <- brem.slice(A,N,K,P,z,s,current,px,model.type="shared",priors)
  current <- b$current
  params[i,] <- c(current)
  lps[i] <- b$olp #brem.lpost.fast(A, N, K, z, s, current, priors)
}


