context("splitmerge")
library(brem)
library(testthat)
source("pkg/R/splitmerge.r")
source("syn-experiment-1/utils.r")

M <- 2000
N <- 30
K <- 5
P <- 13
phi <- array(0,c(P,K,K))
z <- rep(1:K,each=N/K)
set.seed(2)
phi[1,,] <- rnorm(K^2,0,1)
sim <- generate.brem(M,N,phi,z)
y <- table(sim$edgelist[,2],sim$edgelist[,3])

truth <- list(phi=phi,z=z)

# Set priors
priors <- list(alpha=1,phi=list(mu=0,sigma=1),sigma=.1)

# Test lposterior
edgelist <- sim$edgelist
lposterior(phi,z,priors)

fit <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,2,do.sm=TRUE,do.extra=TRUE,niter=20,sigma=1,verbose=TRUE)
  
options(cores=8)
s <- expand.grid(do.sm = c(TRUE,FALSE),
                 do.extra = c(TRUE,FALSE),
                 case  = 1:20,
                 sig   = c(.5,1))
#s <- s[-which(s$do.sm == FALSE & s$sig != .25),]

source("pkg/R/splitmerge.r")
set.seed(2)
niter <- 100
res <- mclapply(1:nrow(s),function(i) {
  priors$sigma <- s$sig[i]
  mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,do.sm=s$do.sm[i],do.extra=s$do.extra[i],niter=niter,sigma=s$sig[i])
})
ix <- which(sapply(res,is.save))

character(res,ix,s,niter,file="synthetic-experiment-2/res.rdata")
