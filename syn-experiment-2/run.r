context("splitmerge")
library(brem)
library(testthat)
source("pkg/R/splitmerge.r")
source("syn-experiment-2/utils.r")


set.seed(1)
M <- 2000
N <- 20
K <- 2
beta <- list("intercept"=matrix(-1,K,K),
             "abba" = matrix(c(1,2,3,4),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(1,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K),
             "dc"=matrix(0,K,K),
             "cc"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
phi <- abind(beta,rev.along=3)

sim <- generate.brem(M,N,phi,z)
s <- new(RemStat,sim$edgelist[,1],sim$edgelist[,2]-1,sim$edgelist[,3]-1,N,M,P)
s$precompute()

# Set priors
priors <- list(alpha=1,phi=list(mu=0,sigma=2),sigma=.1)

# Create lposterior function
truth <- list(phi=phi,z=z)
llk <- sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))

px <- c(1,2,6)
system.time(lp(phi,z,priors))
lp(phi,z,priors)
phi <- sample_phi(phi,z,lp,priors,px=c(1,2,6,12))$phi
lp(phi,z,priors)

phi <- sample_phi(phi,z,lp,priors,px=px)$phi

## phi <- truth$phi
## z <- RemGibbsPc(1:N-1,phi,z-1,s$ptr(),K)$z + 1

priors$sigma <- .5
fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,px=c(1,2,6),do.sm=FALSE,num.extra=3,niter=20,sigma=.1,verbose=TRUE)

phi[c(1,2,6,12),1:2,1:2]
split$phi[c(1,2,6,12),1:2,1:2]
merge$phi[c(1,2,6,12),1:2,1:2]
lp(split$phi,split$z,priors)
lp(merge$phi,merge$z,priors)
  
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
  mcmc.blockmodel(lp,lp_node,priors,N,P,K,do.sm=s$do.sm[i],do.extra=s$do.extra[i],niter=niter)
})
ix <- which(sapply(res,is.character))

save(res,ix,s,niter,file="pkg/experiment/res.rdata")
