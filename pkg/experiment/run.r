context("splitmerge")
library(brem)
library(testthat)
source("pkg/R/splitmerge.r")
source("pkg/experiment/utils.r")

M <- 1000
N <- 9
K <- 3
P <- 13
phi <- array(0,c(P,K,K))
z <- c(1,1,1,2,2,2,3,3,3)
set.seed(2)
phi[1,,] <- rnorm(9,0,1)
sim <- generate.brem(M,N,phi,z)
table(sim$edgelist[,2],sim$edgelist[,3])
truth <- list(phi=phi,z=z)

# Set priors
priors <- list(alpha=1,phi=list(mu=0,sigma=1),sigma=.1)

# Test lposterior
edgelist <- sim$edgelist
lposterior(phi,z,priors)
  
options(cores=8)
s <- expand.grid(do.sm = c(TRUE,FALSE),
                 case  = 1:20,
                 sig   = c(.25,.5,.75,1))#c(.1,.3,.5,.8,1,1.2))
s <- s[-which(s$do.sm == FALSE & s$sig != .1)]

source("pkg/R/splitmerge.r")
set.seed(2)
niter <- 50
res <- mclapply(1:nrow(s),function(i) {
  priors$sigma <- s$sig[i]
  mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,do.sm=s$do.sm[i],niter=niter,sigma=s$sig[i])
})
ix <- which(sapply(res,is.character))

## x <- 0

## x <- x+1
## set.seed(x)
## a <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,do.sm=s$do.sm[i],niter=niter)

save(res,ix,s,file="pkg/experiment/res.rdata")
