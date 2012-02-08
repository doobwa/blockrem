#!Rscript

#' @filename rdata filename with event history A
#' @K number of clusters
#' @niter number of iterations
#' @model.type (shared, full, single, baseline)

source("R/brem.r")
source("R/brem.cpp.r")
source("R/sbm.r")
source("R/utils.r")

load("data/sim.rdata")

# Precompute data structures
N <- max(c(sen,rec))
M <- nrow(A)
P <- 11
s <- new(brem$Stat,A[,1],A[,2]-1,A[,3]-1,N,M,P)
s$precompute()

# fit0 <- sbm.mcmc(sim$A,N,K,niter=niter,z=z,gibbs=TRUE)
fit <- brem.mcmc(A,N,K,s,model.type=model.type,niter=niter,gibbs="fast")
