context("splitmerge")
library(brem)
library(testthat)
source("pkg/R/splitmerge.r")

set.seed(2)
M <- 1000
N <- 10
K <- 2
beta <- list("intercept"=matrix(-1,K,K),
             "abba" = matrix(c(1,0,0,2),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K))

z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
dimnames(beta) <- NULL

# Test simple adding, removing and sampling of priors
priors <- list(phi=list(mu=0,sigma=.1),alpha=1,sigma=.1)
phi <- add_cluster(beta)
expect_that(dim(phi),equals(c(P,3,3)))
phi <- sample_cluster_from_prior(phi,3,priors)
expect_that(phi[,1:2,1:2], equals(beta))
expect_that(all(abs(phi[,1:3,3]) < 1), is_true())
expect_that(all(abs(phi[,3,1:3]) < 1), is_true())
expect_that(remove_cluster(phi,5),throws_error())
expect_that(remove_cluster(phi,0),throws_error())
expect_that(sample_cluster_from_prior(phi,5,priors), throws_error())
expect_that(sample_cluster_from_prior(phi,0,priors), throws_error())

test_that("brem split makes similar but not equal clusters",{
  phi <- add_cluster(beta)
  phi.split <- split_phi(phi,2,3,priors)
  f <- function(x,y) { all(x!=y & cor(x,y) > .95) }
  expect_that(f(phi[,1,2],phi.split[,1,3]),is_true())
  expect_that(f(phi[,2,1],phi.split[,3,1]),is_true())
  expect_that(f(phi[,2,2],phi.split[,2,3]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,2]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,3]),is_true())
})


test_that("brem merge makes similar but not equal clusters",{
  sigma <- .01
  phi.merge <- merge_phi(phi,2,3,priors)
  f <- function(x,y) { all(x!=y & cor(x,y) > .95) }
  a <- (phi.merge[,1,2] + phi.merge[,1,3])/2  # center for new location
  expect_that(f(phi.merge[,1,2],a),is_true())
  ## expect_that(f(phi.merge[,2,1],phi[,3,1]),is_true())
  ## expect_that(f(phi.merge[,2,2],phi[,3,3]),is_true())
})

test_that("transition probabilities for phi work",{
  phi <- beta
  phi <- add_cluster(phi)
  phi <- add_cluster(phi)
  sigma <- priors$phi$sigma
  expect_that(ps2pm(phi,phi,3,4,sigma,priors),equals(0))
  phi <- remove_cluster(phi,4)
  phi <- sample_cluster_from_prior(phi,2,priors)
  phi <- sample_cluster_from_prior(phi,3,priors)
  phi.merge <- phi
  phi <- sample_cluster_from_prior(phi,2,priors)
  phi <- sample_cluster_from_prior(phi,3,priors)
  phi.split <- phi
  exp(ps2pm(phi.split,phi.merge,2,3,.2,priors) - 
      pm2ps(phi.merge,phi.split,2,3,.2,priors))
})



## # TODO: Bug in using the ActorPc function like this.
llk_node <- function(a,phi,z) {
  RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K)
}
lposterior <- function( ) {
  RemLogLikelihoodPc(phi,z-1,s$ptr(),K)
}

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

phi <- truth$phi
z <- truth$z

test_that("restricted gibbs scan finds correct assignments",{
  S <- 4:9
  g <- gibbs_restricted(phi,z,S,2,3,llk_node)
  expect_that(all(g$z==z), is_true())
})

test_that("restricted gibbs can use predetermined assignments",{
  S <- 4:7
  z[5] <- 3
  g <- gibbs_restricted(phi,z,S,2,3,llk_node,prob.only=TRUE)
  # z[5] not gibbs sampled
  expect_that(g$z[5], equals(3))
  # all of transition probability comes from z[5] assigned to 3, not 2
  expect_that(log(g$ys[5,3]), equals(g$transition))
}


# Get edgelist and precompute datastructures
edgelist <- A <- sim$edgelist
px <- rep(0,13)
px[1] <- 1
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
s$precompute()
s$transform()

# Set priors and test lposterior, sampler
priors <- list(alpha=1,phi=list(mu=0,sigma=1),sigma=.1)
lposterior(phi,z,priors)
k1 <- k2 <- 1
lposterior_test <- lposterior
r <- sample_phi(phi,z,lposterior)

## Debug split merge
set.seed(4)
sm <- splitmerge(phi,z,lposterior,llk_node,priors)

source("pkg/R/splitmerge.r")
phi <- truth$phi
z <- truth$z
fit <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,niter=20,verbose=TRUE)

priors$sigma <- 1
fit <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,niter=20,verbose=FALSE)

## Small experiment
          
options(cores=8)
s <- expand.grid(do.sm = c(TRUE,FALSE),
                 case  = 1:50,
                 sig   = c(.1,.3,.5,.8,1,1.2))
s <- s[-which(s$do.sm == FALSE & s$sig != .1)]

set.seed(2)
res <- mclapply(1:5,function(i) {
  priors$sigma <- s$sig[i]
  mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,do.sm=s$do.sm[i])
})
ix <- which(sapply(res,is.character))

res <- res[-ix]
iy <- (1:nrow(s))[-ix]
save(res,iy,file="working.data/res.rdata")
load("working.data/res.rdata")
## res <- list()
## for (i in 1:nrow(s)) {
##   priors$sigma <- s$sig[i]
##   set.seed(s$case[i])
##   res[[i]] <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,K,do.sm=s$do.sm[i],niter=100)
## }


3
## test_that("gibbs runs on small example",{
##   b <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
##   b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
## })

## test_that("no -Inf in loglikelihoods",{
##   beta <- matrix(rnorm(P*K*K),c(P,K,K))
##   b <- RemGibbsPc(1:N-1,beta,z-1,s$ptr(),K)
##   expect_that(brem.lpost.fast(A,N,K,z,s,beta) > -Inf, is_true())
  
#   load("../../../data/synthetic.rdata")
#   s <- new(RemStat,train[,1],train[,2]-1,train[,3]-1,N,M,P)
#   s$precompute()
#   beta <- array(rnorm(P*K*K),c(P,K,K))
#   z <- sample(1:K,N,rep=TRUE)
#   llks <- loglikelihood_fast(beta,z-1,s$ptr(),K)
#   expect_that(all(llks > -Inf), is_true())
#   
#   load("../../../data/eckmann-small.rdata")
#   s <- new(RemStat,train[,1],train[,2]-1,train[,3]-1,N,M,P)
#   s$precompute()
#   beta <- array(rnorm(P*K*K),c(P,K,K))
#   z <- sample(1:K,N,rep=TRUE)
#   llks <- loglikelihood_fast(beta,z-1,s$ptr(),K)
#   expect_that(all(llks > -Inf), is_true())
#})


# 
# # TODO: Test gibbs probabilities
# lrm <- brem$lrm(beta, times, sen-1, rec-1, z-1, N, M, K, P)
# lrm[1,,] <- 1
# for (i in 1:M) diag(lrm[i,,]) <- -Inf
# 
# llks <- rep(0,M)
# for (i in 1:M) {
#   llks[i] <- a[i,sen[i],rec[i]] - (times[i]-times[i-1]) * sum(exp(a[i,,]))
# }
# sum(llks)

# a <- llk_indiv(1,lrm,times,sen-1,rec-1)
# b <- bremf$gibbs(beta,times,sen-1,rec-1,z-1,N,M,K,P,indx)$llks[[1]][[1]]
# expect_that(a,equals(b))
