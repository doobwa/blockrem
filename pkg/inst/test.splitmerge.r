context("splitmerge")
library(brem)
library(testthat)
source("pkg/R/splitmerge.r")
## source("pkg/R/brem.r")
## source("pkg/R/samplers.r")
## source("pkg/R/hmc.r")

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
priors <- list(phi=list(mu=0,sigma=.1),alpha=1)
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
  phi.split <- split_phi(phi,2,3,priors)
  f <- function(x,y) { all(x!=y & cor(x,y) > .99) }
  expect_that(f(phi[,1,2],phi.split[,1,3]),is_true())
  expect_that(f(phi[,2,1],phi.split[,3,1]),is_true())
  expect_that(f(phi[,2,2],phi.split[,2,3]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,2]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,3]),is_true())
})


test_that("brem merge makes similar but not equal clusters",{
  sigma <- .01
  phi.merge <- merge_phi(phi,2,3,sigma,priors)
  f <- function(x,y) { all(x!=y & cor(x,y) > .99) }
  ## expect_that(f(phi.merge[,1,2],phi[,1,3]),is_true())
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
## llk_node <- function(a,phi,z) {
##   RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K)
## }

llk_node_test <- function(a,phi,z) {
  tb <- table(edgelist[,2],edgelist[,3])
  tm <- edgelist[nrow(edgelist),1]
  llk <- 0
  for (i in 1:N) {
    for (j in (1:N)[-i]) {
      lam <- exp(phi[1,z[i],z[j]]) * tm
      llk <- llk + dpois(tb[i,j], lam, log=TRUE)
    }
  }
  return(llk)
}
llk_node <- llk_node_test


lposterior <- function(phi,z,priors) {
  K <- dim(phi)[2]
  phi[-1,,] <- 0
  llk <- llk_node_test(1,phi,z)
  pr.phi <- sum(dnorm(phi[which(phi != 0)],priors$phi$mu,priors$phi$sigma,log=TRUE))
  tb <- table(factor(z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  return(llk + pr.phi + pr.z)
}

# Sample the first dimension of phi
sample_phi <- function(phi,z,lpost,kx=NULL) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  lpost <- function(x,lp=TRUE,lgrad=FALSE) {
    phi[1,k1,k2] <- x
    lposterior_test(phi,z,priors)
  }
  olp <- nlp <- matrix(0,K,K)
  ratios <- c(); r <- 0
  for (k1 in kx) {
    for (k2 in kx) {
      olp[k1,k2] <- lpost(phi[1,k1,k2])
      phi[1,k1,k2] <- slice(phi[1,k1,k2],lpost,m=20)
      nlp[k1,k2] <- lpost(phi[1,k1,k2])
    }
  }
  return(list(phi=phi,olp=olp,nlp=nlp,ratios=ratios))
}


## test_that("Gibbs sampling on small example gives true assignments",{

##   a <- 1
##   truellk <- llk_node(1,phi,z)
##   get.ys <- function(phi,z,llk_node) {
##     z2 <- z
##     ys <- matrix(0,N,K)
##     for (i in 1:N) {
##       for (k in 1:K) {
##         z2 <- z
##         z2[i] <- k
##         ys[i,k] <- llk_node(a,phi,z2)
##       }
##     }
##     return(ys)
##   }

##   ys <- get.ys(phi,z,llk_node)
##   expect_that(apply(ys,1,which.max), equals(z))
##   ys <- get.ys(phi,z,llk_node_test)
##   expect_that(apply(ys,1,which.max), equals(z))
## })

## test_that("restricted gibbs scan finds correct assignments",{
##   S <- 4:9
##   g <- gibbs_restricted(phi,z,S,2,3,llk_node_test)
##   expect_that(all(g$z==z), is_true())
## })


M <- 1000
N <- 9
K <- 3
P <- 13
phi <- array(0,c(P,K,K))
z <- c(1,1,1,2,2,2,3,3,3)
truth <- list(phi=phi,z=z)

set.seed(2)
phi[1,,] <- rnorm(9,0,1)
sim <- generate.brem(M,N,phi,z)
table(sim$edgelist[,2],sim$edgelist[,3])

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
sm <- splitmerge(phi,z,lposterior,llk_node,priors)

source("pkg/R/splitmerge.r")
phi <- truth$phi
z <- truth$z
phi[1,,] <- rnorm(9)
priors$sigma <- 1
z <- sample(1:(dim(phi)[2]),9,rep=T)
fit <- mcmc.blockmodel(lposterior,llk_node,priors,N,P,K)



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
