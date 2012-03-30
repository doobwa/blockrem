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
  phi.split <- split_phi(phi,2,3,sigma=.001)
  f <- function(x,y) { all(x!=y & cor(x,y) > .99) }
  expect_that(f(phi[,1,2],phi.split[,1,3]),is_true())
  expect_that(f(phi[,2,1],phi.split[,3,1]),is_true())
  expect_that(f(phi[,2,2],phi.split[,2,3]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,2]),is_true())
  expect_that(f(phi[,2,2],phi.split[,3,3]),is_true())
})


test_that("brem merge makes similar but not equal clusters",{
  phi.merge <- merge_phi(phi,2,3)
  f <- function(x,y,z) all(x %in% y | x %in% z)# & !(all(x==y) | all(x==z))
  expect_that(f(phi.merge[,1,2],phi[,1,3],phi[,1,2]),is_true())
  expect_that(f(phi.merge[,2,1],phi[,3,1],phi[,2,1]),is_true())
  expect_that(f(phi.merge[,2,2],phi[,3,3],phi[,2,2]),is_true())
})

llk_node <- function(a,phi,z) {
  RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K)
}


M <- 5000
N <- 9
K <- 3
P <- 13
phi <- array(0,c(P,K,K))
z <- c(1,1,1,2,2,2,3,3,3)

set.seed(2)
phi[1,,] <- rnorm(9,0,1)
sim <- generate.brem(M,N,phi,z)
table(sim$edgelist[,2],sim$edgelist[,3])

A <- sim$edgelist
px <- rep(0,13)
px[1] <- 1
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
s$precompute()
s$transform()

llk_node <- function(a,phi,z) {
  RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K)
}

llk_node_test<- function(a,phi,z) {
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

test_that("Gibbs sampling on small example gives true assignments",{
  truellk <- llk_node(1,phi,z)
  get.ys <- function(phi,z,llk_node) {
    z2 <- z
    ys <- matrix(0,N,K)
    for (i in 1:N) {
      for (k in 1:K) {
        z2 <- z
        z2[i] <- k
        ys[i,k] <- llk_node(a,phi,z2)
      }
    }
    return(ys)
  }

  ys <- get.ys(phi,z,llk_node)
  expect_that(apply(ys,1,which.max), equals(z))
  ys <- get.ys(phi,z,llk_node_test)
  expect_that(apply(ys,1,which.max), equals(z))
})

g <- gibbs_restricted(phi,z,2,3,llk_node_test)
g

lposterior_test <- function(phi,z,priors) {
  llk <- llk_node_test(1,phi,z)
  pr.phi <- sum(dnorm(phi[which(phi != 0)],priors$phi$mu,priors$phi$sigma,log=TRUE))
  tb <- table(factor(z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  return(llk + pr.phi + pr.z)
}

priors <- list(alpha=1,phi=list(mu=0,sigma=1))
lposterior_test(phi,z,priors)

source("pkg/R/splitmerge.r")
s <- splitmerge(phi,z,lposterior_test,llk_node_test,priors)

truth <- list(phi=phi,z=z)


# Still bugs in trying this...
for (a in 1:3) {
  s <- splitmerge(phi,z,lposterior_test,llk_node_test,priors)
  cat("init ",z,"\n")
  cat(s$final$type," ",s$final$z,"(prob ",exp(s$final$prob),")\n")
  alpha <- min(1,s$final$prob)
  if (s$final$prob < runif(1)) {
    phi <- s$final$phi
    z <- s$final$z
  }
  K <- dim(phi)[2]
  tb <- table(factor(z,1:K))
  for (k in which(tb==0)) phi <- remove_cluster(phi,k)
}


## sim <- generate.brem(M,N,beta,z)

A <- sim$edgelist
px <- rep(1,13)
px[13] <- 0
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
s$precompute()
s$transform()

## gibbs.collapsed(1:N,beta,z,s$ptr(),nextra=1)
## b <- array(rnorm(P*K*K),c(P,K,K))
## gibbs.collapsed(1:N,b,z,s$ptr(),nextra=1)
## fit <- mcmc(A,N,K,px,slice,beta=b,z=z,niter=1)
## b <- slice(b,lposterior,m=20)
## d


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
