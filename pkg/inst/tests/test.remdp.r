context("remdp testing")
library(brem)
set.seed(1)
M <- 100
N <- 10
P <- 13
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE) 
rec <- sample(1:N,M,replace=TRUE)
ix <- which(sen==rec)
times <- times[-ix]
sen <- sen[-ix]
rec <- rec[-ix]
M <- length(times)
K <- 2
beta <- matrix(0,K,P-1)
gamma <- matrix(0,K,P-1)
eta <- matrix(0,K,K)

colnames(beta) <- colnames(gamma) <- c("intercept","abba","abby","abxa","abxb","abay","abab","sod","rod","sid","rid","dc","cc")[-1]
z <- c(rep(1,N/2),rep(2,N/2))

## test_that("RemDP can be initialized and data structure precomputed",{
##   model <- new(RemDP,times,sen-1,rec-1,N,M,P)
##   model$precompute()
## })

test_that("RemDP constructed via edgelist",{
  edgelist <- data.frame(time=times,s=sen-1,r=rec-1)
  model <- new(RemDP,edgelist,N)
  model$precompute()
  expect_that(edgelist,equals(model$get_edgelist()))
})

test_that("Statistics data structure is accessible",{
  edgelist <- data.frame(time=times,s=sen-1,r=rec-1)
  model <- new(RemDP,edgelist,N)
  model$precompute()

  s <- c(1,0,0,0,0,0,1,1,0,2,1,1,2)
  expect_that(model$get_s(3,2,5),equals(s))
})

test_that("Getting and setting of parameters works",{
  p <- list(eta=eta,beta=beta,gamma=gamma,z=z)
  p$beta[1,2] <- 3
  p$gamma[2,6] <- 3
  p$eta[1,2] <- -1
  p$eta[2,1] <- -1
  edgelist <- data.frame(time=times,s=sen-1,r=rec-1)
  model <- new(RemDP,edgelist,N)
  model$set_params(p)
  expect_that(model$get_params(), equals(p))
})

test_that("Dimension checking works",{
  K <- 2
  N <- 10
  p <- list(eta=matrix(0,K,K),
            beta=matrix(0,K,P),
            gamma=matrix(0,K,P),
            z=rep(0,N))
  edgelist <- data.frame(time=times,s=sen-1,r=rec-1)
  model <- new(RemDP,edgelist,N)
  expect_that(model$Check(), is_false())
  model$set_params(p)
  expect_that(model$Check(), is_true())
  p$eta <- matrix(0,3,3)
  model$set_params(p)
  expect_that(model$Check(), is_false())
  p$eta <- matrix(0,K,K)
  p$z[1] <- K+1
  model$set_params(p)
  expect_that(model$Check(), is_false())  
})

# NB: m,i,j are 0-indexed
RLogLambda <- function(model,m,i,j) { 
  s <- model$get_s(m,i,j)[-1]
  s[7:11] <- log((s[7:11]+1) / (s[12] + N*(N-1)))
  p <- model$get_params()
  lam <- p$eta[p$z[i+1]+1,p$z[j+1]+1]
  lams <- (p$beta[p$z[i+1]+1,] + p$gamma[p$z[j+1]+1,]) * s
  return(sum(lams) + lam)
} 

test_that("LogLambda is correct",{
  m <- 3
  i <- 2
  j <- 5
  p <- list(eta=eta,beta=beta,gamma=gamma,z=z)
  p$beta[1,10] <- 3
  p$gamma[2,6] <- 3
  p$eta[1,2] <- -1
  p$eta[2,1] <- -1
  p$z <- rep(0,N)
  edgelist <- data.frame(time=times,s=sen-1,r=rec-1)
  model <- new(RemDP,edgelist,N)
  model$precompute()
  model$set_params(p)
  lam1 <- RLogLambda(model,m,i,j)
  lam2 <- model$LogLambda(m,i,j)
  expect_that(lam1,equals(lam2))
})

test_that("Likelihood is correct",{
  p <- list(eta=eta,beta=beta,gamma=gamma,z=z)
  p$beta[1,2] <- 3
  p$gamma[2,6] <- 3
  p$eta[1,2] <- -1
  p$eta[2,1] <- -1
  
})

predict.Rcpp_RemDP <- function(model,test,mode="llk") {
  
}

test_that("Prediction on test set",{

  p <- list()
  model -> new(RemDP,train,N)
  model$precompute()
  model$set_params(p)
  res -> predict(model,test,mode="multinomial")  # likelihood on test set
  summary(res)
})

## test_that("Different ways of computing the loglikelihood agree",{
  
##   lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)
  
##   taus <- test_taus(lrm,times,sen-1,rec-1,M,N)
##   taus2 <- test_taus_from_s(times,sen-1,rec-1,N,M,P)
##   expect_that(all.equal(taus,taus2),is_true())
  
##   llks <- RemLogLikelihood(beta,times,sen-1,rec-1,z-1,N,M,K,P)
##   llksa <- RRemLogLikelihoodFromArraySlow(lrm,times,sen-1,rec-1,N,M)
##   llk2 <-  RemLogLikelihoodFromArray(lrm,times,sen-1,rec-1,N,M)
##   llk2a <-  RemLogLikelihoodVecFromArray(lrm,times,sen-1,rec-1,N,M)
##   expect_that(sum(llksa),equals(llk2))
##   expect_that(llksa,equals(llk2a))
  
##   for (m in 1:M) diag(lrm[m,,]) <- -Inf
##   m <- 20
##   a <- RemLogLikelihoodVecFromArray(lrm[1:m,,],times,sen-1,rec-1,N,m)
##   b <- RRemLogLikelihoodFromArrayAlt(lrm[1:m,,],times,sen-1,rec-1,N,m)
##   expect_that(sum(a), equals(sum(b)))
  
##   llk3 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
##   llk3a <- RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,0:(M-1))
##   llk4 <- RRemLogLikelihoodFromArrayAlt(lrm,times,sen-1,rec-1,N,M)
##   eqs <- sapply(1:M,function(m) all.equal(llk3[m],llk4[m]))
##   expect_that(llk3,equals(llk4))
##   expect_that(llk3a,equals(llk4))
  
##   expect_that(sum(llk4),equals(sum(llk2)))
  
## })
