context("online evaluation methods")

set.seed(1)
M <- 1000
N <- 10
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
             "dc"=matrix(.1,K,K),
             "cc"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
set.seed(1)
sim <- generate.brem(M,N,beta,z)
A <- sim$edgelist

train <- A[1:500,]
test <- A

## Precompute data structures
ego <- 0
strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P,ego)
strain$precompute()
stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P,ego)
stest$precompute()

fit <- list(beta=beta,z=z)

## Precompute rate arrays
lrm <- list()
lrm$train <- brem.lrm.fast(strain, fit$z, fit$beta)
lrm$test  <- brem.lrm.fast(stest, fit$z, fit$beta)

test_that("train and test BREM likelihood agree with online version", {
  llk.train <- RemLogLikelihoodVecFromArray(lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
  llk.test <- RemLogLikelihoodVecFromArray(lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

  for (i in 2:10) {
    expect_that(event.llk(train,N,fit,i,strain), equals(llk.train[i]))
  }

})

## train.ix <- 1:500
## test.ix  <- 501:1000
## test <- A[test.ix,]
## p <- eval.online(A,N,train.ix,test.ix,fit,ties.method="first")
## q <- get.pred(train,A,test.ix,fit)

## test_that("train and test multinomial likelihoods agree with online version", {
##   llkm.train <- log(multinomial.score(q$m$train,train))
##   expect_that(sum(p$mllk$train), equals(sum(llkm.train)))
##   llkm.test  <- log(multinomial.score(q$m$test, test))
##   expect_that(sum(p$mllk$test), equals(sum(llkm.test)))
## })

## test_that("train and test ranks agree with online version", {
##   rk.train <- ranks(train,-q$lrm$train,ties.method="first")
##   rk.test  <- ranks(test, -q$lrm$test, ties.method="first")
##   expect_that(p$rks$train, equals(rk.train))
##   expect_that(p$rks$test, equals(rk.test))
## })


