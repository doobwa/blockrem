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
             "dc"=matrix(0,K,K),
             "cc"=matrix(0,K,K))
z <- c(rep(1,N/2),rep(2,N/2))
P <- length(beta)
beta <- abind(beta,rev.along=3)
set.seed(1)
sim <- generate.brem(M,N,beta,z)
A <- sim$edgelist

train <- A[1:500,]
test <- A
train.ix <- 1:500
test.ix  <- 501:1000
test <- A[test.ix,]

## Precompute data structures
ego <- 0
strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P,ego)
strain$precompute()
stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P,ego)
stest$precompute()

fit <- list(params=list(beta=beta,z=z),ego=ego,beta=beta,z=z)

## Precompute rate arrays
lrm <- list()
lrm$train <- brem.lrm.fast(strain, fit$params$z, fit$params$beta)
lrm$test  <- brem.lrm.fast(stest, fit$params$z, fit$params$beta)

test_that("train and test BREM likelihood agree with online version", {
  llk.train <- RemLogLikelihoodVecFromArray(lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
  llk.test <- RemLogLikelihoodVecFromArray(lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

  for (i in 2:10) {
    expect_that(eval.llk(train,N,fit,i,strain), equals(llk.train[i]))
  }
  
})

test_that("Helper functions for multinomial likelihood and ranks work", {
  N <- 10
  lrm <- matrix(1:100,N,N)
  train <- cbind(0:1,1:2,2:1)
  diag(lrm) <- -Inf
  m <- exp(lrm)/sum(exp(lrm))
  i <- 1
  ans <- log(m[train[i,2],train[i,3]])
  expect_that(eval.mult(train,lrm,i), equals(ans))
  i <- 1
  ans <- 81
  expect_that(eval.rank(train,lrm,i,ties.method="first"), equals(ans))
  i <- 2
  ans <- 90
  expect_that(eval.rank(train,lrm,i,ties.method="first"), equals(ans))
})


test_that("BREM models: online mult likelihoods and ranks agree with old code", {
  p <- evaluate(A,N,train.ix,test.ix,fit,ties.method="first")
  q <- get.pred(train,A,N,test.ix,fit)
  llkm.train <- log(multinomial.score(q$m$train,train))
  llkm.test  <- log(multinomial.score(q$m$test, test))
  rk.train <- ranks(train,-q$lrm$train,ties.method="first")
  rk.test  <- ranks(test, -q$lrm$test, ties.method="first")
  expect_that(sum(p$mllk$train), equals(sum(llkm.train)))
  expect_that(sum(p$mllk$test), equals(sum(llkm.test)))
  expect_that(p$rks$train, equals(rk.train))
  expect_that(p$rks$test, equals(rk.test))
})

test_that("BREM models: online mult likelihoods and ranks check with hard coded answers", {
  p <- evaluate(A,N,train.ix,test.ix,fit,ties.method="first")
  ans <- c(6, 25, 77, 74, 3, 10)
  ans <- c(74, 34, 11, 18, 73, 68)
  ans <- c(-4.49980967033027, -4.81081587105776, -5.02537295559027, -4.92589230086563, -3.2896979969751, -4.35176341317661)
  expect_that(head(p$mllk$train), equals(ans))
  ans <- c(-5.04026531620789, -5.09140599388328, -4.15951847928562, -4.13957979834807, -5.13673903647754, -5.2896979969751)
  expect_that(head(p$mllk$test), equals(ans))
  ans <- c(-1, -1.1397952618685, -1.32159558497369, -3.89496853746407, 
0.460317160002888, -0.435237016528844)
  expect_that(head(p$llk$train), equals(ans))
  ans <- c(-2.18733996986865, -3.52452590117665, -0.508645472116811, -0.0568255931138519, -1.08559791999357, -3.1297706427309)
  expect_that(head(p$llk$test), equals(ans))
})

test_that("Baseline models: online mult likelihoods and ranks agree with old", {
  p <- evaluate.baseline(A,N,train.ix,test.ix,model="online",ties.method="first")
  q <- get.pred.baseline(train,A,N,test.ix,model="online")
  llkm.train <- log(multinomial.score(q$m$train,train))
  llkm.test  <- log(multinomial.score(q$m$test, test))
  rk.train <- ranks(train,-q$lrm$train,ties.method="first")
  rk.test  <- ranks(test, -q$lrm$test, ties.method="first")
  expect_that(sum(p$mllk$train), equals(sum(llkm.train)))
  expect_that(sum(p$mllk$test), equals(sum(llkm.test)))
  expect_that(p$rks$train, equals(rk.train))
  expect_that(p$rks$test, equals(rk.test))
})

test_that("Baseline models: online mult likelihoods and ranks check with hard coded answers", {
  p <- evaluate.baseline(A,N,train.ix,test.ix,model="online",ties.method="first")
  ans <- c(6,17,75,73,28,38)
  expect_that(head(p$rks$train), equals(ans))
  ans <- c(58,47,80,41,38,15)
  expect_that(head(p$rks$test), equals(ans))
  ans <- c(-4.49980967033027, -4.51085950651685, -4.52178857704904, -4.53259949315326, -4.54329478227, -4.55387689160054)
  expect_that(head(p$mllk$train), equals(ans))
  ans <- c(-4.58836306767171, -4.59005654817804, -4.99721227376411, -4.43928424994241, -4.44096917030733, -4.19133682820941)
  expect_that(head(p$mllk$test), equals(ans))
  # TODO: llk check for evaluate.baseline
})

test_that("BREM likelihood same between Pc method and evaluate()", {

  llk.pc  <- RemLogLikelihoodPc(fit$beta,fit$z-1,strain$ptr(),K)
  llk.train <- RemLogLikelihoodVecFromArray(lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))

  res <- sapply(1:nrow(train),function(i) {
    eval.brem(train,lrm$train[i,,],i)
  })
  expect_that(sum(res), equals(sum(llk.pc)))
  expect_that(res, equals(llk.train))

  llk.test <- RemLogLikelihoodVecFromArray(lrm$test,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A))[test.ix]
  res <- sapply(1:1000,function(i) {
    eval.brem(A,lrm$test[i,,],i)
  })
  expect_that(res[test.ix], equals(llk.test))
  
})
