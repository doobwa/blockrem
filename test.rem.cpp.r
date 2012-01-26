source("rem.cpp.r")
library(testthat)
test_that("llk works on different subsets of dyads",{
  
  set.seed(1)
  M <- 5000
  N <- 100
  P <- 6
  times <- sort(runif(M,0,100))
  sen <- sample(1:N,M,replace=TRUE) - 1
  rec <- sample(1:N,M,replace=TRUE) - 1
  beta <- rnorm(P)
  ix <- 0:(N-1)
  px <- rep(1,7)
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M,P)
  llk
  N <- 10
  ix <- 1:5 - 1
  jx <- 5:10 - 1
  sen <- sample(ix,M,replace=TRUE) - 1
  rec <- sample(jx,M,replace=TRUE) - 1
  llk <- drem$llk(beta,times,sen,rec,ix,jx,px,N,M,P)
})

test_that("log rate matrix computation works",{
  lrm <- drem$lrm(beta,times,sen,rec,ix,jx,px,N,M,P)
})

test_that("calculating lambda runs",{
  drem$testIntConv(c(3,2,3))
  beta <- rnorm(7)
  expect_that(drem$computeLambdaOld(3,2,2,3,beta),
              equals(drem$computeLambda(3,2,2,3,beta,rep(1L,7)))
})