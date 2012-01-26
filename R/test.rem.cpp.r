source("rem.cpp.r")
library(testthat)

set.seed(1)
M <- 100
N <- 10
P <- 6
times <- sort(runif(M,0,100))
sen <- sample(1:N,M,replace=TRUE) - 1
rec <- sample(1:N,M,replace=TRUE) - 1
beta <- rnorm(P)
ix <- 0:(N-1)
px <- rep(1,7)

test_that("llk works on different subsets of dyads",{
  
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M,P)
  llk
  N <- 10
  ix <- 1:5 - 1
  jx <- 5:10 - 1
  sen <- sample(ix,M,replace=TRUE)
  rec <- sample(jx,M,replace=TRUE)
  llk <- drem$llk(beta,times,sen,rec,ix,jx,px,N,M,P)
})

test_that("log rate matrix computation correct",{
  lrm <- drem$lrm(beta,times,sen,rec,ix,ix,px,N,M,P)
  expect_that(lrm[2,1,1],equals(beta[1]))
  expect_that(lrm[2,rec[1]+1,sen[1]+1],equals(beta[1] + beta[2]))  # AB-BA
})

test_that("calculating lambda runs",{
  drem$testIntConv(c(3,2,3))
  beta <- rnorm(7)
  expect_that(drem$computeLambdaOld(3,2,2,3,beta),
              equals(drem$computeLambda(3,2,2,3,beta,rep(1L,7)))
})