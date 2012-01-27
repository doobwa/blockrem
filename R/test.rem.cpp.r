source("R/rem.cpp.r")
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
})
test_that("llk works on small example",{
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
  beta <- rnorm(7)
  drem$computeLambda(3,2,2,3,beta,rep(1L,7))
})

test_that("sampling without replacement examples are correct",{
  x <- drem$sampleWithoutReplacement(1:20,10)
  expect_that(length(x),equals(10))
  expect_that(max(x) <= 20,is_true())
  
  x <- drem$sampleWithoutReplacement(6:10,5)
  expect_that(x,equals(6:10))
  
  x <- drem$sampleWithoutReplacement(50:100,5)
  expect_that(all(x <=100 & x>=50),is_true())
})

test_that("approximate likelihood correct",{
  N <- 100
  ix <- 1:N - 1
  sen <- sample(ix,M,replace=TRUE)
  rec <- sample(ix,M,replace=TRUE)
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M,P)
  allk <- drem$allk(beta,times,sen,rec,ix,ix,px,N,M,P,N)
  expect_that(llk,equals(allk))
  
  allks <- sapply((N/2):N,function(i) {
    drem$allk(beta,times,sen,rec,ix,ix,px,N,M,P,i)
  })
  plot(allks,type="l")
})