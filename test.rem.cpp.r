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
  system.time(llk <- drem$llk(beta,times,sen,rec,ix,ix,N,M,P))
  
  N <- 10
  ix <- 1:5
  jx <- 5:10
  sen <- sample(ix,M,replace=TRUE) - 1
  rec <- sample(jx,M,replace=TRUE) - 1
  llk <- drem$llk(beta,times,sen,rec,ix,jx,N,M,P)
})