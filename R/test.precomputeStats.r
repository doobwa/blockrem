
source("R/brem.cpp.r")
library(testthat)
source("R/utils.r")
library(abind)

M <- 7
N <- 5
P <- 11
times <- seq(0,.6,by=.1)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()

x <- s$get_all_s()
sij <- x[[3]][[1]]
v <- s$get_all_v()
s$ptr()

test_that("a few statistics vectors are correct",{
  expect_that(sij[[1]],equals(rep(0,11)))
  ans <- rep(0,P)
  ans[2]  <- 1 # abba
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  expect_that(sij[[2]],equals(ans))
  ans <- rep(0,P)
  ans[7]  <- 1 # abab
  ans[8]  <- 1
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  ans[11] <- 1
  expect_that(sij[[3]],equals(ans))
})
test_that("get_v gets vectors as expected",{
  a <- s$get_all_s()[[3]][[1]]
  b <- s$get_v(3-1,1-1)
  d <- s$get_w(3-1,1-1)
  expect_that(b,equals(c(0,1,2,3,5,6)))
  expect_that(d,equals(c(0,1,2,3,3,4,5)))
  expect_that(length(a), equals(length(b)))
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i!=j) {
        a <- s$get_v(i-1,j-1)
        b <- s$get_all_s()[[i]][[j]]
        expect_that(length(a),equals(length(b)))
      }
    }
  }
})

test_that("taus from get_tau match with R version",{
  
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
               "rid"=matrix(0,K,K))
  z <- c(1,1,1,2,2)#c(rep(1,N/2),rep(2,N/2))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  times <- seq(0,.6,by=.1)
  taus <- sapply(0:(M-1),function(m) s$get_tau(m,3-1,1-1))
  expect_that(s$get_tau(1,3-1,1-1),equals(0))
  expect_that(s$get_tau(M-1,3-1,1-1),equals(.5))
  
  lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  taus <- test_taus(lrm,times,sen-1,rec-1)
  taus2 <- test_taus_from_s(times,sen-1,rec-1,N,M,P)
  expect_that(all.equal(taus,taus2),is_true())
})
