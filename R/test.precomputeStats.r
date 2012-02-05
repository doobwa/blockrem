
source("R/brem.cpp.r")
library(testthat)
source("R/utils.r")
M <- 7
N <- 5
P <- 11
times <- seq(0,.6,by=.1)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()

x <- s$get_all_s()
v <- s$get_all_v()
s$ptr()

test_that("a few statistics vectors are correct",{
  sij <- s$get_all_s()[[3]][[1]]
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

times <- seq(0,.6,by=.1)
taus <- sapply(0:(M-1),function(m) s$get_tau(m,3-1,1-1))
expect_that(s$get_tau(1,3-1,1-1),equals(.5))
expect_that(s$get_tau(M-1,3-1,1-1),equals(.5))
