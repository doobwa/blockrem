
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
  ans[1] <- 1 # intercept
  ans[2]  <- 1 # abba
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  expect_that(sij[[2]],equals(ans))
  ans <- rep(0,P)
  ans[1] <- 1 #intercept
  ans[7]  <- 1 # abab
  ans[8]  <- 1
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  ans[11] <- 1
  expect_that(sij[[3]],equals(ans))
})
test_that("get_v gets vectors as expected",{
  a <- x[[3]][[1]]
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
#   source("R/brem.r")
#   source("R/brem.cpp.r")
#   source("R/utils.r")
#   library(testthat)
#   require(abind)
#   set.seed(1)
  set.seed(1)
  M <- 100
  N <- 10
  times <- sort(runif(M,0,1))
  sen <- sample(1:N,M,replace=TRUE)
  rec <- sample(1:N,M,replace=TRUE)
  ix <- which(sen==rec)
  times <- times[-ix]
  sen <- sen[-ix]
  rec <- rec[-ix]
  M <- length(times)
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
  z <- c(rep(1,N/2),rep(2,N/2))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  
  lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  taus <- test_taus(lrm,times,sen-1,rec-1)
  taus2 <- test_taus_from_s(times,sen-1,rec-1,N,M,P)
  expect_that(all.equal(taus,taus2),is_true())
  
  llks <- llk_slow(lrm,times,sen-1,rec-1)
  llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  llk3 <- brem$llkfast(beta,z-1,s$ptr(),K)
  llk4 <- llk_fast(lrm,times,sen-1,rec-1)
  expect_that(llk3,equals(llk4))
})


test_stats_from_s <- function(times,sen,rec,N,M,P) {
  source("R/brem.cpp.r")
  s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  r <- brem$initializeStatistics(N,P);
  z <- z-1
  for (m in 0:(M-1)) {
    a <- sen[m+1]-1
    b <- rec[m+1]-1
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i != j) {
          x <- s$get_s(m,i,j)
          y <- r[,i+1,j+1]
          expect_that(x,equals(y))
        }
      }
    }
    r <- brem$updateStatistics(r,a,b,N,P)
  }
}