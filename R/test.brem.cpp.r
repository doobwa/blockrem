source("brem.r")
source("brem.cpp.r")
library(testthat)

require(abind)

# Example 1
N <- 4
s <- list("intercept" = matrix(1,N,N),
          "abba" = matrix(0,N,N),
          "abby" = matrix(0,N,N),
          "abxa" = matrix(0,N,N),
          "abxb" = matrix(0,N,N),
          "abay" = matrix(0,N,N),
          "abab" = matrix(0,N,N),
          "sod" = matrix(0,N,N),
          "rod" = matrix(0,N,N),
          "sid" = matrix(0,N,N),
          "rid" = matrix(0,N,N))
P <- length(s)
s <- abind(s,rev.along=3)
i <- 1
j <- 2
a <- 2
b <- 1
s <- brem$updateStatistics(s,a-1,b-1,N,P)

test_that("statistics creation works",{
  i <- 1
  j <- 2
  dimnames(s) <- list(NULL,NULL,NULL)
  expect_that(s[2,i,j],equals(1))  # abba
  expect_that(s[3,i,j],equals(0)) # abby
  expect_that(s[10,i,j],equals(1)) # sid
  expect_that(s[8,i,j],equals(0)) # sod
  expect_that(s[9,i,j],equals(1)) # rod
})


test_that("computeLambda correct for small example",{
  times <- c(1,2,3,4)
  sen <- c(1,3,3,1)
  rec <- c(3,1,1,3)
  M <- 4
  N <- 4
  K <- 1
  z <- rep(1,N)
  beta <- list("intercept"=matrix(-1,1,1),
               "abba" = matrix(3,1,1),
               "abby"=matrix(0,1,1),
               "abxa"=matrix(0,1,1),
               "abxb"=matrix(0,1,1),
               "abay"=matrix(0,1,1),
               "abab"=matrix(0,1,1),
               "sod"=matrix(0,1,1),
               "rod"=matrix(0,1,1),
               "sid"=matrix(0,1,1),
               "rid"=matrix(1,1,1))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  
  tmp <- matrix(0,N,N)
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      tmp[i+1,j+1] <- brem$computeLambda(i,j,0,0,s,beta,N,K,P)
    }
  }
  ans <- matrix(beta[1,1,1],N,N)
  ans[-b,b] <- ans[-b,b] + beta[11,1,1]  # rec. indegree effect
  ans[b,a] <- ans[b,a] + beta[2,1,1]
  expect_that(ans,equals(tmp))
})
test_that("lrm and llk functions work on small example for K=1",{
  # Set up example
  set.seed(1)
  M <- 4
  N <- 5
  times <- c(1,2,3,4)
  sen <- c(1,3,3,1)
  rec <- c(3,1,1,3)
  K <- 1
  beta <- list("intercept"=matrix(1,1,1),
               "abba" = matrix(1,1,1),
               "abby"=matrix(0,1,1),
               "abxa"=matrix(0,1,1),
               "abxb"=matrix(0,1,1),
               "abay"=matrix(0,1,1),
               "abab"=matrix(0,1,1),
               "sod"=matrix(0,1,1),
               "rod"=matrix(0,1,1),
               "sid"=matrix(0,1,1),
               "rid"=matrix(0,1,1))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  
  a <- 1
  b <- 3
  s <- array(0,c(P,N,N))
  s <- brem$updateStatistics(s,a-1,b-1,N,P)
  brem$computeLambda(1,0,0,0,s,beta,N,K,P)  # 
  
  # Constract log rate matrix by hand and compare to drem$lrm
  a <- array(1,c(M,N,N))
  a[1,,] <- matrix(0,N,N)
  a[2,3,1] <- 2
  a[3,1,3] <- 2
  a[4,1,3] <- 2
  z <- rep(1,N)
  K <- 1
  lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  expect_that(lrm,equals(a))
  
  # Compute log likelihood by hand.  
  diag(a[1,,]) <- diag(a[2,,]) <- diag(a[3,,]) <- diag(a[4,,]) <- -Inf
  llks <- c(a[1,sen[1],rec[1]],
            a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
            a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
            a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])) )
  sum(llks)
  
  # Compare to drem$llk2
  llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  llk3 <- brem.llk.slow(lrm,times,sen,rec,z,N,M)
  expect_that(sum(llks),equals(llk3))
  
  # Compare to drem$llk
  system.time(llk4 <- brem$llk(beta,times,sen-1,rec-1,z-1,N,M,K,P))
  browser()
  expect_that(sum(llks),equals(llk3))
  
  # Compare to drem$allk
#   px <- c(1,1,0,0,0,0,0)  
#   allk <- drem$allk(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M,N)
#   expect_that(sum(llks),equals(allk))
})

test_that("llk runs for larger example",{
  
  source("R/brem.cpp.r")
  require(abind)
  set.seed(1)
  M <- 1000
  N <- 100
  times <- sort(runif(M,0,1))
  sen <- sample(1:N,M,replace=TRUE)
  rec <- sample(1:N,M,replace=TRUE)
  K <- 1
  beta <- list("intercept"=matrix(1,1,1),
               "abba" = matrix(1,1,1),
               "abby"=matrix(0,1,1),
               "abxa"=matrix(0,1,1),
               "abxb"=matrix(0,1,1),
               "abay"=matrix(0,1,1),
               "abab"=matrix(0,1,1),
               "sod"=matrix(0,1,1),
               "rod"=matrix(0,1,1),
               "sid"=matrix(0,1,1),
               "rid"=matrix(0,1,1))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  z <- rep(1,N)
  system.time(llk1 <- brem$llk(beta,times,sen-1,rec-1,z-1,N,M,K,P))
  system.time(llk2 <- brem$llkp(beta,times,sen-1,rec-1,z-1,N,M,K,P))
  
})

test_that("lrm and llk functions work on small example for K=2",{
  # Set up example
  set.seed(1)
  M <- 5
  N <- 6
  K <- 2
  times <- c(1,2,3,4,5)
  sen <- c(1,5,1,5,1)
  rec <- c(3,1,5,6,3)
  
  beta <- list("intercept"=matrix(1,K,K),
               "abba" = matrix(c(1,2,3,4),K,K),
               "abby"=matrix(0,K,K),
               "abxa"=matrix(0,K,K),
               "abxb"=matrix(0,K,K),
               "abay"=matrix(0,K,K),
               "abab"=matrix(0,K,K),
               "sod"=matrix(0,K,K),
               "rod"=matrix(0,K,K),
               "sid"=matrix(0,K,K),
               "rid"=matrix(0,K,K))
  z <- c(1,1,1,2,2,2)
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  
  # Constract log rate matrix by hand and compare to drem$lrm
  a <- array(1,c(M,N,N))
  a[1,,] <- matrix(0,N,N)
  a[2,3,1] <- 1 + 1
  a[4,5,1] <- 1 + 2
  a[3,1,5] <- 1 + 3
  a[5,6,5] <- 1 + 4
  lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  expect_that(lrm,equals(a))
  
  # Compute log likelihood by hand.  
  diag(a[1,,]) <- diag(a[2,,]) <- diag(a[3,,]) <- diag(a[4,,]) <- diag(a[5,,]) <-  -Inf
  llks <- c(a[1,sen[1],rec[1]],
            a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
            a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
            a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])),
            a[5,sen[5],rec[5]] - (times[5]-times[5-1]) * sum(exp(a[5,,])))
  sum(llks)
  
  # Compare to drem$llk2
  llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  # Test R interface
  A <- cbind(times,sen,rec)
  expect_that(sum(llks),equals(brem.llk(A,N,z,beta)))
  
  # make sure lpost runs
  brem.lpost(A,N,K,z,beta)
})