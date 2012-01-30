#source("R/brem.r")
source("R/brem.cpp.r")
library(testthat)

N <- 4
s <- list("abba" = matrix(0,N,N),
          "abby" = matrix(0,N,N),
          "abxa" = matrix(0,N,N),
          "abxb" = matrix(0,N,N),
          "abay" = matrix(0,N,N),
          "abab" = matrix(0,N,N),
          "sod" = matrix(0,N,N),
          "rod" = matrix(0,N,N),
          "sid" = matrix(0,N,N),
          "rid" = matrix(0,N,N))

i <- 1
j <- 2
a <- 2
b <- 1
s <- brem$updateStatistics(s,a-1,b-1,N)
times <- c(1.1,2,3,4)
sen <- c(1,3,3,1)
rec <- c(3,1,1,3)
M <- 4
beta <- list("intercept"=-1,
             "abba" = 3,
             "abby"=0,
             "abxa"=0,
             "abxb"=0,
             "abay"=0,
             "abab"=0,
             "sod"=0,
             "rod"=0,
             "sid"=0,
             "rid"=1)
lrm <- brem$lrm(beta,times,sen-1,rec-1,N,M)
llk2 <- brem$llk2(lrm,times,sen-1,rec-1,N,M)

test_that("computeLambda correct for small example",{
  
  beta <- list("intercept"=-1,
               "abba" = 3,
               "abby"=0,
               "abxa"=0,
               "abxb"=0,
               "abay"=0,
               "abab"=0,
               "sod"=0,
               "rod"=0,
               "sid"=0,
               "rid"=1)
  tmp <- matrix(0,N,N)
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      tmp[i+1,j+1] <- brem$computeLambda(i,j,s,beta)
    }
  }
  ans <- matrix(beta$intercept,N,N)
  ans[-b,b] <- ans[-b,b] + beta$rid  # rec. indegree effect
  ans[b,a] <- ans[b,a] + beta$abba
  expect_that(ans,equals(tmp))
})
test_that("statistics creation",{
  i <- 1
  j <- 2
  expect_that(s$abba[i,j],equals(1))
  expect_that(s$abby[i,j],equals(0))
  expect_that(s$sid[i,j],equals(1))
  expect_that(s$sod[i,j],equals(0))
  expect_that(s$rod[i,j],equals(1))
})

test_that("lrm and llk functions work on small example",{
  # Set up example
  set.seed(1)
  M <- 4
  N <- 4
  P <- 7
  times <- c(1,2,3,4)
  sen <- c(1,3,3,1)
  rec <- c(3,1,1,3)
  beta <- list("intercept"=1,
               "abba" = 1,
               "abby"=0,
               "abxa"=0,
               "abxb"=0,
               "abay"=0,
               "abab"=0,
               "sod"=0,
               "rod"=0,
               "sid"=0,
               "rid"=0)
  ix <- 1:N
  px <- rep(1,7)
  
  # Constract log rate matrix by hand and compare to drem$lrm
  a <- array(1,c(M,N,N))
  a[1,,] <- matrix(0,N,N)
  a[2,3,1] <- 2
  a[3,1,3] <- 2
  a[4,1,3] <- 2
  lrm <- brem$lrm(beta,times,sen-1,rec-1,N,M)
  expect_that(lrm,equals(a))
  
  # Compute log likelihood by hand.  
  diag(a[1,,]) <- diag(a[2,,]) <- diag(a[3,,]) <- diag(a[4,,]) <- -Inf
  llks <- c(a[1,sen[1],rec[1]],
            a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
            a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
            a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])) )
  sum(llks)
  
  # Compare to drem$llk2
  llk2 <- brem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
#   # Compare to drem$llk
#   px <- c(1,1,0,0,0,0,0)
#   llk3 <- drem$llk(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M)
#   expect_that(sum(llks),equals(llk3))
  
  # Compare to drem$allk
#   px <- c(1,1,0,0,0,0,0)  
#   allk <- drem$allk(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M,N)
#   expect_that(sum(llks),equals(allk))
})
# 
# test_that("log rate matrix computation correct",{
#   # Make sure example is in the first group
#   sen[1] <- 1
#   rec[1] <- 3
#   lrm <- drem$lrm(beta,times,sen,rec,ix,ix,px,N,M)
#   expect_that(dim(lrm),equals(c(100,10,10)))
#   expect_that(lrm[2,5,1],equals(beta[1]))
#   expect_that(lrm[2,rec[1]+1,sen[1]+1],equals(beta[1] + beta[2]))  # AB-BA
#   
#   # diagnoals should be -15
#   #expect_that(all(diag(lrm[5,,]) == -15),is_true())
#   
#   # only the (ix,ix) portion should be populated
#   ix <- 1:5
#   lrm <- drem$lrm(beta,times,sen,rec,ix-1,ix-1,px,N,M)
#   expect_that(sum(lrm[,-ix,-ix]),equals(0))  
# })

