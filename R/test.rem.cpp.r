source("R/brem.r")
source("R/rem.cpp.r")
library(testthat)

test_that("lrm and llk functions work on small example",{
  # Set up example
  set.seed(1)
  M <- 4
  N <- 3
  P <- 7
  times <- c(1,2,3,4)
  sen <- c(1,3,3,1)
  rec <- c(3,1,1,3)
  beta <- c(1,1,0,0,0,0,0)
  ix <- 1:N
  px <- rep(1,7)
  
  # Constract log rate matrix by hand and compare to drem$lrm
  a <- array(1,c(M,N,N))
  a[1,,] <- matrix(0,3,3)
  a[2,3,1] <- 2
  a[3,1,3] <- 2
  a[4,1,3] <- 2
  lrm <- drem$lrm(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M)
  expect_that(lrm,equals(a))
  
  # Compute log likelihood by hand.  
  diag(a[1,,]) <- diag(a[2,,]) <- diag(a[3,,]) <- diag(a[4,,]) <- -Inf
  llks <- c(a[1,sen[1],rec[1]],
            a[2,sen[2],rec[2]] - (times[2]-times[2-1]) * sum(exp(a[2,,])),
            a[3,sen[3],rec[3]] - (times[3]-times[3-1]) * sum(exp(a[3,,])),
            a[4,sen[4],rec[4]] - (times[4]-times[4-1]) * sum(exp(a[4,,])) )
  sum(llks)
  
  # Compare to drem$llk2
  llk2 <- drem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  # Compare to drem$llk
  px <- c(1,1,0,0,0,0,0)
  llk3 <- drem$llk(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M)
  expect_that(sum(llks),equals(llk3))
  
  # Compare to drem$allk
  px <- c(1,1,0,0,0,0,0)  
  allk <- drem$allk(beta,times,sen-1,rec-1,ix-1,ix-1,px,N,M,N)
  expect_that(sum(llks),equals(allk))
})

N <- 10
M <- 100
ix <- 1:5 - 1
jx <- 5:10 - 1
times <- sort(runif(M,0,10))
sen <- sample(ix,M,replace=TRUE)
rec <- sample(jx,M,replace=TRUE)
px <- rep(1,7)
beta <- rnorm(7)

test_that("llk on subset of dyads runs on small example",{
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M)
})

test_that("log rate matrix computation correct",{
  # Make sure example is in the first group
  sen[1] <- 1
  rec[1] <- 3
  lrm <- drem$lrm(beta,times,sen,rec,ix,ix,px,N,M)
  expect_that(dim(lrm),equals(c(100,10,10)))
  expect_that(lrm[2,5,1],equals(beta[1]))
  expect_that(lrm[2,rec[1]+1,sen[1]+1],equals(beta[1] + beta[2]))  # AB-BA
  
  # diagnoals should be -15
  #expect_that(all(diag(lrm[5,,]) == -15),is_true())
  
  # only the (ix,ix) portion should be populated
  ix <- 1:5
  lrm <- drem$lrm(beta,times,sen,rec,ix-1,ix-1,px,N,M)
  expect_that(sum(lrm[,-ix,-ix]),equals(0))  
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

# test_that("get speedup over brute force llk computation",{
#   N <- 50
#   M <- 1000
#   z <- rep(1,N)
#   px <- c(1,1,1,1,1,1,1)
#   beta <- array(rnorm(length(px)),c(1,1,length(px)))
#   sim <- simulate.brem(M,N,z,beta,px)
#   ix <- 1:N - 1
#   llk <- drem$llk(beta,times,sen-1,rec-1,ix,ix,px,N,M)
#   allk <- drem$allk(beta,times,sen-1,rec-1,ix,ix,px,N,M,N)
#   expect_that(llk,equals(allk))
# })

test_that("approximate likelihood correct",{
  N <- 200
  M <- 100
  ix <- 1:N - 1
  sen <- sample(ix,M,replace=TRUE)
  rec <- sample(ix,M,replace=TRUE)
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M)
  allk <- drem$allk(beta,times,sen,rec,ix,ix,px,N,M,N)
  expect_that(llk,equals(allk))
  
  allks <- sapply(2:N,function(i) {
    drem$allk(beta,times,sen,rec,ix,ix,px,N,M,i)
  })
  plot(allks[50:100],type="l")
  abline(h=llk)
})
test_that("Check scaling properities of app. llk",{
  M <- 5000
  N <- 1000
  ix <- 1:N - 1
  times <- sort(runif(M,0,100))
  sen <- sample(ix,M,replace=TRUE)
  rec <- sample(ix,M,replace=TRUE)
  
  system.time(drem$llk(beta,times,sen,rec,ix,ix,px,N,M))
  system.time(drem$allk(beta,times,sen,rec,ix,ix,px,N,M,250))
  system.time(drem$allk(beta,times,sen,rec,ix,ix,px,N,M,50))
  
  llk <- drem$llk(beta,times,sen,rec,ix,ix,px,N,M)
  allks <- sapply(seq(50,N,by=100),function(i) {
    print(i)
    drem$allk(beta,times,sen,rec,ix,ix,px,N,M,i)
  })
  plot(allks,type="l")
  abline(h=llk,col="red")
})
  