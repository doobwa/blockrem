source("../brem.r")
source("../brem.cpp.r")
source("../utils.r")

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
  times <- c(0,2,3,4)
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
  times <- c(0,2,3,4)
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
  K <- 1
  lrm <- brem$lrm(beta,times,sen-1,rec-1,z-1,N,M,K,P)
#   expect_that(lrm,equals(a))
  
  # Compute log likelihood by hand.  
  llks <- llk_slow(lrm,times,sen-1,rec-1)
  
  # Compare to drem$llk2
  llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)
  
  s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  llk3 <- brem$llkfast(beta,z-1,s$ptr(),K)
  lrm2 <- lrm_slow(beta,z-1,s,M,N,K,P)

  
  true.fast <- c(1,
                 2 - 2*(13*exp(1) + exp(2)), 
                 1 - (13*exp(1) + exp(2)), 
                 2 -  (times[4] - times[1])*6*exp(1))

  llk4 <- llk_fast(lrm,times,sen-1,rec-1)

  expect_that(sum(llks),equals(llk2))
  expect_that(sum(llks),equals(sum(llk3)))
  expect_that(sum(llks),equals(sum(llk4)))
  expect_that(sum(true.fast),equals(sum(llk4)))
})

test_that("lrm and llk functions work on small example for K=2",{
  # Set up example
  set.seed(1)
  M <- 1000
  N <- 100
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
  llks <- llk_slow(lrm,times,sen-1,rec-1)
  
  # Compare to drem$llk2
  llk2 <-  brem$llk2(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  llk3 <- brem$llkfast(beta,z-1,s$ptr(),K)
  
  llk4 <- llk_fast(lrm,times,sen-1,rec-1)
  expect_that(sum(llk3),equals(sum(llk4)))
  
  expect_that(sum(llks) + 1,equals(sum(llk4)))  # Off by intercept term bug
  
#   x <- llk_fast_last(lrm,times,sen-1,rec-1)
#   y <- brem$test_last(beta,z-1,s$ptr(),K)
#   expect_that(all.equal(x$taus,y$taus),is_true())
})