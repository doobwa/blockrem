library(testthat)
source("R/utils.r")
source("R/brem.cpp.r")
M <- 7
N <- 5
P <- 11
times <- c(1,2,3,4,5,6,7)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

tau <- bremf$precomputeTauDyad(times,sen-1,rec-1,N,M)
s <- bremf$precomputeStats(times,sen-1,rec-1,N,M,P)

test_that("Correct number of statistics vectors",{
  for (i in 1:N) {
    for (j in 1:N) {
      a <- length(tau[[i]][[j]])
      b <- length(s[[i]][[j]])
      expect_that(a,equals(b))
    }
  }
})

ans <- rep(0,P)
ans[2]  <- 1 # abba
ans[9]  <- 1 # receiver out-degree
ans[10] <- 1 # sender in-degree
expect_that(s[[3]][[1]][[2]],equals(ans))

# (3,1) changes again when (5,1)
ans <- rep(0,P)
ans[7]  <- 1 # abab
ans[8]  <- 1
ans[9]  <- 2 # receiver out-degree
ans[10] <- 2 # sender in-degree
ans[11] <- 1
expect_that(s[[3]][[1]][[4]],equals(ans))
s[[3]][[1]]

# Time test
set.seed(1)
M <- 3000
N <- 300
P <- 11
times <- sort(runif(M,0,1))
sen <- sample(1:N,M,replace=TRUE)
rec <- sample(1:N,M,replace=TRUE)
s <- bremf$precomputeStats(times,sen-1,rec-1,N,M,P)
system.time(bremf$precomputeStats(times,sen-1,rec-1,N,M,P))
lsos()