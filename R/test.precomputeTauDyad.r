library(testthat)
source("R/utils.r")
source("R/brem.cpp.r")
M <- 7
N <- 5
times <- c(1,2,3,4,5,6,7)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

tau <- brem$precomputeTauDyad(times,sen-1,rec-1,N,M)
tau

expect_that(tau[[2]][[4]],equals(c(0,4,6)))
expect_that(tau[[3]][[1]],equals(c(0,0,1,2,3,5)))

ms <- c(0,2,4,6)
expect_that(brem$getTau(ms,6),equals(4))
expect_that(brem$getTau(ms,5),equals(4))
expect_that(brem$getTau(ms,4),equals(2))
expect_that(brem$getTau(ms,3),equals(2))
expect_that(brem$getTau(ms,2),equals(0))
expect_that(brem$getTau(ms,1),equals(0))
expect_that(brem$getTau(ms,0),equals(0))
