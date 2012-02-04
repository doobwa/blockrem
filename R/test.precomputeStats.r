library(testthat)
source("R/utils.r")
source("R/brem.cpp.r")
M <- 7
N <- 5
P <- 11
times <- c(1,2,3,4,5,6,7)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

s <- new(bremf$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()

x <- s$get_all_s()
v <- s$get_all_v()
s$ptr()

A <- cbind(times,sen,rec)
indx <- get.indices(A,N)
s$get_all_u()

ms <- c(0,2,4,6,7,8)


s$get_prev_index(ms,6)
s$get_prev_index(ms,5)
s$get_prev_index(ms,4)
s$get_prev_index(ms,3)
s$get_prev_index(ms,2)
s$get_prev_index(ms,1)
expect_that(s$get_prev_index(ms,6),equals(2))
expect_that(s$get_prev_index(ms,5),equals(2))
expect_that(s$get_prev_index(ms,4),equals(1))
expect_that(s$get_prev_index(ms,3),equals(1))
expect_that(s$get_prev_index(ms,2),equals(0))
expect_that(s$get_prev_index(ms,1),equals(0))
expect_that(s$get_prev_index(ms,0),equals(0))
s$get_s(2,3-1,1-1)
x[[3]][[1]]


ms <- c(0,2,4,6)
expect_that(s$get_prev(ms,6),equals(4))
expect_that(s$get_prev(ms,5),equals(4))
expect_that(s$get_prev(ms,4),equals(2))
expect_that(s$get_prev(ms,3),equals(2))
expect_that(s$get_prev(ms,2),equals(0))
expect_that(s$get_prev(ms,1),equals(0))
expect_that(s$get_prev(ms,0),equals(0))

times <- seq(.1,.7,by=.1)
expect_that(s$get_tau(6,2,0),equals())

# test_that("Correct number of statistics vectors",{
#   for (i in 1:N) {
#     for (j in 1:N) {
#       a <- length(tau[[i]][[j]])
#       b <- length(x[[i]][[j]])
#       expect_that(a,equals(b))
#     }
#   }
# })

ans <- rep(0,P)
ans[2]  <- 1 # abba
ans[9]  <- 1 # receiver out-degree
ans[10] <- 1 # sender in-degree
expect_that(s$get_s(2,3-1,1-1),equals(ans))

# (3,1) changes again when (5,1)
ans <- rep(0,P)
ans[7]  <- 1 # abab
ans[8]  <- 1
ans[9]  <- 2 # receiver out-degree
ans[10] <- 2 # sender in-degree
ans[11] <- 1
expect_that(s$get_s(4,3,1),equals(ans))
x[[3]][[1]]

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