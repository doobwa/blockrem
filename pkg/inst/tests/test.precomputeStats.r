context("using precomputed statistics")

M <- 8
N <- 5
P <- 15
times <- seq(0,.6,by=.1)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)
ego <- P
s <- new(RemStat,times,sen-1,rec-1,N,M,ego)
s$precompute()
x <- s$get_all_s()

test_that("Can use precomputed data structures",{
  
  expect_that(length(x),equals(5))
  
  sij <- x[[3]][[1]]
  expect_that(length(sij),equals(7))
  expect_that(length(sij[[1]]),equals(P))
  
  v <- s$get_all_v()
  
  expect_that( s$ptr(), is_a("externalptr") )
})

test_that("a few statistics vectors are correct",{
  sij <- x[[3]][[1]]
  expect_that(sij[[1]],equals(rep(0,P)))
  ans <- rep(0,P)
  ans[1] <- 1 # intercept
  ans[2]  <- 1 # abba
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  ans[13] <- 0 # event index of the starting changepoint
  ans[14] <- 1 # rrs: (1,3) occurred last changepoint
#  ans[15] <- 0 # rss
  expect_that(sij[[2]],equals(ans))
  ans <- rep(0,P)
  ans[1] <- 1 #intercept
  ans[7]  <- 1 # abab
  ans[8]  <- 1
  ans[9]  <- 1 # receiver out-degree
  ans[10] <- 1 # sender in-degree
  ans[11] <- 1 # receiver in degree
  ans[12] <- 1 # dyad count
  ans[13] <- 1 # event count
  ans[14] <- 1 # rrs: (1,3) only incoming event for a=3
  ans[15] <- 1 # rss: (3,1) occurred last changepoint
  expect_that(sij[[3]],equals(ans))
  sij <- x[[4]][[2]]
  ans <- rep(0,P)
  ans[1] <- 1 #intercept
  ans[4] <- 1 # ab-xa
  ans[9]  <- 1 # receiver out-degree
  ans[13] <- 4 # global event count
#  ans[14] <- 
  expect_that(sij[[2]],equals(ans))
})

test_that("get_v gets vectors as expected",{
  a <- x[[3]][[1]]
  b <- s$get_v(3-1,1-1)
  d <- s$get_w(3-1,1-1)
  expect_that(b,equals(c(0,0,1,2,3,5,6)))
  expect_that(d,equals(c(1,2,3,4,4,5,6)))
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



test_stats_from_s <- function(times,sen,rec,N,M,P) {
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
  s <- new(RemStat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  r <- InitializeStatisticsArray(N,P)
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
    r <- UpdateStatisticsArray(r,a,b,N,P)
  }
}

test_that("Transform works as expected",{
  s0 <- s1 <- s$get_all_s()[[5]][[3]][[7]]
  s1[8:12] <- log((s0[8:12]+1) / (s0[13] + N*(N-1)))
  s$transform()
  s2 <- s$get_all_s()[[5]][[3]][[7]]
  expect_that(s1,equals(s2))
})

test_that("Transform works for recency effects", {

})

test_that("Ego restriction updates fewer dyads", {
  ego <- 1
  s <- new(RemStat,times,sen-1,rec-1,N,M,ego)
  s$precompute()
  x <- s$get_all_s()
  
  sij <- do.call(rbind,x[[5]][[2]])
  ans <- c(0,4,5,7)
  expect_that(sij[,13], equals(ans))

  sij <- do.call(rbind,x[[4]][[1]])
  ans <- c(0,6,7)
  expect_that(sij[,13], equals(ans))
})

test_that("LogLambdaPc uses degree effects properly",{
  s <- new(RemStat,times,sen-1,rec-1,N,M,ego)
  s$precompute()
  i <- 1
  j <- 3
  zi <- 1
  zj <- 1
  K <- 2
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
               "rid"=matrix(0,K,K),
               "dc"=matrix(0,K,K),
               "cc"=matrix(0,K,K),
               "rrs"=matrix(0,K,K),
               "rss"=matrix(0,K,K))
  beta <- abind(beta,rev.along=3)
  # Add in some degree effects
  beta[8,1,1] <- 4
  beta[9,1,1] <- 3
  beta[10,1,1] <- 2
  beta[11,1,1] <- 1
  beta[12,1,1] <- .5
  
  sij <- s$get_s(3,i-1,j-1)
  sij[8:12] <- log((sij[8:12]+1)/(sij[13] + N*(N-1)))
  sij <- sij[-13]
  b <- beta[-13,zi,zj]
  ans <- as.numeric(b %*% sij)

  s$transform()
  sij <- s$get_s(3,i-1,j-1)
  lam <- LogLambdaPc(i-1,j-1,zi-1,zj-1,sij,beta,N,K,P)
  
  expect_that(lam,equals(ans))
})

