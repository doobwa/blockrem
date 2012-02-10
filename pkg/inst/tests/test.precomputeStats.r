context("using precomputed statistics")

M <- 7
N <- 5
P <- 11
times <- seq(0,.6,by=.1)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)
s <- new(Stat,times,sen-1,rec-1,N,M,P)
s$precompute()
x <- s$get_all_s()

test_that("Can use precomputed data structures",{
  
  expect_that(length(x),equals(5))
  
  sij <- x[[3]][[1]]
  expect_that(length(sij),equals(7))
  
  v <- s$get_all_v()
  
  expect_that( s$ptr(), is_a("externalptr") )
})

test_that("a few statistics vectors are correct",{

  sij <- x[[3]][[1]]
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

test_that("taus from get_tau match with R version",{
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
               "rid"=matrix(0,K,K))
  z <- c(rep(1,N/2),rep(2,N/2))
  P <- length(beta)
  beta <- abind(beta,rev.along=3)
  
  lrm <- log_intensity_array(beta,times,sen-1,rec-1,z-1,N,M,K,P)
  
  s <- new(Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  
  lrm2 <- lrm_slow(beta,z-1,s,M,N,K,P)
  
  taus <- test_taus(lrm,times,sen-1,rec-1,M,N)
  taus2 <- test_taus_from_s(times,sen-1,rec-1,N,M,P)
  expect_that(all.equal(taus,taus2),is_true())
  
  llks <- llk_slow(lrm,times,sen-1,rec-1,M,N)
  llk2 <-  loglikelihood_from_lrm(lrm,times,sen-1,rec-1,N,M)
  expect_that(sum(llks),equals(llk2))
  
  llk3 <- loglikelihood_fast(beta,z-1,s$ptr(),K)
  llk4 <- llk_fast(lrm,times,sen-1,rec-1,M,N)
  expect_that(llk3,equals(llk4))
  
  # First event of non-fast version has a mistake
  expect_that(sum(llk3),equals(llk2 + 1))
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
  s <- new(Stat,times,sen-1,rec-1,N,M,P)
  s$precompute()
  r <- initialize_statistics(N,P)
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
    r <- update_statistics(r,a,b,N,P)
  }
}

test_that("computeLambdaFast uses degree effects properly",{
  s <- new(Stat,times,sen-1,rec-1,N,M,P)
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
               "rid"=matrix(0,K,K))
  P <-11
  beta <- abind(beta,rev.along=3)
  # Add in some degree effects
  beta[8,1,1] <- 4
  beta[9,1,1] <- 3
  beta[10,1,1] <- 2
  beta[11,1,1] <- 1
  
  sij <- s$get_s(3,i-1,j-1)
  lam <- compute_lambda_fast(i-1,j-1,zi-1,zj-1,sij,beta,N,K,P)
  
  sij[8:11] <- log(sij[8:11] + 1)
  ans <- as.numeric(beta[,zi,zj] %*% sij)
  expect_that(lam,equals(ans))
  
})

