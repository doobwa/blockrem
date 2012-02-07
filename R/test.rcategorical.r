source("R/brem.cpp.r")

K <- 5

test_that("rcategorical does not create illogical values",{

  ps <- rep(10,K)
  x <- sapply(1:1000,function(i) brem$rcategorical(log(ps)))

  expect_that(all(x >=0 & x < K), is_true())
})

test_that("draws make sense with positive values of log probability",{
  ps <- rep(10,K)
  x <- sapply(1:10000,function(i) brem$rcategorical(log(ps)))
  
  expect_that(all(abs(table(x) - 2000) < 200),is_true())
})
test_that("draws make sense with negative values of log probability",{
  lps <- rep(-10,K)
  x <- sapply(1:10000,function(i) brem$rcategorical(lps))

  expect_that(all(abs(table(x) - 2000) < 200),is_true())
})

test_that("draws make sense with typical values log probability",{
  lps <- c(-500,-400,-100)
  x <- sapply(1:10000,function(i) brem$rcategorical(lps))
  
  expect_that(all(x == 2),is_true())
  
  lps <- c(882,100)
  x <- sapply(1:10000,function(i) brem$rcategorical(lps))
  expect_that(all(x==0), is_true())
  lps <- c(100,882)
  x <- sapply(1:10000,function(i) brem$rcategorical(lps))
  expect_that(all(x==1), is_true())
  
})