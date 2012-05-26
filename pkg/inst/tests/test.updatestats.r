

test_that("statistics creation works",{
  
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
            "rid" = matrix(0,N,N),
            "dc"=matrix(0,N,N),
            "cc"=matrix(0,N,N))
  P <- length(s)
  s <- abind(s,rev.along=3)
  i <- 1
  j <- 2
  a <- 2
  b <- 1
  m <- 1
  s <- UpdateStatisticsArray(s,1,a-1,b-1,N,P)
  i <- 1
  j <- 2
  dimnames(s) <- list(NULL,NULL,NULL)
  expect_that(s[2,i,j],equals(1))  # abba
  expect_that(s[3,i,j],equals(0)) # abby
  expect_that(s[10,i,j],equals(1)) # sid
  expect_that(s[8,i,j],equals(0)) # sod
  expect_that(s[9,i,j],equals(1)) # rod
})