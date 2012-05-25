context("computing likelihoods and intensity functions")

set.seed(1)
M <- 4
N <- 5
times <- c(0,2,3,4)
sen <- c(1,3,3,1)
rec <- c(3,1,1,3)
z <- rep(1,N)
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
             "rid"=matrix(0,1,1),
             "dc"=matrix(0,1,1),
             "cc"=matrix(0,1,1))
P <- length(beta)
beta <- abind(beta,rev.along=3)
# 
# 
# test_that("lrm and llk functions work on small example for K=1",{
#   # Set up example
#   
#   # Constract log rate matrix by hand and compare to drem$lrm
#   K <- 1
#   lrm <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)
#   
#   # Compute log likelihood by hand.  
#   llks <- RRemLogLikelihoodFromArrayAlt(lrm,times,sen-1,rec-1,N,M)
#   
#   s <- new(RemStat,times,sen-1,rec-1,N,M,P)
#   s$precompute()
#   llk3 <- loglikelihood_fast(beta,z-1,s$ptr(),K)
#   lrm2 <- (beta,z-1,s,M,N,K,P)
# 
#   # log intensity array functions act as expected and agree with previous versions
#   lrm3 <- log_intensity_array_fast(beta,z-1,s$ptr(),K)
#   lrm4 <- log_intensity_array_fast_subset(beta,z-1,s$ptr(),K,0:(M-1))
#   expect_that(lrm3, equals(lrm4))
#   expect_that(all(lrm4==lrm),is_true())
#   
#   true.fast <- c(1,
#                  2 - 2*(19*exp(1) + exp(2)), 
#                  1 - (19*exp(1) + exp(2)), 
#                  2 -  (times[4] - times[1])*6*exp(1))
# 
#   llk4 <- RemLogLikelihoodVecFromArray(lrm,times,sen-1,rec-1,N,M)
# 
#   expect_that(sum(llks),equals(llk2))
#  # expect_that(sum(true.fast),equals(sum(llk4)))  # TODO: check this example again
#   expect_that(sum(llks)-2,equals(sum(llk3)))  # TODO: BUG where old likelihood is off by intercept
#   expect_that(sum(llks)-2,equals(sum(llk4)))
# })
