#gctorture(TRUE)
# 
# inc <- paste(readLines("pkg/inst/tests/tmp.cpp"))
# fx <- cxxfunction(,"",includes=inc, plugin="Rcpp")
# 
# brem <- Module("brem",getDynLib(fx))
# 
# load("data/synthetic.rdata")
# K <- 5
# load("results/synthetic/full.10.rdata")
# set.seed(2)
# s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,M,P)
# s$precompute()
# a <- s$ptr()
# b <- loglikelihood_fast(res$beta,res$z-1,a,K)
# b <- gibbs(res$beta,res$z-1,a,K)


# for (i in 1:500) {
# print(i)
#     b <- brem$gibbs_with_init(A[,1],A[,2]-1,A[,3]-1,res$beta,res$z-1,N,M,P,K)
# }
# 
# for (i in 1:500) {
# print(i)
#     b <- brem$gibbs_with_init(A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,res$beta,res$z-1,N,M,P,K)
# }
# 
# load("data/eckmann-small.rdata")
# K <- 5
# z <- sample(1:K,N,replace=TRUE)
# beta <- array(rnorm(K*K*P),c(P,K,K))
# b <- brem$gibbs_with_init(train[,1],train[,2]-1,train[,3]-1,beta,z-1,N,M,P,K)
# s <- new(brem$RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
# s$precompute()
# b <- brem$gibbs(beta,z-1,s$ptr(),K)
# 
# 
# #
