require(Rcpp)
require(inline)
require(testthat)
#gctorture(TRUE)

inc <- paste(readLines("pkg/inst/tests/tmp.cpp"))
fx <- cxxfunction(,"",includes=inc, plugin="Rcpp")

brem <- Module("brem",getDynLib(fx))
load("data/synthetic.rdata")
K <- 5
load("results/synthetic/full.10.rdata")

for (i in 1:500) {
print(i)
    b <- brem$gibbs_with_init(A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,res$beta,res$z-1,N,M,P,K)
}

set.seed(2)
s <- new(brem$RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,M,P)
s$precompute()
a <- s$ptr()
b <- brem$loglikelihood_fast(res$beta,res$z-1,a,K)
b <- brem$gibbs(res$beta,res$z-1,a,K)




#
