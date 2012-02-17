
library(ggplot2)
library(brem)

N <- 10
K <- 2
beta <- list("intercept"=matrix(c(0,-1,-1,0),K,K),
             "abba" = matrix(c(3,0,0,0),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(c(0,0,0,2),K,K),
             "abab"=matrix(c(0,-1,-1,0),K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(c(2,0,0,0),K,K),
             "rid"=matrix(0,K,K),
#             "dc" =matrix(c(1,0,0,2),K,K),
             "dc" =matrix(0,K,K),
             "cc" =matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)

M <- 2000
set.seed(1)
z <- c(rep(1,5),rep(2,5))

sim <- simulate.brem(M,N,z,beta)
mat <- table(sim$A[,2],sim$A[,3])
mat <- melt(as.matrix(mat))
colnames(mat) <- c("X1","X2","value")
plotmat(mat)
A <- train <- sim$A
plot(A[,1],type="l")

test <- simulate.brem(5000,N,z,beta)$A

# Compute log likelihood and log posterior of data under true model
times <- sim$A[,1]
sen <- sim$A[,2]
rec <- sim$A[,3]
s <- new(brem:::RemStat,times,sen-1,rec-1,N,M,P)
s$precompute()
llk.true <- loglikelihood_fast(beta,z-1,s$ptr(),K)  
true.lpost <- brem.lpost.fast(A,N,K,z,s,beta)

save(A,sim,N,K,P,M,z,beta,train,test,true.lpost,file="data/synthetic.rdata")
#ggsave("figs/syn/mat.pdf",width=3,height=3)


test_that("simulated lrm agrees with brem.lrm",{
  tmp <- brem.lrm(sim$A,N,z,beta)
  for (i in 1:M) diag(tmp[i,,]) <- -Inf
  diag(sim$lrm[1,,]) <- -Inf
  expect_that(all.equal(tmp[-1,,],sim$lrm[-1,,]),is_true())
})


## OLD CODE
# Look at distance between true and estimated parameter vectors
fit <- fit2
dist <- function(x,y) sqrt(sum((x-y)^2))
ds <- sapply(1:niter,function(i) {
  c(dist(beta[2:6,1,1],fit$param[i,1,1,2:6]),# - beta[1,1,1]),
    dist(beta[2:6,2,2],fit$param[i,2,2,2:6]),#, - beta[1,2,2]),
    dist(beta[2:6,2,1],fit$param[i,2,1,2:6]),#, - beta[1,2,1]),
    dist(beta[2:6,1,2],fit$param[i,1,2,2:6]))# - beta[1,1,2]))
})
ds <- melt(t(ds))
qplot(X1,value,data=ds,geom="line",colour=factor(X2)) + theme_bw() + labs(x="iteration",colour="block",y="Euclidean distance to truth")
ggsave("figs/syn/bias.pdf",width=5,height=4)
