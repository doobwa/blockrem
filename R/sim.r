library(releventhier)
library(ggplot2)
require(abind)
require(testthat)
source("R/brem.cpp.r")
source("R/brem.r")
source("R/sbm.r")
N <- 10
K <- 2
beta <- list("intercept"=matrix(c(2,1,1,2),K,K),
             "abba" = matrix(c(3,0,0,0),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(c(0,0,0,2),K,K),
             "abab"=matrix(c(0,-1,-1,0),K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K))
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
ggsave("figs/syn/mat.pdf",width=3,height=3)

# # Temporary testing
# lrm <- brem.lrm(sim$A,N,z,beta.init)
# llk <- brem.llk(sim$A,N,z,beta)
# beta.init <- beta + rnorm(P*K^2,0,.1)
# beta.init[7:11,,] <- 0
# brem.llk(sim$A,N,z,beta.init)
# A <- sim$A
# brem.lpost(A,N,K,z,beta)
# times <- A[,1]
# sen <- A[,2]-1
# rec <- A[,3]-1
# lrm <- brem$lrm(beta,times,sen,rec,z-1,N,M,K,P)
# brem$llk2(lrm,times,sen,rec,N,M)
# brem.llk(sim$A,N,z,beta)
# beta.init <- beta + rnorm(P*K^2,0,.01)
# brem.llk(sim$A,N,z,beta.init)

true.lpost <- brem.lpost(sim$A,N,K,z,beta)
true.lpost

test_that("simulated lrm agrees with brem.lrm",{
  tmp <- brem.lrm(sim$A,N,z,beta)
  for (i in 1:M) diag(tmp[i,,]) <- -Inf
  diag(sim$lrm[1,,]) <- -Inf
  expect_that(all.equal(tmp,sim$lrm),is_true())
})

set.seed(4)
niter <- 300
px <- rep(1,11)
px[7:11] <- 0
fit0 <- sbm.mcmc(sim$A,N,K,niter=niter,z=z,gibbs=FALSE)
fit1 <- brem.mcmc(sim$A,N,K,model.type="diag.rem",niter=niter,z=z,gibbs=FALSE)
fit2 <- brem.mcmc(sim$A,N,K,model.type="full",niter=niter,z=z,gibbs=FALSE)
fit3 <- brem.mcmc(sim$A,N,1,model.type="full",niter=niter,gibbs=FALSE,mcmc.sd=.05)
save(sim,true.lpost,fit0,fit1,fit2,fit3,file="data/syn/fits.rdata")

llks <- melt(list(base=fit0$llks,diag=fit1$llks,full=fit2$llks,sing=fit3$llks))
llks$iter <- 1:niter
qplot(iter,value,data=subset(llks,iter>50),geom="line",colour=factor(L1)) + geom_abline(intercept=true.lpost,slope=0) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
ggsave("figs/syn/logposterior.pdf",height=4,width=5)

# Compute dyad counts for each pshift
source("R/utils.r")
df <- dyad.ps(sim$A,N)
df <- melt(df)
df$i <- z[df$X1]
df$j <- z[df$X2]
qplot(X3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="shift type",y="count for a given dyad") + opts(axis.text.x=theme_text(angle=90))
ggsave("figs/syn/counts.pdf",width=6,height=4)

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

# Prediction experiment on test data: precision
M <- 5000
test <- simulate.brem(M,N,z,beta)
table(test$A[,2],test$A[,3])
lrms <- list(unif = array(1,c(M,N,N)),
             true = brem.lrm(test$A,N,z,beta),
             base = sbm.lrm(test$A,N,fit0$z,fit0$beta),
             diag = brem.lrm(test$A,N,fit1$z,fit1$beta),
             full = brem.lrm(test$A,N,fit2$z,fit2$beta),
             sing = brem.lrm(test$A,N,fit3$z,fit3$beta))
ps <- lapply(lrms,function(lrm) {
  recall(ranks(test$A,-lrm,ties.method="random"),top=1:100)
})
res <- melt(ps,id.vars=c("k"),measure.vars="recall")
qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
qplot(k,value,data=subset(res,k<50),geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
ggsave("figs/syn/test-recall.pdf",width=5,height=4)

# Compute out of sample log posterior
lposts <- list(true = brem.lpost(test$A,N,K,z,beta),
               base = sbm.lpost(test$A,N,K,fit0$z,fit0$beta),
               diag = brem.lpost(test$A,N,K,fit1$z,fit1$beta),
               full = brem.lpost(test$A,N,K,fit2$z,fit2$beta),
               sing = brem.lpost(test$A,N,K,fit3$z,fit3$beta))
unlist(lposts)
save(lposts,file="data/syn/lpost.rdata")

# Look at bias of estimates
true <- brem.lrm(sim$A,N,z,beta)
full <- brem.lrm(sim$A,N,fit2$z,fit2$beta)
dimnames(beta) <- NULL
b <- melt(list(beta,fit2$beta))
qplot(X1,value,data=b,geom="point",colour=factor(L1)) + facet_grid(X2~X3) + theme_bw()
ggsave("figs/syn/bias.pdf",width=5,height=5)