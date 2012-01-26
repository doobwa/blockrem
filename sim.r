library(releventhier)
library(ggplot2)
source("rem.cpp.r")
source("fns.r")
source("sbm.r")
N <- 10
K <- 2
P <- 7
beta <- array(0,c(K,K,P))
beta[1,1,] <- c(2,3,0,0,0,0,0)
beta[2,1,] <- c(1,0,0,0,0,0,1)
beta[1,2,] <- c(1,0,0,0,0,0,1) 
beta[2,2,] <- c(2,0,0,0,0,2,0)  # AB-AY


M <- 2000
set.seed(1)
z <- c(rep(1,5),rep(2,5))
sim <- simulate.brem(M,N,z,beta)
mat <- table(sim$A[,2],sim$A[,3])
mat <- melt(as.matrix(mat))
colnames(mat) <- c("X1","X2","value")
require(releventhier)
plotmat(mat)
#ggsave("figs/syn-mat.pdf",width=3,height=3)


px <- rep(1,7)
brem.llk(sim$A,N,z,beta,px)
brem.lpost(sim$A,N,z,beta,px)

set.seed(4)
niter <- 100
beta.init <- beta + rnorm(length(beta),0,1)
px0 <- px1 <- rep(1,7)
px2 <- c(0,1,1,1,1,1,1)
fit0 <- sbm.mcmc(sim$A,N,K)
fit1 <- brem.mcmc(sim$A,N,K,P,px,model.type="diag.rem",niter=niter,init=beta.init)
fit2 <- brem.mcmc(sim$A,N,K,P,px,model.type="full",niter=niter,z=z,beta=beta)

# Compare llk and lpost of true and fit
fit <- fit2
brem.llk(sim$A,N,z,beta,px)
brem.llk(sim$A,N,fit$z,fit$beta,px)
brem.lpost(sim$A,N,z,beta,px)
brem.lpost(sim$A,N,fit$z,fit$beta,px)

pdf("figs/syn-llk.pdf",width=4,height=4)
plot(fit$llk[1:100],type="l",ylab="loglikelihood",xlab="iteration")
dev.off()

# Compute dyad counts for each pshift
df <- dyad.ps(sim$A,N)
df <- melt(df)
df$i <- z[df$X1]
df$j <- z[df$X2]
qplot(X3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="shift type",y="count for a given dyad") + opts(axis.text.x=theme_text(angle=90))
ggsave("figs/syn-counts.pdf",width=6,height=4)

# Look at distance between true and estimated parameter vectors
dist <- function(x,y) sqrt(sum((x-y)^2))
ds <- sapply(1:niter,function(i) {
  c(dist(beta[1,1,],fit$param[i,1,1,]),
    dist(beta[2,2,],fit$param[i,2,2,]))
})
ds <- sapply(1:niter,function(i) {
  c(dist(beta[1,1,2:6],fit$param[i,1,1,2:6] - beta[1,1,1]),
    dist(beta[2,2,2:6],fit$param[i,2,2,2:6] - beta[2,2,1]),
    dist(beta[2,1,2:6],fit$param[i,2,1,2:6] - beta[2,1,1]),
    dist(beta[1,2,2:6],fit$param[i,1,2,2:6] - beta[1,2,1]))
})
ds <- melt(t(ds))
qplot(X1,value,data=ds,geom="line",colour=factor(X2)) + theme_bw()


# Prediction experiment: precision
library(releventhier)
test <- simulate.brem(M,N,z,beta)
lrms <- list(unif = array(1,c(M,N,N)),
             true = brem.lrm(test$A,N,z,beta,px0),
             base = sbm.lrm(test$A,N,fit0$z,fit0$beta),
             diag = brem.lrm(test$A,N,fit1$z,fit1$beta,px1),
             full = brem.lrm(test$A,N,fit2$z,fit2$beta,px2))
ps <- lapply(lrms,function(lrm) {
  precision(ranks(test$A,-lrm,ties.method="random"))
})
res <- melt(ps,id.vars=c("k"),measure.vars="precision")
qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw()

# Compute out of sample log posterior
lposts <- list(true = brem.lpost(test$A,N,z,beta,px0),
               base = sbm.lpost(test$A,N,K,fit0$z,fit0$beta),
               diag = brem.lpost(test$A,N,fit1$z,fit1$beta,px1),
               full = brem.lpost(test$A,N,fit2$z,fit2$beta,px2))
lposts
