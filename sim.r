# 
# beta <- c(1,0)
# effects <- names(beta) <- c("PSAB-BA","PSAB-BY")
# A <- simulate.rem(1000,10,effects,beta)
# 
library(ggplot2)
source("rem.cpp.r")
source("fns.r")
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
plotmat(mat)
ggsave("figs/syn-mat.pdf",width=3,height=3)

brem.llk(sim$A,N,z,beta)
brem.lpost(sim$A,N,z,beta)

set.seed(4)
niter <- 100
beta.init <- beta + rnorm(length(beta),0,1)
fit <- brem.mcmc(sim$A,N,K,P,niter=niter,init=beta.init)

# Compare llk and lpost of true and fit
brem.llk(sim$A,N,z,beta)
brem.llk(sim$A,N,fit$z,fit$beta)
brem.lpost(sim$A,N,z,beta)
brem.lpost(sim$A,N,fit$z,fit$beta)

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
