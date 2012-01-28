library(releventhier)
library(ggplot2)
source("R/rem.cpp.r")
source("R/brem.r")
source("R/sbm.r")
N <- 10
K <- 2
P <- 7
beta <- array(0,c(K,K,P))
beta[1,1,] <- c(2,3,0,0,0,0,0)
beta[2,1,] <- c(1,0,0,0,0,0,-1)
beta[1,2,] <- c(1,0,0,0,0,0,-1) 
beta[2,2,] <- c(2,0,0,0,0,2,0)  # AB-AY


M <- 2000
set.seed(1)
z <- c(rep(1,5),rep(2,5))
px <- c(1,1,1,1,1,1,1)
sim <- simulate.brem(M,N,z,beta,px)
mat <- table(sim$A[,2],sim$A[,3])
mat <- melt(as.matrix(mat))
colnames(mat) <- c("X1","X2","value")
plotmat(mat)
ggsave("figs/syn/mat.pdf",width=3,height=3)

brem.llk(sim$A,N,K,z,beta,px)
brem.lpost(sim$A,N,K,z,beta,px)
llk.true <- brem.lpost(sim$A,N,K,z,beta,px)

# Make sure simulated lrm agrees with brem.lrm
tmp <- brem.lrm(sim$A,N,K,z,beta,px)
tmp[which(tmp == -15)] <- -Inf
all.equal(tmp,sim$lrm)

set.seed(4)
niter <- 500
#fit0 <- sbm.mcmc(sim$A,N,1,niter=niter)
fit0 <- sbm.mcmc(sim$A,N,K,niter=niter,z=z)
fit1 <- brem.mcmc(sim$A,N,K,px,model.type="diag.rem",niter=niter,z=z,gibbs=FALSE)
fit2 <- brem.mcmc(sim$A,N,K,px,model.type="full",niter=niter,z=z,gibbs=FALSE)
fit3 <- brem.mcmc(sim$A,N,1,px,model.type="full",niter=niter,gibbs=FALSE)
save(fit0,fit1,fit2,fit3,file="data/syn-fits.rdata")

llks <- melt(list(base=fit0$llks,diag=fit1$llks,full=fit2$llks,sing=fit3$llks))
llks$iter <- 1:niter
qplot(iter,value,data=llks,geom="line",colour=factor(L1)) + geom_abline(intercept=llk.true)

# Compare llk and lpost of true and fit
fit <- fit2
brem.llk(sim$A,N,z,beta,px)
brem.llk(sim$A,N,fit$z,fit$beta,px)
brem.lpost(sim$A,N,z,beta,px)
brem.lpost(sim$A,N,fit$z,fit$beta,px)

pdf("figs/syn/llk.pdf",width=4,height=4)
plot(fit$llk[1:300],type="l",ylab="loglikelihood",xlab="iteration")
dev.off()

# Compute dyad counts for each pshift
df <- dyad.ps(sim$A,N)
df <- melt(df)
df$i <- z[df$X1]
df$j <- z[df$X2]
qplot(X3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="shift type",y="count for a given dyad") + opts(axis.text.x=theme_text(angle=90))
ggsave("figs/syn/counts.pdf",width=6,height=4)

# Look at distance between true and estimated parameter vectors
dist <- function(x,y) sqrt(sum((x-y)^2))
ds <- sapply(1:niter,function(i) {
  c(dist(beta[1,1,2:6],fit$param[i,1,1,2:6] - beta[1,1,1]),
    dist(beta[2,2,2:6],fit$param[i,2,2,2:6] - beta[2,2,1]),
    dist(beta[2,1,2:6],fit$param[i,2,1,2:6] - beta[2,1,1]),
    dist(beta[1,2,2:6],fit$param[i,1,2,2:6] - beta[1,2,1]))
})
ds <- melt(t(ds))
qplot(X1,value,data=ds,geom="line",colour=factor(X2)) + theme_bw() + labs(x="iteration",colour="block",y="Euclidean distance to truth")
ggsave("figs/syn/bias.pdf",width=5,height=4)

# Prediction experiment on test data: precision
source("R/utils.r")
M <- 10000
test <- simulate.brem(M,N,z,beta,px)
lrms <- list(unif = array(1,c(M,N,N)),
             true = brem.lrm(test$A,N,K,z,beta,px),
             base = sbm.lrm(test$A,N,fit0$z,fit0$beta),
             diag = brem.lrm(test$A,N,K,fit1$z,fit1$beta,px),
             full = brem.lrm(test$A,N,K,fit2$z,fit2$beta,px),
             sing = brem.lrm(test$A,N,1,fit3$z,fit3$beta,px))
ps <- lapply(lrms,function(lrm) {
  recall(ranks(test$A,-lrm,ties.method="random"),top=1:100)
})
res <- melt(ps,id.vars=c("k"),measure.vars="recall")
qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
qplot(k,value,data=subset(res,k<20),geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
ggsave("figs/syn/test-recall.pdf",width=5,height=4)

# Compute out of sample log posterior
lposts <- list(true = brem.lpost(test$A,N,K,z,beta,px),
               base = sbm.lpost(test$A,N,K,fit0$z,fit0$beta),
               diag = brem.lpost(test$A,N,K,fit1$z,fit1$beta,px),
               full = brem.lpost(test$A,N,K,fit2$z,fit2$beta,px),
               sing = brem.lpost(test$A,N,K,fit3$z,fit3$beta,px))
lposts



true = brem.lrm(sim$A,N,K,z,beta,px)
full = brem.lrm(sim$A,N,K,fit2$z,fit2$beta,px)

sim2 <- function(lrm) {
  M <- dim(lrm)[1]
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  for (i in 1:(M-1)) {
    lambda <- lrm[i,,]
    diag(lambda) <- -Inf
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    i <- cells[drawcell,1]
    j <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    A <- rbind(A,c(time,i,j))
  }
  return(list(A=A))
}