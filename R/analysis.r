
llks <- melt(list(base=fit0$llks,diag=fit1$llks,full=fit2$llks,sing=fit3$llks))
llks$iter <- 1:niter
qplot(iter,value,data=subset(llks,iter>50),geom="line",colour=factor(L1)) + geom_abline(intercept=true.lpost,slope=0) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
ggsave("figs/syn/logposterior2.pdf",height=4,width=5)

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
