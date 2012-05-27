load(paste("data/synthetic-1.rdata",sep=""))
outfile <- c("results/synthetic-1/full.rdata")
library(brem)
# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
ego <- 1
K <- 5
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P,ego)
s$precompute()

## Set priors
effects <- c("intercept","abba","abby","abay")
priors <- list(alpha=1,sigma.proposal=.1,phi=list(mu=0,sigma=1),mu=list(mu=0,sigma=1),sigma=list(alpha=3,beta=1))
fit <- brem(train,N,K,effects)

save(fit,file=outfile)

train.ix <- 1:nrow(train)
test.ix <- 2001:7000
p <- evaluate(A,N,train.ix,test.ix,fit,ties.method="random")
q <- evaluate.baseline(A,N,train.ix,test.ix,"online",ties.method="random")

par(mfrow=c(2,2))
plot(p$llk$train,q$llk$train,xlab="model",ylab="baseline"); abline(0,1)
plot(p$llk$test,q$llk$test,xlab="model",ylab="baseline"); abline(0,1)
plot(p$mllk$train,q$mllk$train,xlab="model",ylab="baseline"); abline(0,1)
plot(p$mllk$test,q$mllk$test,xlab="model",ylab="baseline"); abline(0,1)

# Profile model fitting
Rprof()
#fit <- brem(train,N,K,effects)
Rprof(NULL)
summaryRprof()

# Mixing
library(reshape2)
library(ggplot2)
tmp <- melt(lapply(fit$samples,function(x) x$phi))
qplot(L1, value, data=tmp, geom="line", colour=factor(Var1)) +facet_grid(Var2~Var3,scales="free")
#tmp <- subset(tmp,Var2==k1 & Var3==k2)
#qplot(L1, value, data=tmp, geom="line") +facet_grid(Var1~.,scales="free")

plot.posterior(fit,priors)
