

load("results/eckmann-small/llks/online.rdata")
llk.test.online <- llk.test
llk.train.online <- llk.train
mllk.test.online <- llkm.test
mllk.train.online <- llkm.train

load("results/eckmann-small/llks/full.3.rdata")
llk.test.full.3 <- llk.test
llk.train.full.3 <- llk.train
mllk.test.full.3 <- llkm.test
mllk.train.full.3 <- llkm.train
load("results/eckmann-small/llks/full.2.rdata")
llk.test.full.2 <- llk.test
llk.train.full.2 <- llk.train
mllk.test.full.2 <- llkm.test
mllk.train.full.2 <- llkm.train
load("results/eckmann-small/llks/countsonly.rdata")

plot(llk.train.online,llk.train.full.3,col="#666666FF")
abline(0,1)
abline(h=0)
abline(v=0)
plot(llk.test.online,llk.test.full.3)
abline(0,1)
abline(h=0)
abline(v=0)
plot(llk.test.online - llk.test.full.2)
plot(mllk.test.online, mllk.test.full.3)
abline(0,1)
mean(llk.test.online - llk.test.full.3)       

df <- data.frame(i=1:nrow(test),online=mllk.test.online,full=mllk.test.full.3,block=paste(res$z[test[,2]],res$z[test[,3]]))
qplot(i,online-full,colour=block,data=df) + scale_color_brewer() + theme_bw()

# Look at how fast our model says the overall rate is
opts$model <- "full.3"
f <- paste("results/",opts$dataset,"/",opts$model,".rdata",sep="")
load(f)
cat("precomputing\n")
P <- 13
strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
strain$precompute()
stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
stest$precompute()
lrm.train <- brem.lrm.fast(nrow(train), strain, res$z, res$beta)
test.ix <- (1:nrow(A))[-(1:nrow(train))]
lrm.test  <- brem.lrm.fast2(stest, res$z, res$beta,test.ix)

lambda.hat.train <- sapply(1:nrow(train),function(i) sum(exp(lrm.train[i,,])))
lambda.hat.test <- sapply(1:length(test.ix),function(i) sum(exp(lrm.test[i,,])))
nrow(train)/train[nrow(train),1]

ix <- 100:2000
deltas <- c(diff(train[ix,1]),0)
plot(lambda.hat.train[ix],1/deltas/length(ix),asp=1)

# Fit submodels to the first portion of data and compare llks to online baseline.  Do we win on multinomial?  If so, is our misspecification of time screwing us up?
library(brem)
load("data/eckmann-small.rdata")
M <- 1500#nrow(train)
P <- 13
K <- 1#opts$numclusters
s <- new(RemStat,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M,P)
s$precompute()
px <- rep(0,13)
px[c(2,3,12)] <- 1
opts=list(dataset="eckmann-small",numclusters=K,model.type="full",gibbs="none",numiterations=10,slice=TRUE)
fit <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                 outfile=NULL)

#lrm <- brem.lrm.fast(M, s, fit$z, fit$beta)
lrm <- brem.lrm(A,N,fit$z,fit$beta)
lamhat <- sapply(1:M,function(i) sum(exp(lrm[i,,])))

m.train <- exp(lrm)
for (i in 1:M) {
  m.train[i,,] <- m.train[i,,]/sum(m.train[i,,],na.rm=TRUE)
}
mllk.train <- log(multinomial.score(m.train,train[1:M,]))
llk.train <- loglikelihood_vec_from_lrm(lrm,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M)
plot(mllk.train,mllk.train.online[1:M])
abline(0,1)
mean(mllk.train.online)
mean(mllk.train)
mean(llk.train.online)
mean(llk.train)
# 
# qqplot(lambda.hat.train,c(1/2000/diff(train[,1]),.5))
# lambda.hat.train <- sapply(1:nrow(train),function(i) sum(exp(lrm.train[i,,])))

pdf("time-incorrect.pdf",width=5,height=5)
plot(cumsum(1/lamhat),type="l")
lines(train[,1],type="l")
dev.off()
