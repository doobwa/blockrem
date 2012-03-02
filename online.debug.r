
get.pred.baseline <- function(train,A,test.ix,model="online") {
  
  get.ms <- switch(model,
                   "uniform" = function(x) array(1,c(nrow(x),N,N)),
                   "online"  = function(x) ratemat.online(x,N),
                   "marg"    = function(x) {
                     b <- table(factor(train[, 2], 1:N), factor(train[, 3], 1:N))
                     rowrates <- rowSums(b)
                     colrates <- colSums(b)
                     r <- rowrates %*% t(colrates)
                     ma <- array(0, c(nrow(x), N, N))
                     for (i in 1:nrow(x)) ma[i, , ] <- r
                     return(ma)
                   })
  
  eps <- 1  # smoothing
  cat("lambdas (train)\n")
  m.train <- get.ms(train)
  m.train[which(m.train==-Inf)] <- 0
  for (i in 1:nrow(train)) {
    lam <- m.train[i,,] + eps
    m.train[i,,] <- lam/sum(lam)
  }
  m.test <- get.ms(A)
  m.test <- m.test[test.ix,,]
  m.test[which(m.test==-Inf)] <- 0
  for (i in 1:length(test.ix)) {
    lam <- m.test[i,,] + eps
    m.test[i,,] <- lam/sum(lam)
  }
  # Get lambda estimates using global rate
  lrm.train <- log(m.train * nrow(train)/train[nrow(train),1])
  lrm.test  <- log(m.test  * nrow(train)/train[nrow(train),1])
  return(list(lrm=list(train=lrm.train,test=lrm.test),m=list(train=m.train,test=m.test)))
}
get.pred <- function(train,A,test.ix,fit) {
  cat("precomputing\n")
  P <- 13
  strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
  strain$precompute()
  stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
  stest$precompute()
  lrm <- list(train = brem.lrm.fast(nrow(train), strain, fit$z, fit$beta),test=NULL)
#              test  = brem.lrm.fast2(stest, res$z, res$beta,test.ix))
  # Compute multinomial likelihoods
  m.train <- exp(lrm$train)
  m.test <- NULL
#   m.test  <- exp(lrm$test)
  for (i in 1:nrow(train)) {
    m.train[i,,] <- m.train[i,,]/sum(m.train[i,,],na.rm=TRUE)
  }
#   for (i in 1:length(test.ix)) {
#     m.test[i,,] <- m.test[i,,]/sum(m.test[i,,],na.rm=TRUE)
#   }
  list(lrm=lrm,m=list(train=m.train,test=m.test))
}

# Compare fit of abba, abby, dyad count to baseline on eckmann
opts <- list(dataset="eckmann-small",numclusters=K,model.type="full",gibbs="none",numiterations=10,slice=TRUE)
library(brem)
load("data/eckmann-small.rdata")
test.ix <- 2001:nrow(A)
M <- nrow(train)
P <- 13
K <- 1#opts$numclusters
s <- new(RemStat,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M,P)
s$precompute()
px <- rep(0,13)
px[c(2,3,12)] <- 1
fit <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                 outfile=NULL)
save(fit,file="tmp.rdata")

preds <- get.pred(train,A,test.ix,fit)
mllk.train <- log(multinomial.score(preds$m$train,train))
llk.train <- loglikelihood_vec_from_lrm(preds$lrm$train,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M)

preds.online <- get.pred.baseline(train,A,test.ix,"online")
mllk.train.online <- log(multinomial.score(preds.online$m$train,train[1:M,]))
llk.train.online <- loglikelihood_vec_from_lrm(preds.online$lrm$train,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M)

pdf("figs/onlinedebug.eckmann.pdf",width=8,height=4)
par(mfrow=c(1,2))
plot(mllk.train,mllk.train.online,asp=1,pch=".")
abline(0,1)
plot(llk.train,llk.train.online,asp=1,pch=".")
abline(0,1)
dev.off()
mean(mllk.train.online)
mean(mllk.train)
mean(llk.train.online)
mean(llk.train)

load("data/synthetic.rdata")
test.ix <- 2001:nrow(A)
M <- nrow(train)
fit <- list(beta=beta,z=z)

preds <- get.pred(train,A,test.ix,fit)
mllk.train <- log(multinomial.score(preds$m$train,train))
llk.train <- loglikelihood_vec_from_lrm(preds$lrm$train,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M)

preds.online <- get.pred.baseline(train,A,test.ix,"online")
mllk.train.online <- log(multinomial.score(preds.online$m$train,train[1:M,]))
llk.train.online <- loglikelihood_vec_from_lrm(preds.online$lrm$train,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M)

pdf("figs/onlinedebug.synthetic.pdf",width=8,height=4)
par(mfrow=c(1,2))
plot(mllk.train,mllk.train.online)
abline(0,1)
plot(llk.train,llk.train.online)
abline(0,1)
dev.off()
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

