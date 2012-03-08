
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
multinomial.score <- function(m,x) {
  M <- nrow(x)
  r <- rep(0, M)
  for (i in 1:M) {
    r[i] <- m[i,x[i,2],x[i,3]]
  }
  return(r)
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
#px[c(2,3,12)] <- 1
px[c(12)] <- 1
fit <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                 outfile=NULL)
save(fit,file="tmp.rdata")

load("tmp.rdata")
fit$beta[2:3] <- 0
fit$beta[12] <- 1
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

################
# Make plot of rate of popular dyad over time
i <- 75; j <- 85
i <- 21; j <- 20
plot(preds$m$train[1:300,i,j],type="S")
lines(preds.online$m$train[1:300,i,j],type="S",col="red")

############################################
# After meeting with Smyth

# Compute p_ij(m)
m.train <- ratemat.online(train,N)
m.train[which(m.train==-Inf)] <- 0
eps <- 1
for (i in 1:nrow(train)) {
  lam <- m.train[i,,] + eps
  m.train[i,,] <- lam/sum(lam)
}
m.train.online <- m.train

# Compute log lambda_ij(m)
fit$beta[12] <- 1
P <- 13
lrm <- brem.lrm.fast(nrow(train), strain, fit$z, fit$beta)
m.train <- lrm
for (i in 1:nrow(train)) {
  m.train[i,,] <- exp(m.train[i,,])/sum(exp(m.train[i,,]))
}

# Compare rates for given dyads
i <- 75#sample(1:N,1)
j <- 85#sample(1:N,1)
summary(m.train[,i,j])
plot(m.train[,i,j],type="S")
lines(m.train.online[,i,j],type="S",col="red")

i <- 3#sample(1:N,1)
j <- 5#sample(1:N,1)
summary(m.train[,i,j])
summary(m.train.online[,i,j])
plot(m.train[,i,j],type="S")
lines(m.train.online[,i,j],type="S",col="red")

# Compare prob for observed dyads
obs.lambda <- sapply(1:nrow(train),function(i) m.train[i,train[i,2],train[i,3]])
obs.lambda.online <- sapply(1:nrow(train),function(i) m.train.online[i,train[i,2],train[i,3]])
plot(obs.lambda,obs.lambda.online,pch=".",asp=1)
abline(0,1)

# Now fit model

### OLD
# check lrm is correct
x <- s$get_s(200,72,0)
b <- c(fit$beta)
compute_lambda_fast(72,0,0,0,x,b,N,K,P) == log((9+1)/(199+N*(N-1))) * b[12]
# compare "pij" to "sij"
obs.count.online <- sapply(1:nrow(train),function(i) m.train.online[i,train[i,2],train[i,3]])
obs.count <- sapply(1:nrow(train),function(i) s$get_s(i-1,train[i,2]-1,train[i,3]-1)[12])
obs.count.alt <- sapply(1:nrow(train),function(i) s$get_s(i-1,train[i,2]-1,train[i,3]-1)[12])

m.train <- ratemat.online(train,N)
m.train[200,73,1]
lrm.train[200,73,1]
log(m.train.online[200,73,1])

obs.lambda <- sapply(1:nrow(train),function(i) lrm.train[i,train[i,2],train[i,3]])
obs.lambda.online <- sapply(1:nrow(train),function(i) log(m.train.online[i,train[i,2],train[i,3]]))
plot(obs.lambda,obs.lambda.online,pch=".");abline(0,1)

# Compute lambda_ij / sum_ij lambda_ij
m.train <- exp(lrm.train)/exp(1)
for (i in 1:nrow(train)) {
  lam <- m.train[i,,]
  m.train[i,,] <- lam/sum(lam)
}

obs.lambda <- sapply(1:nrow(train),function(i) m.train[i,train[i,2],train[i,3]])
obs.lambda.online <- sapply(1:nrow(train),function(i) m.train.online[i,train[i,2],train[i,3]])
plot(obs.lambda/exp(1),obs.lambda.online,pch=".");abline(0,1)




########################################
preds <- get.pred(train,A,test.ix,fit)
obs.count.online <- sapply(1:nrow(train),function(i) tmp[i,train[i,2],train[i,3]])
obs.count <- sapply(1:nrow(train),function(i) s$get_s(i-1,train[i,2]-1,train[i,3]-1)[12])
plot(obs.count,obs.count.online)
obs.lambda.online <- sapply(1:nrow(train),function(i) m.train[i,train[i,2],train[i,3]])
obs.lambda.online2 <- sapply(1:nrow(train),function(i) preds.online$m$train[i,train[i,2],train[i,3]])
obs.lambda <- sapply(1:nrow(train),function(i) lrm$train[i,train[i,2],train[i,3]])# preds$m$train[i,train[i,2],train[i,3]])
obs.lambda2 <- sapply(1:nrow(train),function(i) lrm.train[i,train[i,2],train[i,3]])
plot(obs.lambda.online,obs.lambda,pch=".")
plot(obs.lambda2,obs.lambda,pch=".")
plot(obs.lambda.online2,obs.lambda2,pch=".");abline(0,1)

cbind(train[1:30,],obs.lambda2[1:30],obs.lambda.online2[1:30])

P <- 13
strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
strain$precompute()
lrm <- list(train = brem.lrm.fast(nrow(train), strain, fit$z, fit$beta),test=NULL)
lrm.train <- brem.lrm(train,N,fit$z,fit$beta)

m.train <- ratemat.online(train,N)
m.train[which(m.train==-Inf)] <- 0
for (i in 1:nrow(train)) {
  lam <- m.train[i,,] + eps
  m.train[i,,] <- lam/sum(lam)
}



# high rate @ m=1835: 68,3
chosen <- 1835
x <- s$get_s(1835,67,2)
compute_lambda_fast(67,2,0,0,x,b,N,K,P)
x <- s$get_s(1835,74,84)
compute_lambda_fast(74,84,0,0,x,b,N,K,P)







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

