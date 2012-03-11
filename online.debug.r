
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
  lrm <- list()
  lrm$train <- brem.lrm.fast(nrow(train), strain, fit$z, fit$beta)
  lrm$test  <- brem.lrm.fast(nrow(A), stest, fit$z, fit$beta)
  lrm$test  <- lrm$test[test.ix,,]
  
  
  # Compute multinomial likelihoods
  m.train <- exp(lrm$train)
  m.test <- exp(lrm$test)
  for (i in 1:nrow(train)) {
    m.train[i,,] <- m.train[i,,]/sum(m.train[i,,])
  }
  for (i in 1:length(test.ix)) {
    m.test[i,,] <- m.test[i,,]/sum(m.test[i,,],na.rm=TRUE)
  }
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
K <- 1
opts <- list(dataset="eckmann-small",numclusters=K,model.type="full",gibbs="none",numiterations=10,slice=TRUE)
library(brem)
load("data/eckmann-small.rdata")
source("pkg/R/brem.r")
test.ix <- 2001:nrow(A)
M <- nrow(train)
P <- 13
K <- 1#opts$numclusters
s <- new(RemStat,train[1:M,1],as.integer(train[1:M,2])-1,as.integer(train[1:M,3])-1,N,M,P)
s$precompute()

# Test that llks agree using lrm (pc)
beta <- array(0,c(P,K,K))
beta[12] <- 1
z <- rep(1,N)

lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)

times <- train[,1];sen <- train[,2];rec <- train[,3]

llk.a <- RemLogLikelihood(beta,times,sen-1,rec-1,z-1,N,M,K,P)
llk.b <- RRemLogLikelihoodFromArraySlow(lrm,times,sen-1,rec-1,N,M)
llk.c <-  RemLogLikelihoodFromArray(lrm,times,sen-1,rec-1,N,M)
llk.d <-  RemLogLikelihoodVecFromArray(lrm,times,sen-1,rec-1,N,M)
llk.e <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
beta[12] <- .15
llk.e2 <- RemLogLikelihoodPc(beta,z-1,s$ptr(),K)
sum(llk.e2)

px <- rep(0,13)
px[c(12)] <- 1
fit.1 <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                 outfile=NULL)

lrm.1 <- LogIntensityArrayPc(fit.1$beta,z-1,s$ptr(),K)

px[c(2,3,12)] <- 1
fit.2 <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL)
fit.0 <- fit.1
fit.0$beta[12] <- 1
fit.0$beta[-12] <- 0
px[c(1,2,3,12)] <- 1
fit.3 <- brem.mcmc(train[1:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL,fix.first=FALSE)
fit.3=fit.2
save(fit.0,fit.1,fit.2,fit.3,file="tmp.rdata")


load("tmp.rdata")

preds <- list(get.pred.baseline(train,A,test.ix,"online"),
              get.pred(train,A,test.ix,fit.0),
              get.pred(train,A,test.ix,fit.1),
              get.pred(train,A,test.ix,fit.2),
              get.pred(train,A,test.ix,fit.3))

# Compute brem llk for multinomial from model multiplied by baserate
preds[[5]] <- preds[[4]]
preds[[5]]$lrm <- list(train=log(preds[[4]]$m$train * nrow(train)/train[nrow(train),1]),
                       test =log(preds[[4]]$m$test * nrow(train)/train[nrow(train),1]))

mllk.trains <- lapply(preds,function(p) log(multinomial.score(p$m$train,train)))
llk.trains <- lapply(preds,function(p) RemLogLikelihoodVecFromArray(p$lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train)))

mllk.tests <- lapply(preds,function(p) log(multinomial.score(p$m$test,test)))
llk.tests <- lapply(preds,function(p) RemLogLikelihoodVecFromArray(p$lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test)))

names(mllk.trains) <- names(llk.trains) <- names(mllk.tests) <- names(llk.tests) <- c("baseline","beta.12 fixed","beta.12","beta.1,2,12","beta.0,1,2,12")#,"beta.0,1,2,12 alt")


# Compute likelihoods
rbind(mllk.train=sapply(mllk.trains,sum),
      llk.train =sapply(llk.trains,sum),
      mllk.test =sapply(mllk.tests,sum),
      llk.test  =sapply(llk.tests,sum))


pdf("figs/onlinedebug.eckmann.pdf",width=12,height=10)
par(mfrow=c(2,4))
for (i in 1:4) { plot(mllk.trains[[5]],mllk.trains[[i]],asp=1,pch=".",xlab="baseline",ylab="model",main="mult. llk")
abline(0,1)
}
for (i in 1:4) {
  plot(llk.trains[[5]], llk.trains[[i]],asp=1,pch=".",xlab="baseline",ylab="model",main="brem llk")
  abline(0,1)
}
dev.off()


# TODO: Idea: Devide lambdas by lambda_12, and get llk of lam.hat * lam
# rm2 <- exp(preds[[3]]$lrm$train)
# for (m in 1:nrow(train)) rm2[i,,] <- rm2[i,,] / rm2[i,1,2]
# lam.hat <- nrow(train)/train[nrow(train),1]
# lrm2 <- log( lam.hat * rm2/sum(rm2))
# sum(loglikelihood_vec_from_lrm(lrm2,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M))

# Plot loglikelihood of events over time
par(mfrow=c(1,1))
plot(llk.trains[[3]],pch=".")
points(llk.trains[[4]],pch=".",col="red")

# Examine misfit of temporal component
total.lambda <- sapply(1:nrow(train),function(i) sum(exp(preds[[5]]$lrm$train[i,,])))
par(mfrow=c(1,1))
time.hat <- cumsum(1/total.lambda)
time.fixed <- cumsum(rep(train[nrow(train),1]/nrow(train),nrow(train)))
time.obs <- train[,1]
plot(time.obs,type="l")
lines(time.hat)
lines(time.fixed)


fit.3 <- brem.mcmc(train[100:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL)
preds.3 <- get.pred(train,A,test.ix,fit.3)
total.lambda.3 <- sapply(1:nrow(train),function(i) sum(exp(preds.3$lrm$train[i,,])))
time.hat <- cumsum(1/total.lambda.3)
lines(time.hat,col="red")

# Make plot of rate of popular dyad over time
i <- 75; j <- 85
plot(exp(preds[[4]]$lrm$train[1:500,i,j]),type="S",ylab="rate",xlab="m")
v <- sapply(1:nrow(train), function(m) i %in% train[m,] | j %in% train[m,])
points(which(v),rep(.2,sum(v)),pch=17)
lines(exp(preds[[5]]$lrm$train[1:500,i,j]),type="S",col="red")


