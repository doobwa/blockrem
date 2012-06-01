
# Compare fit of abba, abby, dyad count to baseline on eckmann
opts <- list(dataset="eckmann-small",numclusters=K,model.type="full",gibbs=TRUE,numiterations=10,slice=TRUE)
library(brem)
load("data/eckmann-small.rdata")
A[,1] <- 10 * A[,1]
train[,1] <- 10 * train[,1]
source("pkg/R/brem.r")
test.ix <- 2001:nrow(A)
M <- nrow(train)
P <- 13
K <- 1#opts$numclusters

# Test that llks agree using lrm (pc)
beta <- array(0,c(P,K,K))
beta[12] <- 1
z <- rep(1,N)
times <- train[,1];sen <- train[,2];rec <- train[,3]

lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)

lrm2 <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)
all.equal(lrm,lrm2)


# Check if various ways of computing llk are equal
beta2 <- beta; beta2[12] <- .15  # see if we improve by setting beta manually
llks <- list(RemLogLikelihood(beta,times,sen-1,rec-1,z-1,N,M,K,P),
             RRemLogLikelihoodFromArraySlow(lrm,times,sen-1,rec-1,N,M),
             RemLogLikelihoodFromArray(lrm,times,sen-1,rec-1,N,M),
             RemLogLikelihoodVecFromArray(lrm,times,sen-1,rec-1,N,M),
             RemLogLikelihoodPc(beta,z-1,s$ptr(),K),
             RemLogLikelihoodPc(beta2,z-1,s$ptr(),K))
sapply(llks,sum)

priors <- list(beta=list(mu=0,sigma=3))
px <- rep(0,13)
px[c(12)] <- 1
fit.1 <- brem.mcmc(train[1:M,],N,K,model.type=opts$model.type,mh=FALSE,slice=TRUE,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                 outfile=NULL,priors=priors)
fit.0 <- fit.1
fit.0$beta[12] <- 1
fit.0$beta[-12] <- 0

px <- rep(0,13)
px[c(1,12)] <- 1
fit.2 <- brem.mcmc(train[1:M,],N,K,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL,priors=priors,skip.intercept=FALSE)

px <- rep(0,13)
px[c(1,2,3,12)] <- 1
fit.3 <- brem.mcmc(train[1:M,],N,K,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,px=px,
                   outfile=NULL,skip.intercept=FALSE,priors=priors,beta=fit.2$beta)

px <- rep(0,13)
px[c(1)] <- 1
fit.4 <- brem.mcmc(train[1:M,],N,K,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL,skip.intercept=FALSE,priors=priors)

save(fit.0,fit.1,fit.2,fit.3,fit.4,file="tmp.rdata")


# Examine misfit of temporal component
lrm4 <- LogIntensityArrayPc(fit.4$beta,fit.4$z-1,s$ptr(),K)
for (m in 1:2000) diag(lrm4[m,,]) <- -Inf
lrm <- LogIntensityArrayPc(fit.3$beta,fit.3$z-1,s$ptr(),K)
for (m in 1:2000) diag(lrm[m,,]) <- -Inf
total.lambda4 <- sapply(1:nrow(train),function(i) sum(exp(lrm4[i,,])))
total.lambda  <- sapply(1:nrow(train),function(i) sum(exp(lrm[i,,])))
pdf("figs/online-debug-time.pdf",width=7,height=7)
par(mfrow=c(1,1))
time.hat4 <- cumsum(1/total.lambda4)
time.hat  <- cumsum(1/total.lambda)
time.fixed <- cumsum(rep(train[nrow(train),1]/nrow(train),nrow(train)))
time.obs <- train[,1]
plot(time.obs,type="l")
lines(time.hat,col="red")
lines(time.fixed,col="green")
lines(time.hat4,col="blue")
dev.off()

load("tmp.rdata")

preds <- list(get.pred.baseline(train,A,test.ix,"online"),
              get.pred(train,A,test.ix,fit.0),
              get.pred(train,A,test.ix,fit.1),
              get.pred(train,A,test.ix,fit.2),
              get.pred(train,A,test.ix,fit.3),
              get.pred(train,A,test.ix,fit.4))

# Check that Pc gives same answer
sum(RemLogLikelihoodVecFromArray(preds[[5]]$lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train)))
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
s$precompute()
sum(RemLogLikelihoodPc(fit.3$beta,fit.3$z-1,s$ptr(),K))


# # TODO: Possible bug since intensity for fit0 and baseline don't agree
# # Compare online multinomial with fit.0
# a <- RemLogLikelihoodVecFromArray(preds[[1]]$m$train,times,sen-1,rec-1,N,M)
# b <- RemLogLikelihoodVecFromArray(preds[[1]]$lrm$train,times,sen-1,rec-1,N,M)
# beta[1] <- 8; beta[12] <- 1
# lrm2 <- LogIntensityArray(beta,times,sen-1,rec-1,z-1,N,M,K,P)
# d <- RemLogLikelihoodVecFromArray(lrm2,times,sen-1,rec-1,N,M) 
# sum(a)
# sum(b)
# sum(d)

# 
# # Compute brem llk for multinomial from model multiplied by baserate
# preds[[5]] <- preds[[4]]
# preds[[5]]$lrm <- list(train=log(preds[[4]]$m$train * nrow(train)/train[nrow(train),1]),
#                        test =log(preds[[4]]$m$test * nrow(train)/train[nrow(train),1]))

# Compute likelihoods
mllk.trains <- lapply(preds,function(p) log(multinomial.score(p$m$train,train)))
llk.trains <- lapply(preds,function(p) RemLogLikelihoodVecFromArray(p$lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train)))

mllk.tests <- lapply(preds,function(p) log(multinomial.score(p$m$test,test)))
llk.tests <- lapply(preds,function(p) RemLogLikelihoodVecFromArray(p$lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test)))

names(mllk.trains) <- names(llk.trains) <- names(mllk.tests) <- names(llk.tests) <- c("baseline","beta.12 fixed","beta.12","beta.0.12","beta.0,1,2,12","beta.0")#,"beta.0,1,2,12 alt")

ans <- rbind(mllk.train=sapply(mllk.trains,mean),
      llk.train =sapply(llk.trains,mean),
      mllk.test =sapply(mllk.tests,mean),
      llk.test  =sapply(llk.tests,mean))


pdf("figs/onlinedebug.eckmann.pdf",width=12,height=10)
par(mfrow=c(2,4))
for (i in 2:5) { plot(mllk.trains[[1]],mllk.trains[[i]],asp=1,pch=".",xlab="baseline",ylab="model",main="mult. llk")
abline(0,1)
}
for (i in 2:5) {
  plot(llk.trains[[1]], llk.trains[[i]],asp=1,pch=".",xlab="baseline",ylab="model",main="brem llk")
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



fit.3 <- brem.mcmc(train[100:M,],N,K,s,model.type=opts$model.type,mh=!opts$slice,
                   niter=opts$numiterations,gibbs=opts$gibbs,beta=NULL,px=px,
                   outfile=NULL)
preds.3 <- get.pred(train,A,test.ix,fit.3)
total.lambda.3 <- sapply(1:nrow(train),function(i) sum(exp(preds.3$lrm$train[i,,])))
time.hat <- cumsum(1/total.lambda.3)
lines(time.hat,col="red")

# Make plot of rate of popular dyad over time
pdf("onlinedebug-75-85.pdf")
i <- 75; j <- 85
ymax <- max(c(exp(preds[[1]]$lrm$train[1:500,i,j]),
              exp(preds[[5]]$lrm$train[1:500,i,j])))
plot(exp(preds[[1]]$lrm$train[1:500,i,j]),type="S",ylab="rate",xlab="m",ylim=c(0,ymax))
v <- sapply(1:nrow(train), function(m) i %in% train[m,] | j %in% train[m,])
points(which(v),rep(.2,sum(v)),pch=17)
lines(exp(preds[[2]]$lrm$train[1:500,i,j]),type="S",col="red")
# lines(exp(preds[[3]]$lrm$train[1:500,i,j]),type="S",col="green")
# lines(exp(preds[[4]]$lrm$train[1:500,i,j]),type="S",col="blue")
dev.off()
