opts <- list(dataset="synthetic-1",numclusters=2,model.type="full",numiterations=20)

load(paste("data/",opts$dataset,".rdata",sep=""))
library(brem)
source("pkg/R/brem.r")
source("pkg/R/splitmerge.r")
source("pkg/R/brem.alt.r")

# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
ego <- 1
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P,ego)
s$precompute()

# Fit model

source("pkg/R/brem.alt.r")


## Set priors
effects <- c("intercept","abba","abby","abay")

priors <- list(alpha=1,sigma=.1,phi=list(mu=0,sigma=1))
fit <- brem(train,N,K,effects)

train.ix <- 1:nrow(train)
p <- evaluate(A,N,train.ix,test.ix,fit,ties.method="random")
q <- evaluate.baseline(A,N,train.ix,test.ix,"online",ties.method="random")

par(mfrow=c(2,2))
plot(p$llk$train,q$llk$train,xlab="model",ylab="baseline"); abline(0,1)
plot(p$llk$test,q$llk$test,xlab="model",ylab="baseline"); abline(0,1)
plot(p$mllk$train,q$mllk$train,xlab="model",ylab="baseline"); abline(0,1)
plot(p$mllk$test,q$mllk$test,xlab="model",ylab="baseline"); abline(0,1)

# Profile model fitting
Rprof()
fit <- brem(train,N,K,effects)
Rprof(NULL)
summaryRprof()

# Investigate role of prior


#fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,px=px,do.sm=opts$splitmerge,num.extra=opts$numextra,niter=as.numeric(opts$numiterations),verbose=TRUE,outfile=outfile)
Rprof(NULL)
summaryRprof()

px <- rep(1,13)
px[13] <- 0
px[7] <- 0
#px[8]  <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,
                 outfile=NULL)
fit2 <- fit

# look at mixing for both
library(coda)
plot(mcmc(fit1$param[,,1,1]))
plot(mcmc(fit2$param[,,1,1]))


# Use larger K, learn z, learn beta.  Initialize with K=1 fit.
K <- 2
opts$gibbs <- "fast"
opts$slice <- TRUE
px <- rep(1,13)
px[8:13] <- 0
px[7] <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,
                 outfile=NULL)
fit1 <- fit

px <- rep(1,13)
px[13] <- 0
px[7] <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,#beta=beta,
                 outfile=NULL)
fit2 <- fit

# Use shared model
opts$model.type <- "shared"
px <- rep(1,13)
px[13] <- 0
px[7] <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,z=z,
                 outdir=NULL)
fit3 <- fit


# Examine trace plot for block (1,1)
plot(mcmc(fit1$param[,,1,1]))
plot(mcmc(fit2$param[,,1,1]))
plot(mcmc(fit3$param[,,1,1]))

# compare estimated rates for both at a given timepoint
load(paste("data/",opts$dataset,".rdata",sep=""))
lrm <- brem.lrm(train,N,z,beta)
lrm1 <- brem.lrm(train,N,fit1$z,fit1$beta)
lrm2 <- brem.lrm(train,N,fit2$z,fit2$beta)
m <- 50
plot(c(lrm1[m,,]),c(lrm2[m,,]),xlab="without degree effects",ylab="with degree effects")
abline(0,1)
plot(fit1$beta,fit2$beta)
for (i in 1:2000) {
  #lrm1[i,,] <- exp(lrm1[i,,])
  #lrm2[i,,] <- exp(lrm2[i,,])
  diag(lrm[i,,]) <- diag(lrm2[i,,]) <- diag(lrm1[i,,]) <- -5
}
plotmat(melt(exp(lrm[30,,])))
