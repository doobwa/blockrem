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

## Set priors
priors <- list(alpha=1,sigma=.1,phi=list(mu=0,sigma=1))
px <- rep(1,13)
px[8:13] <- 0
px[7] <- 0

## Test basic likelihood functions
K <- 2
k1 <- k2 <- 1
current <- array(rnorm(P*K^2,priors$phi$mu,priors$phi$sigma),c(P,K,K))
olp <- brem.lpost.fast(A,N,K,z,s,beta)

priors <- list(alpha=1,sigma=.1,phi=list(mu=0,sigma=1))
llk <- sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
llk <- sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,1:5))
llk <- sum(RemLogLikelihoodActorPc(3,beta,z-1,s$ptr(),K))
k1nodes <- which(z==k1)
k2nodes <- which(z==k2)
llk <- sum(RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,beta,z-1,s$ptr(),K))

# Fit model
source("pkg/R/brem.alt.r")
lp(beta,z,priors)
px <- c(1,2,3,4,5,6)
K <- 5; P <- 13
fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,do.sm=FALSE,num.extra=10,niter=15,verbose=TRUE)
lp(fit$samples[[5]]$phi,fit$samples[[5]]$z,priors)

fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,10,do.sm=TRUE,num.extra=10,niter=20,verbose=TRUE)

brem <- function(train,N,K=2,effects=c("intercept","abba","abby","abay"),ego=TRUE,do.sm=FALSE,num.extra=2,niter=20,verbose=TRUE) {
  M <- nrow(train)
  P <- 13
  ego <- ego*1  # RemStat doesn't want boolean
  s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P,ego)
  s$precompute()
  px <- which(effects %in% c("intercept","abba","abby","abxa","abxb","abay","abab","sen_outdeg","rec_outdeg","sen_indeg","rec_indeg","dyad_count","changepoint_count"))
  fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,px=px,do.sm=do.sm,num.extra=num.extra,niter=niter,verbose=TRUE)

  fit$ego <- ego
  fit$beta <- fit$samples[[niter]]$phi
  fit$z <- fit$samples[[niter]]$z
  class(fit) <- "brem"
  return(fit)
}

## Set priors
effects <- c("intercept","abba","abby","abay")
priors <- list(alpha=1,sigma=.1,phi=list(mu=0,sigma=1))
fit <- brem(train,N,K,effects)

train.ix <- 1:nrow(train)
p <- eval.online(A,N,train.ix,test.ix,fit,ties.method="random")
q <- eval.online.baseline(A,N,train.ix,test.ix,"online",ties.method="random")



# Profile model fitting
Rprof()
fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,px=c(1,2,3,4,5,6),do.sm=FALSE,num.extra=2,niter=20,verbose=TRUE)
Rprof(NULL)
summaryRprof()

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
