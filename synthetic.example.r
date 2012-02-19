opts <- list(dataset="synthetic",numclusters=1,model.type="full",gibbs="fast",numiterations=20,slice=TRUE)

load(paste("data/",opts$dataset,".rdata",sep=""))
library(brem)

# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

# Degree vs. no degree effects, slice sampling, K=1 full

px <- rep(1,13)
px[8:13] <- 0
px[7] <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,
                 outdir=NULL)
fit1 <- fit

px <- rep(1,13)
px[13] <- 0
px[7] <- 0
#px[8]  <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,
                 outdir=NULL)
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
                 outdir=NULL)
fit1 <- fit

px <- rep(1,13)
px[13] <- 0
px[7] <- 0
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,px=px,#beta=beta,
                 outdir=NULL)
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
