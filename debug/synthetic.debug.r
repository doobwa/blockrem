library(brem)
opts=list(dataset="synthetic",numclusters=2,model.type="full",gibbs=TRUE,numiterations=100,slice=TRUE,initialize=FALSE,fixz=FALSE,skip.intercept=FALSE)
priors <- list(beta=list(mu=0,sigma=1))
outfile <- NULL
beta <- NULL
M <- nrow(train)
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

beta[c(2,3,4,5,6,8,9,10,11,12),,] <- rnorm(10*K*K,0,3)
beta[c(1,2,3,4,5,6),,] <- rnorm(6*K*K,0,3)
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=TRUE,slice=TRUE,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=beta,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)


outfile <- paste("results/",opts$dataset,"/",opts$model.type,".",K,".rdata",sep="")

# Fixing degree effects and intercepts works
b <- beta
px <- rep(1,13)
px[7:13] <- 0; px[1] <- 0
b[2:6,,] <- rnorm(5*K*K,0,3)

fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=TRUE,slice=TRUE,
                 niter=opts$numiterations,gibbs=FALSE,beta=b,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)


# Setting all  degree effects to 0
b <- beta
px <- rep(1,13)
px[7:13] <- 0
b[1:6,,] <- rnorm(6*K*K,0,1)

fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=TRUE,slice=TRUE,
                 niter=opts$numiterations,gibbs=FALSE,beta=b,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)

# Learning all but intercepts
chosen <- c(2,3,4,5,8,9,10,11,12)
b <- beta
px <- rep(1,13); px[-chosen] <- 0
b[chosen,,] <- rnorm(length(chosen)*K*K,0,1)
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=b,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)

# Learning intercepts and pshifts
b <- beta
px <- rep(1,13); px[7:13] <- 0
b[c(1:6),,] <- rnorm(6*K*K,0,1)
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=b,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)

# Just intercepts
b <- beta
px <- rep(0,13); px[1] <- 1
b[c(1),,] <- rnorm(1*K*K,0,1)
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=opts$mh,slice=opts$slice,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=b,px=px,z=z,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)

# Get statistics for actors 1-5 


#load("results/synthetic/full.2.rdata")
load(paste("data/",opts$dataset,".rdata",sep=""))

# Check llk functions
M <- nrow(train)
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
s$precompute()
tmp1 <- RemLogLikelihoodPc(res$beta, res$z - 1, s$ptr(), K)
tmp2 <- RemLogLikelihood(res$beta, train[,1],train[,2]-1,train[,3]-1,res$z-1,N,nrow(train),K,P)
sum(tmp2)
sum(brem.llk(train,N,res$z,res$beta))
brem.lpost(train,N,K,res$z,res$beta,priors)

# Look at logposterior as a function of one of the intercepts.  Is it flat?
b <- res$beta
xs <- seq(-3,3,by=.1)
lpost.int11 <- sapply(xs,function(x) {
  b[1,1,1] <- x
  brem.lpost(train,N,K,z,b,priors)
})
plot(xs,lpost.int11,type="l")

# For true model
b <- beta
xs <- seq(-2,2,by=.1)
lpost.int11 <- sapply(xs,function(x) {
  b[1,1,1] <- x
  brem.lpost(train,N,K,res$z,b,priors)
})
plot(xs,lpost.int11,type="l")

# Logposterior as we vary (intercept,pshifts) all at once.  Flat?
b <- res$beta
xs <- seq(-1,1,by=.1)
lpost.int11 <- sapply(xs,function(x) {
  b[2:6,1,1] <- b[2:6,1,1] + x
  b[1,1,1] <- b[1,1,1] - x
  brem.lpost(train,N,K,res$z,b,priors)
})
plot(xs,lpost.int11,type="l")

# Twodims
xs <- seq(-.5,.5,by=.05)
df <- expand.grid(x1=xs,x2=xs)
lpost.int11 <- ddply(df,.(x1,x2),function(x) {
  b <- res$beta
  b[1,1,1] <- b[1,1,1] + x$x1
  b[3,1,1] <- b[3,1,1] + x$x2
  brem.lpost(train,N,K,res$z,b,priors)
})
ggplot(subset(lpost.int11,x1 < .5 & x1 > -.5 & x2 < .5 & x2 > -.5), aes(x=x1, y=x2, z = V1)) + stat_contour() 