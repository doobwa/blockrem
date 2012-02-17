opts=list(dataset="eckmann-small",numclusters=1,model.type="full",gibbs=TRUE,numiterations=500)
px <- rep(0,13)
px[12] <- 1

library(brem)

load(paste("data/",opts$dataset,".rdata",sep=""))

# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()
priors <- list(beta=list(mu=0,sigma=1))

# Initialize with K=1 solution, if available
f <- paste("results/",opts$dataset,"/full.1.rdata",sep="")
if (K > 1 & file.exists(f)) {
  load(f)
  beta <- array(res$beta,c(P,K,K))
} else {
  beta <- array(0,c(P,K,K))
}

px <- rep(0,13)
px[8:12] <- 1
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,px=px,mh=FALSE,
                 niter=opts$numiterations,gibbs=opts$gibbs,beta=beta,
                 outdir=paste("results/",opts$dataset,"/",sep=""))
