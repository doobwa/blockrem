# Get estimates for full 2 when we manually set the groups

load("data/twitter-small.rdata")

tb <- table(factor(c(train[,2],train[,3]),1:N))
chosen <- names(tb)[which(tb > 20)]
chosen <- as.numeric(chosen)

z <- rep(1,N)
z[chosen] <- 2


opts <- list(slice=FALSE,gibbs=FALSE,model.type="full",dataset="twitter-small",numiterations=500,numclusters=2)


# Precompute data structures
M <- nrow(train)
P <- 13
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

px <- rep(1,13)
px[13] <- 0
beta <- NULL

fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations, gibbs=opts$gibbs,  beta=beta, px=px, z=z,   
                 outdir=paste("tmp/",sep=""))

load("results/twitter-small/full.2.fixed.rdata")
z.fixed <- res$z
load("results/twitter-small/full.2.rdata")
z <- res$z
table(z,z.fixed)

# Start with the z's that were learned by full.3 and see if we get more distinct estimates
load("results/twitter-small/full.3.rdata")
z <- res$z
K <- 3

fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
                 niter=opts$numiterations, gibbs=FALSE,  beta=beta, px=px, z=z,   
                 outdir=paste("results/",opts$dataset,"/tmp",sep=""))


df <- data.frame(count=tb,z=z)
qplot(log(count.Freq),data=df,geom="histogram") +facet_grid(z~.)