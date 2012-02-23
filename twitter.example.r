# Get estimates for full 2 when we manually set the groups
library(brem)
load("data/twitter-small.rdata")
source("pkg/R/brem.r")

tb <- table(factor(c(A[,2],A[,3]),1:N))
chosen <- names(tb)[which(tb > 20)]
chosen <- as.numeric(chosen)
z <- rep(1,N)
z[chosen] <- 2

opts <- list(slice=FALSE,gibbs=FALSE,model.type="full",dataset="twitter-small",numiterations=500,numclusters=2,initialize=TRUE)

# Precompute data structures
M <- nrow(train)
P <- 13
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

px <- rep(1,13)
px[13] <- 0
px[7]  <- 0
f <- paste("results/",opts$dataset,"/full.1.rdata",sep="")
if (K > 1 & file.exists(f) & opts$initialize) {
  load(f)
  beta <- array(res$beta,c(P,K,K))
} else {
  beta <- NULL
}

# fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice, mcmc.sd=.05,m=5,
#                  niter=opts$numiterations, gibbs=opts$gibbs,  beta=beta, px=px, z=z,   
#                  outfile="tmp.slice.m5.rdata")
# 

opts$gibbs <- "fast"
opts$slice <- TRUE
# fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice, mcmc.sd=.05,m=20,
#                  niter=opts$numiterations, gibbs=opts$gibbs,  beta=beta, px=px, z=z,   
#                  outfile="tmp.slice.gibbs.rdata")


fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice, mcmc.sd=.05,m=20,
                 niter=opts$numiterations, gibbs=opts$gibbs,  beta=beta, px=px,  
                 outfile="tmp.slice.gibbs.learn2.rdata")

# Visualzie model paramters

load("tmp.slice.rdata")
r <- melt(res$param)
r <- subset(r,X1<10)
q7 <- qplot(X1,value,data=r, colour=factor(X2),geom="line") + labs(colour="parameters for\n 1x1 block",x="iteration") + theme_bw() + facet_grid(X3~X4,scales="free")
q7
load("tmp.slice.gibbs.rdata")
r <- melt(res$param)
r <- subset(r,X1<10)
q7 <- qplot(X1,value,data=r, colour=factor(X2),geom="line") + labs(colour="parameters for\n 1x1 block",x="iteration") + theme_bw() + facet_grid(X3~X4,scales="free")
q7

load("tmp.slice.gibbs.learn2.rdata")
r <- melt(res$param)
r <- subset(r,X1<10)
q7 <- qplot(X1,value,data=r, colour=factor(X2),geom="line") + labs(colour="parameters for\n 1x1 block",x="iteration") + theme_bw() + facet_grid(X3~X4,scales="free")
q7

# Look at whether learned betas are statistically significant
b <- melt(fit.lrn$param)

df <- data.frame(count=tb,z=z)
qplot(log(count.Freq),data=df,geom="histogram") +facet_grid(z~.)

load("results/twitter-small/full.3.rdata")
M <- nrow(train)
P <- 13
K <- 3
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()
beta <- res$beta
brem.lpost.fast(train,N,K,res$z,s,res$beta)

beta[c(8,11),3,1] <- 10
brem.lpost.fast(train,N,K,res$z,s,beta)
