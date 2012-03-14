#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("brem"))
source("pkg/R/brem.r")

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of /data/[dataset].rdata file containing an event history matrix named train.  Saves results to /results/[dataset]/."),
  make_option(c("-k", "--numclusters"), type="integer",default=2,
              help="Number of latent classes to use [default %default]."),
  make_option(c("-n", "--numiterations"), type="integer", default=100,
              help="Number of MCMC iterations [default %default]"),
  make_option(c("-t","--model.type"), default="full",
              help="Type of model to fit.  Options: \"baserate\", \"shared\", \"full\"."),
  make_option(c("-g","--gibbs"), default=TRUE,
              help="Slow or fast version of gibbs."),
  make_option(c("-m","--mh"), default=FALSE,
              help="Slow or fast version of gibbs."),
  make_option(c("-s","--slice"), default=FALSE,
              help="Slice sample instead of MH."),
  make_option(c("-i","--initialize"), default=FALSE,
              help="initialize with K=1 solution"),
  make_option(c("-f","--fixz"), default=FALSE,
              help="if twitter, fix z to predetermined values")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

load(paste("data/",opts$dataset,".rdata",sep=""))


# Precompute data structures
# N should be loaded by dataset
M <- nrow(train)
P <- 13
K <- opts$numclusters

# Initialize with K=1 solution, if available
if (K==1) opts$gibbs <- FALSE

f <- paste("results/",opts$dataset,"/full.1.rdata",sep="")
if ((K > 1) & file.exists(f) & opts$initialize) {
 load(f)
 beta <- array(res$beta,c(P,K,K))
} else {
 beta <- NULL
}

px <- rep(1,13)
px[13] <- 0
px[7]  <- 0

if (opts$dataset=="synthetic") {
  px[8:13] <- 0
}

if (opts$dataset=="twitter-small" & opts$fixz) {
  tb <- table(factor(c(train[,2],train[,3]),1:N))
  chosen <- names(tb)[which(tb > 20)]
  chosen <- as.numeric(chosen)
  z <- rep(1,N)
  z[chosen] <- 2
  opts$gibbs <- FALSE
} else {
  z <- NULL
}

outfile <- paste("results/",opts$dataset,"/",opts$model.type,".",K,".rdata",sep="")

priors <- list(beta=list(mu=0,sigma=1))

fit <- brem.mcmc(train,N,K,model.type=opts$model.type,slice=opts$slice,gibbs=opts$gibbs,mh=opts$mh,
                 niter=opts$numiterations,beta=beta,px=px,z=z,m=10,
                 outfile=outfile,priors=priors,skip.intercept=FALSE)
