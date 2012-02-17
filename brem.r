#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of /data/[dataset].rdata file containing an event history matrix named train.  Saves results to /results/[dataset]/."),
  make_option(c("-k", "--numclusters"), type="integer",default=2,
              help="Number of latent classes to use [default %default]."),
  make_option(c("-n", "--numiterations"), type="integer", default=100,
              help="Number of MCMC iterations [default %default]"),
  make_option(c("-m","--model.type"), default="full",
              help="Type of model to fit.  Options: \"baserate\", \"shared\", \"full\"."),
  make_option(c("-g","--gibbs"), default="fast",
              help="Slow or fast version of gibbs."),
  make_option(c("-s","--slice"), default=TRUE,
              help="Slice sample instead of MH.")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

load(paste("data/",opts$dataset,".rdata",sep=""))

library(brem)

# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
K <- opts$numclusters
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

# Initialize with K=1 solution, if available
f <- paste("results/",opts$dataset,"/full.1.rdata",sep="")
if (K > 1 & file.exists(f)) {
  load(f)
  beta <- array(res$beta,c(P,K,K))
} else {
  beta <- NULL
}

source("pkg/R/brem.r")
fit <- brem.mcmc(train,N,K,s,model.type=opts$model.type,mh=!opts$slice,
		 niter=opts$numiterations,gibbs=opts$gibbs,beta=beta,
                 outdir=paste("results/",opts$dataset,"/",sep=""))
