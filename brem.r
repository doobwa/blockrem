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
              help="Slow or fast version of gibbs.")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

load(paste("data/",opts$dataset,".rdata",sep=""))

library(brem)

# Precompute data structures
N <- max(c(train[,2],train[,3]))
M <- nrow(train)
P <- 13
s <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,M,P)
s$precompute()

fit <- brem.mcmc(train,N,opts$numclusters,s,model.type=opts$model.type,
		 niter=opts$numiterations,gibbs=opts$gibbs,
                 outdir=paste("results/",opts$dataset,"/",sep=""))
