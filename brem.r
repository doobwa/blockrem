#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("brem"))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of /data/[dataset].rdata file containing an event history matrix named train.  Saves results to /results/[dataset]/."),
  make_option(c("-k", "--numclusters"), type="integer",default=2,
              help="Number of latent classes to use [default %default]."),
  make_option(c("-n", "--numiterations"), type="integer", default=100,
              help="Number of MCMC iterations [default %default]"),
  make_option(c("-s", "--splitmerge"), default=FALSE,
              help="Perform split merge moves."),
  make_option(c("-g", "--degrees"), default=FALSE,
              help="Include degree effects."),
  make_option(c("-t", "--transform"), default=TRUE,
              help="Transform degree effects using log"),
  make_option(c("-b", "--negbinom"), default=FALSE,
              help="Use negative binomial prior instead of DP"),
  make_option(c("-e", "--numextra"), type="integer", default=2,
              help="Number of extra clusters to sample from the prior")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))
#library(brem); opts <- list(dataset="twitter-small",numclusters=20,numiterations=100,splitmerge=FALSE,numextra=5,model.type="full")
#library(brem); opts <- list(dataset="realitymining-small",numclusters=20,numiterations=500,splitmerge=FALSE,numextra=5,model.type="full")
#library(brem); opts <- list(dataset="eckmann-small",numclusters=10,numiterations=500,splitmerge=FALSE,numextra=5,negbinom=FALSE,degrees=FALSE)

load(paste("data/",opts$dataset,".rdata",sep=""))

# Precompute data structures
# N should be loaded by dataset
P <- 13
opts$model.type <- paste("kinit",opts$numclusters,".sm",opts$splitmerge*1,".nb",opts$negbinom*1,".deg",opts$degrees*1,".trans",opts$transform*1,sep="")

dir.create(paste("results/",opts$dataset,sep=""),showWarn=FALSE)
dir.create(paste("results/",opts$dataset,"/fits/",sep=""),showWarn=FALSE)
outfile <- paste("results/",opts$dataset,"/fits/",opts$model.type,".rdata",sep="")

if (opts$degrees) {
  effects <- c("intercept","abba","abby","abay","sen_outdeg","sen_indeg","dyad_count")
} else {
  effects <- c("intercept","abba","abby","abay")
}

priors <- list(alpha=.1,sigma.proposal=.1,phi=list(mu=0,sigma=1),mu=list(mu=0,sigma=1),sigma=list(alpha=2,beta=2))

K <- opts$numclusters
if (opts$negbinom) priors$negbinom <- list(shape = 3, mu=N/K)

# Fit and save model
fit <- brem(train,N,priors,K,effects,transform=opts$transform,do.sm=opts$splitmerge,num.extra=opts$numextra,niter=opts$numiterations,outfile=outfile)
save(fit,file=outfile)
