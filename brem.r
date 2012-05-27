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
  make_option(c("-m", "--model.type"), 
              help="Perform split merge moves."),
  make_option(c("-e", "--numextra"), type="integer", default=2,
              help="Number of extra clusters to sample from the prior")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))
#library(brem); opts <- list(dataset="twitter-sm",numcluster=20,numiterations=100,splitmerge=FALSE,numextra=5,model.type="full")

load(paste("data/",opts$dataset,".rdata",sep=""))

# Precompute data structures
# N should be loaded by dataset
P <- 13
K <- opts$numclusters

outfile <- paste("results/",opts$dataset,"/",opts$model.type,".",K,".rdata",sep="")
effects <- c("intercept","abba","abby","abay","dyad_count")
priors <- list(alpha=1,sigma.proposal=.1,phi=list(mu=0,sigma=1),mu=list(mu=0,sigma=1),sigma=list(alpha=3,beta=1))

# Fit and save model
fit <- brem(train,N,K,effects,outfile=outfile)
save(fit,outfile)
