#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-f", "--inputfile"), 
              help="Name of .rdata file containing an event history matrix named A."),
  make_option(c("-k", "--numclusters"), type="integer",default=2,
              help="Number of latent classes to use [default %default]."),
  make_option(c("-n", "--numiterations"), type="integer", default=100,
              help="Number of MCMC iterations [default %default]"),
  make_option(c("-m","--model.type"), default="full",
              help="Type of model to fit.  Options: \"baserate\", \"shared\", \"full\"."),
  make_option(c("-g","--gibbs"), default="fast",
              help="Slow or fast version of gibbs."),
  make_option(c("-d","--outdir"), help="Directory to save output from the model.")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

load(opts$inputfile)

library(brem)

# Precompute data structures
N <- max(c(A[,2],A[,3]))
M <- nrow(A)
P <- 11
s <- new(RemStat,A[,1],A[,2]-1,A[,3]-1,N,M,P)
s$precompute()

fit <- brem.mcmc(A,N,opts$numclusters,s,model.type=opts$model.type,
		             niter=opts$numiterations,outdir=opts$outdir,gibbs=opts$gibbs)

