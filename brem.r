#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("brem"))
suppressPackageStartupMessages(library("MCMCpack"))

option_list <- list(
  make_option("--dataset", 
              help="Name of /data/[dataset].rdata file containing an event history matrix named train.  Saves results to /results/[dataset]/."),
  make_option("--kinit", type="integer",default=2,
              help="Number of initial latent classes [default %default]."),
  make_option("--kmax", type="integer",default=3,
              help="Maximum number of latent classes [default %default]"),
  make_option("--numiterations", type="integer", default=100,
              help="Number of MCMC iterations [default %default]"),
  make_option("--splitmerge", default=FALSE,
              help="Perform split merge moves."),
  make_option("--pshifts", default=FALSE,
              help="Include pshift effects."),
  make_option("--degrees", default=FALSE,
              help="Include degree effects."),
  make_option("--transform", default=TRUE,
              help="Transform degree effects using log"),
  make_option("--type", default="crp",
              help=""),
  make_option("--negbinom", default=FALSE,
              help="Use negative binomial prior instead of CRP"),
  make_option("--collapse", default=FALSE,
              help="Integrate out sigma when sampling from beta"),
  make_option("--crp.hyper", type="numeric",default=.1,
              help="Sigma hyperparameter"),
  make_option("--sigma.hyper.alpha", type="numeric",default=1,
              help="Sigma hyperparameter"),
  make_option("--sigma.hyper.beta", type="numeric",default=1,
              help="Sigma hyperparameter"),
  make_option("--force", default=FALSE,
              help="Fit model even if rdata exists"),
  make_option("--numextra", type="integer", default=2,
              help="Number of extra clusters to sample from the prior")
  )
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#library(MCMCpack);library(brem); opts <- list(dataset="classroom-16",numiterations=500,splitmerge=FALSE,numextra=5,negbinom=FALSE,degrees=FALSE,transform=TRUE,collapse=TRUE,numextra="5",force=TRUE,pshifts=TRUE,degrees=TRUE,kinit=2,kmax=3,crp.hyper=1,sigma.hyper=.1,type="crp");outfile=NULL

load(paste("data/",opts$dataset,".rdata",sep=""))

# Precompute data structures
# N should be loaded by dataset
P <- 13
opts$model.type <- paste("kinit",opts$kinit,".kmax",opts$kmax,".sm",opts$splitmerge*1,".nb",opts$negbinom*1,".pshift",opts$pshifts*1,".deg",opts$degrees*1,".trans",opts$transform*1,".collapse",opts$collapse*1,".xsigalpha",opts$sigma.hyper.alpha*1000,".xsigbeta",opts$sigma.hyper.beta*1000,sep="")#".crphyper",opts$crp.hyper,".sigmahyper",opts$sigma.hyper,sep="")
dir.create(paste("results/",opts$dataset,sep=""),showWarn=FALSE)
dir.create(paste("results/",opts$dataset,"/fits/",sep=""),showWarn=FALSE)
outfile <- paste("results/",opts$dataset,"/fits/",opts$model.type,".rdata",sep="")


if (file.exists(outfile) & !opts$force) {
  load(outfile)
  if (fit$iter == fit$niter) stop(paste("Already fit this:",outfile))
}

effects <- c("intercept")
if (opts$pshifts) {
  effects <- c(effects,"abba","abby","abay")
}
if (opts$degrees) {
  effects <- c(effects,"sen_outdeg","sen_indeg","dyad_count")
}

opts$sigma.hyper <- as.numeric(opts$sigma.hyper)
priors <- list(alpha=as.numeric(opts$crp.hyper),type=opts$type,
               sigma=list(alpha=opts$sigma.hyper.alpha,beta=opts$sigma.hyper.beta),
               sigma.proposal=.1,phi=list(mu=0,sigma=1),mu=list(mu=0,sigma=1))

kinit <- as.numeric(opts$kinit)
kmax  <- as.numeric(opts$kmax)
if (kmax < kinit) kinit <- kmax
if (opts$negbinom) priors$negbinom <- list(shape = 3, mu=N/kinit)

# Fit and save model
fit <- brem(train,N,priors,kinit=kinit,kmax=kmax,effects=effects,transform=opts$transform,do.sm=opts$splitmerge,num.extra=opts$numextra,niter=opts$numiterations,collapse.sigma=opts$collapse,outfile=outfile)
save(fit,file=outfile)
