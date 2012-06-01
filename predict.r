#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option("--force", default=FALSE),
  make_option("--baselines", default=FALSE),
  make_option("--niters", type="integer", default=10)  )
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#
library(brem);library(reshape2);opts <- list(dataset="eckmann-small",baselines=FALSE,force=TRUE,niters=1)
source("utils.r")

options(verbose=FALSE)

dataset <- opts$dataset
load(paste("data/",dataset,".rdata",sep=""))
train.ix <- 1:nrow(train)
test.ix <- (1:nrow(A))[-train.ix]
file.remove("progress.rdata")
dir.create(paste("results/",dataset,"/preds/",sep=""),showWarnings=FALSE)

for (model in c("online","uniform","marg")) {
  progress(dataset,model,"started")
  pred <- evaluate.baseline(A,N,train.ix,test.ix,model)
  save(pred,opts,file=paste("results/",dataset,"/preds/",model,".rdata",sep=""))
  progress(dataset,model,"completed")  
}

folder <- paste("results/",opts$dataset,"/fits",sep="")
models <- modelnames(folder)
for (model in models) {
  progress(dataset,model,"started")
  f <- paste(folder,"/",model,".rdata",sep="")
  load(f)
  pred <- evaluate(A,N,train.ix,test.ix,fit,niters=as.numeric(opts$niters))
  save.pred(pred,dataset,model)
  progress(dataset,model,"completed")
}

if (opts$dataset == "synthetic-1") {
  niter <- 500
  fit <- list(samples=list(list(beta=beta,z=z)),llks=rep(true.lpost,niter),niter=niter,zs=rep(list(z),iter=1),params=list(beta=beta,z=z),ego=1,transform=TRUE)  # true values
  pred <- evaluate(A,N,train.ix,test.ix,fit,niters=NULL)
  save.pred(pred,dataset,"truth")
}

