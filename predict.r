#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(multicore))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option("--niters", type="integer", default=10)  )
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#
#library(brem);library(reshape2);library(multicore);opts <- list(dataset="twitter-small",niters=10)
source("utils.r")

options(verbose=FALSE)

dataset <- opts$dataset
load(paste("data/",dataset,".rdata",sep=""))
train.ix <- 1:nrow(train)
test.ix <- (1:nrow(A))[-train.ix]

dir.create(paste("results/",dataset,"/preds/",sep=""),showWarnings=FALSE)
pf <- "progress.rdata"


folder <- paste("results/",opts$dataset,"/fits",sep="")
models <- modelnames(folder)
options(cores=8)

# TEMP
models <- c("kinit2.kmax1.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha5000.xsigbeta1000",
            "kinit2.kmax2.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha5000.xsigbeta1000",
            "kinit2.kmax3.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha5000.xsigbeta1000",
            "kinit2.kmax10.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha5000.xsigbeta1000",
            "kinit2.kmax10.sm0.nb0.pshift0.deg0.trans1.collapse1.xsigalpha5000.xsigbeta1000")
for (model in models) progress(pf,c(dataset,model,"started",""))

x <- mclapply(models,function(model) {
#for (model in models) {
  st <- proc.time()
  f <- paste(folder,"/",model,".rdata",sep="")
  if (!file.exists(f)) {
    progress(pf,c(dataset,model,"failed","no model file"))
  } else {
    load(f)
    if (opts$niters >= fit$iter) {
      progress(pf,c(dataset,model,"failed","not enough samples"))
    } else {
      pred <- evaluate(A,N,train.ix,test.ix,fit,niters=as.numeric(opts$niters))
      save(pred,opts,file=paste("results/",dataset,"/preds/",model,".rdata",sep=""))
      progress(pf,c(dataset,model,"complete",proc.time()[3] - st[3]))
    }
  }
  return(NULL)
})

models <- c("online","uniform","marg")
x <- mclapply(models,function(model) {
#for (model in c("online","uniform","marg")) {
  st <- proc.time()  
  progress(pf,c(dataset,model,"started",""))
  pred <- evaluate.baseline(A,N,train.ix,test.ix,model)
  save(pred,opts,file=paste("results/",dataset,"/preds/",model,".rdata",sep=""))
  progress(pf,c(dataset,model,"completed",proc.time()[3] - st[3]))
  return(NULL)
})


if (opts$dataset == "synthetic-1") {
  niter <- 500
  fit <- list(samples=list(list(beta=beta,z=z)),llks=rep(true.lpost,niter),niter=niter,zs=rep(list(z),iter=1),params=list(beta=beta,z=z),ego=1,transform=TRUE)  # true values
  pred <- evaluate(A,N,train.ix,test.ix,fit,niters=NULL)
  save.pred(pred,dataset,"truth")
}

