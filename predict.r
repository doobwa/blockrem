#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-t", "--model"), default=NULL,
              help=""))
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#opts <- list(dataset="synthetic-1",model="full")

options(verbose=FALSE)

# Pull data from each saved file and grab the name of the fit
load(paste("data/",opts$dataset,".rdata",sep=""))
if (opts$dataset == "synthetic" & opts$model == "truth") {
  niter <- 500
  fit <- list(beta=beta,z=z,llks=rep(true.lpost,niter),niter=niter,zs=rep(list(z),niter),param=beta,ego=1)  # true values
  save(fit,file=paste("results/",opts$dataset,"/",opts$model,".rdata",sep=""))
}

train.ix <- 1:nrow(train)
test.ix <- (1:nrow(A))[-(1:nrow(train))]

if (opts$model %in% c("online","uniform","marg")) {
  pred <- evaluate.baseline(A,N,train.ix,test.ix,opts$model)
} else {
  load(paste("results/",opts$dataset,"/",opts$model,".rdata",sep=""))
  pred <- evaluate(A,N,train.ix,test.ix,fit)
}

# Temporary: to fit in with rest of pipeline
llkm.train <- pred$mllk$train
llkm.test <-  pred$mllk$test
llk.train <- pred$llk$train
llk.test  <- pred$llk$test
rk.train  <- pred$rks$train
rk.test <- pred$rks$test

# Save
dir.create(paste("results/",opts$dataset,"/llks/",sep=""),showWarnings=FALSE)
save(llkm.train,llkm.test,llk.train,llk.test,opts,file=paste("results/",opts$dataset,"/llks/",opts$model,".rdata",sep=""))
dir.create(paste("results/",opts$dataset,"/ranks/",sep=""),showWarnings=FALSE)
save(rk.train,rk.test,opts,file=paste("results/",opts$dataset,"/ranks/",opts$model,".rdata",sep=""))


