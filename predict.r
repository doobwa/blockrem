#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-m", "--model"), default=NULL,
              help="if not specified, baseline will be computed"),
  make_option(c("-b", "--baseline"), default="uniform",
              help="uniform, online, marg"))
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

options(verbose=FALSE)

# Pull data from each saved file and grab the name of the fit
load(paste("data/",opts$dataset,".rdata",sep=""))

if (!is.null(opts$model)) {
  f <- paste("results/",opts$dataset,"/",opts$model,".rdata",sep="")
  load(f)
  cat("precomputing\n")
  strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
  strain$precompute()
  stest <- new(RemStat,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test),P)
  stest$precompute()
  cat("lambdas (train)\n")
  lrm.train <- brem.lrm.fast(nrow(train), strain, res$z, res$beta)
  cat("lambdas (test)\n")
  lrm.test  <- brem.lrm.fast(nrow(test),   stest, res$z, res$beta)
} else {
  get.lrm <- switch(opts$baseline,
                    "uniform" = function(x) array(1,c(nrow(x),N,N)),
                    "online"  = function(x) ratemat.online(x,N),
                    "marg" = function(x) ratemat.from.marginals(train,x,N))
  cat("lambdas (train)\n")
  lrm.train <- get.lrm(train)
  cat("lambdas (test)\n")
  lrm.test  <- get.lrm(test)
}

cat("ranks (train)\n")
rk.train <- ranks(train,-lrm.train,ties.method="random")
cat("ranks (test)\n")
rk.test  <- ranks(test, -lrm.test, ties.method="random")

dir.create(paste("results/",opts$dataset,"/ranks/",sep=""))
save(rk.train,rk.test,opts,file=paste("results/",opts$dataset,"/ranks/",opts$model,".rdata",sep=""))
