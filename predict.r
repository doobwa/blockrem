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

options(verbose=FALSE)

# Pull data from each saved file and grab the name of the fit
load(paste("data/",opts$dataset,".rdata",sep=""))


test.ix <- (1:nrow(A))[-(1:nrow(train))]

if (opts$model %in% c("online","uniform","marg")) {
  pred <- get.pred.baseline(train,A,test.ix,opts$model)
} else {
  if (opts$dataset == "synthetic" & opts$model == "truth") {
    fit <- list(beta=beta,z=z,llks=rep(true.lpost,500),niter=500,zs=rep(list(z),500),param=beta)  # true values
    save(fit,file=paste("results/",opts$dataset,"/",opts$model,".rdata",sep=""))
  }
  load(paste("results/",opts$dataset,"/",opts$model,".rdata",sep=""))
  pred <- get.pred(train,A,test.ix,fit)
}

llkm.train <- log(multinomial.score(pred$m$train,train))
llkm.test  <- log(multinomial.score(pred$m$test, test))

# Compute loglikelihood of each observation
llk.train <- RemLogLikelihoodVecFromArray(pred$lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
llk.test <- RemLogLikelihoodVecFromArray(pred$lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

cat("\nllkm.train:",mean(llkm.train))
cat("\nllkm.test:",mean(llkm.test))
cat("\nllk.train:",mean(llk.train))
cat("\nllk.test:",mean(llk.test))

# Compute rank of each observation
cat("\nranks (train)\n")
rk.train <- ranks(train,-pred$lrm$train,ties.method="random")
cat("ranks (test)\n")
rk.test  <- ranks(test, -pred$lrm$test, ties.method="random")

# Save
dir.create(paste("results/",opts$dataset,"/llks/",sep=""),showWarnings=FALSE)
save(llkm.train,llkm.test,llk.train,llk.test,opts,file=paste("results/",opts$dataset,"/llks/",opts$model,".rdata",sep=""))
dir.create(paste("results/",opts$dataset,"/ranks/",sep=""),showWarnings=FALSE)
save(rk.train,rk.test,opts,file=paste("results/",opts$dataset,"/ranks/",opts$model,".rdata",sep=""))


