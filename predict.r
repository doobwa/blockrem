#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-m", "--model"), default=NULL,
              help=""))
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

options(verbose=FALSE)

# Pull data from each saved file and grab the name of the fit
load(paste("data/",opts$dataset,".rdata",sep=""))

test.ix <- (1:nrow(A))[-(1:nrow(train))]


llkm.train <- log(multinomial.score(m.train,train))
llkm.test  <- log(multinomial.score(m.test, test))

# Compute loglikelihood of each observation
llk.train <- loglikelihood_vec_from_lrm(lrm.train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
llk.test <- loglikelihood_vec_from_lrm(lrm.test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

# Compute rank of each observation
cat("ranks (train)\n")
rk.train <- ranks(train,-lrm.train,ties.method="random")
cat("ranks (test)\n")
rk.test  <- ranks(test, -lrm.test, ties.method="random")

# Save
dir.create(paste("results/",opts$dataset,"/llks/",sep=""))
save(llkm.train,llkm.test,llk.train,llk.test,opts,file=paste("results/",opts$dataset,"/llks/",opts$model,".rdata",sep=""))
dir.create(paste("results/",opts$dataset,"/ranks/",sep=""))
save(rk.train,rk.test,opts,file=paste("results/",opts$dataset,"/ranks/",opts$model,".rdata",sep=""))


