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

load(paste("data/",opts$dataset,".rdata",sep=""))

# Get indices for events that had dyads which were previously unobserved
ctrain <- counts.of.observed(train,N)
ctest <- counts.of.observed(test,N)

save(ctrain,ctest,opts,file=paste("results/",opts$dataset,"/counts.rdata",sep=""))
