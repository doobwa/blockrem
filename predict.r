#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."))
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#library(brem);opts <- list(dataset="eckmann-small")

options(verbose=FALSE)

dataset <- opts$dataset
load(paste("data/",dataset,".rdata",sep=""))

train.ix <- 1:nrow(train)
test.ix <- (1:nrow(A))[-(1:nrow(train))]

## Temporary: to fit in with rest of pipeline
save.pred <- function(pred,dataset,model) {
  llkm.train <- pred$mllk$train
  llkm.test <-  pred$mllk$test
  llk.train <- pred$llk$train
  llk.test  <- pred$llk$test
  rk.train  <- pred$rks$train
  rk.test <- pred$rks$test

  dir.create(paste("results/",dataset,"/llks/",sep=""),showWarnings=FALSE)
  save(llkm.train,llkm.test,llk.train,llk.test,opts,file=paste("results/",dataset,"/llks/",model,".rdata",sep=""))
  dir.create(paste("results/",dataset,"/ranks/",sep=""),showWarnings=FALSE)
  save(rk.train,rk.test,opts,file=paste("results/",dataset,"/ranks/",model,".rdata",sep=""))
}

filenames <- function(folder) {
  as.vector(sapply(dir(folder),function(x) strsplit(x,".rdata")[[1]][1]))
}

for (model in c("online","uniform","marg")) {
  cat(model,"\n")
  pred <- evaluate.baseline(A,N,train.ix,test.ix,model,ties.method="random")
  save.pred(pred,dataset,model)
}

results.dir <- paste("results/",opts$dataset,"/fits",sep="")
models <- filenames(results.dir)
for (model in models) {
  cat(model,"\n")
  load(paste(results.dir,"/",model,".rdata",sep=""))
  pred <- evaluate(A,N,train.ix,test.ix,fit,niters=NULL,ties.method="random")
  save.pred(pred,dataset,model)
}

if (dataset == "synthetic-1") {
  niter <- 500
  fit <- list(params=list(beta=beta,z=z),llks=rep(true.lpost,niter),niter=niter,zs=rep(list(z),niter),param=beta,ego=1)  # true values
  pred <- evaluate(A,N,train.ix,test.ix,fit)
  save.pred(pred,dataset,"truth")
}

# Load all loglikelihoods
  load(paste("data/",opts$dataset,".rdata",sep=""))
  folder <- paste("results/",opts$dataset,"/llks",sep="")
  fs <- dir(folder,full.names=TRUE)
  llks.test <- lapply(fs,function(f) {
    load(f)
    return(llk.test)
  })
  llks.train <- lapply(fs,function(f) {
    load(f)
    return(llk.train)
  })
  mllks.test <- lapply(fs,function(f) {
    load(f)
    return(llkm.test)
  })
  mllks.train <- lapply(fs,function(f) {
    load(f)
    return(llkm.train)
  })
  models <- filenames(folder)
  names(llks.test) <- names(llks.train) <- names(mllks.test) <- names(mllks.train) <- models
  
  llks.train <- melt(llks.train)
  llks.train$event <- 1:nrow(train)
  llks.test <- melt(llks.test)
  llks.test$event <- 1:nrow(test)
  mllks.train <- melt(mllks.train)
  mllks.train$event <- 1:nrow(train)
  mllks.test <- melt(mllks.test)
  mllks.test$event <- 1:nrow(test)
  
  # Make results table
library(plyr)
  df <- rbind(
        cbind(likelihood="rem",type="train",
              ddply(llks.train,.(L1),summarise,value=mean(value))),
        cbind(likelihood="rem",type="test",
              ddply(llks.test,.(L1),summarise,value=mean(value))),
        cbind(likelihood="mult",type="train",
              ddply(mllks.train,.(L1),summarise,value=mean(value))),
        cbind(likelihood="mult",type="test",
              ddply(mllks.test,.(L1),summarise,value=mean(value)))
              )
 # TODO: Fix this up 
#  df$L1 <- factor(df$L1,c("uniform","marg","online","full.1","full.2","full.3","dp","truth"))
  df$dataset <- dataset

# Get mean and sd of K over each fit
folder <- paste("results/",dataset,"/fits/",sep="")
models <- filenames(folder)
ks <- lapply(models,function(model) {
  load(paste(folder,model,".rdata",sep=""))
  K <- sapply(fit$samples,function(s) max(s$z))
  data.frame(model=model,mean.K=mean(K),sd.K=sd(K))
})
names(ks) <- models

for (m in c("uniform","marg","online")) {
  ks[[m]] <- data.frame(model=m,mean.K=NA,sd.K=NA)
}

df$mean.k <- sapply(as.character(df$L1),function(m) ks[[m]]$mean.K)
df$sd.k <- sapply(as.character(df$L1),function(m) ks[[m]]$sd.K)

  save(df,file=paste("results/",dataset,"/final/results.rdata",sep=""))

## Recall plots
load(paste("data/",dataset,".rdata",sep=""))
folder <- paste("results/",dataset,"/ranks/",sep="")
rks <- lapply(dir(folder,full.names=TRUE),function(f) {
  load(f)
  return(list(train=rk.train,test=rk.test))
})

models <- filenames(folder)
names(rks) <- models

ds <- lapply(rks,function(r) {
  d <- list(train=list(recall.30  = recall(r$train,top=1:30),
                       recall.200 = recall(r$train,top=1:100),
                       recall.all = recall(r$train,top=seq(1,N^2,length.out=100))),
            test =list(recall.30  = recall(r$test,top=1:30),
                       recall.200 = recall(r$test,top=1:100),
                       recall.all = recall(r$test,top=seq(1,N^2,length.out=100))))
  return(melt(d,id.vars="k",measure.vars="recall"))
})
for (i in 1:length(ds)) ds[[i]]$model <- models[i]
ds <- do.call(rbind,ds)
rownames(ds) <- c()
save(ds,file=paste("results/",dataset,"/final/recall.rdata",sep=""))

load(paste("results/",dataset,"/final/recall.rdata",sep=""))

# Recall at k
q3 <- qplot(k,value,data=subset(ds,L2=="recall.30"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q3)

q4 <- qplot(k,value,data=subset(ds,L2=="recall.200"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q4)
#ggsave(paste("figs/",opts$dataset,"/recall.200.pdf",sep=""),width=10,height=8)

q5 <- qplot(k,value,data=subset(ds,L2=="recall.all"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q5)
