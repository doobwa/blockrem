#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-p", "--predictions"), default=FALSE,
              help="compute prediction experiments?"),
  make_option(c("-s", "--save.figs"),default=FALSE,
              help="Save each individual figure to figs/[dataset]/ [default %default]."))
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#' Perform analysis for a given dataset given the the current results of model fits.  
#' @dataset name of dataset.  Searches /results/[dataset]/ for results, and saves figures to /figs/[dataset]/
#' 
#' Each "res" object should have:
#' @z latent class assignments
#' @beta values of parameters in the last iteration
#' @llks log posterior at each MCMC iteration
#' @param array of parameter values at each MCMC iteration
#' @zs list of latent class assignment vectors at each MCMC iteration
#' The synthetic dataset should also have true.lpost variable available as well as the true beta parameters.
#' 
options(verbose=FALSE)
library(brem)
library(ggplot2)

# Pull data from each saved file and grab the name of the fit
results.dir <- paste("results/",opts$dataset,sep="")
load(paste("data/",opts$dataset,".rdata",sep=""))
fs <- list.files(results.dir,full.names=TRUE)
fs <- fs[-grep("ranks",fs)]
fs <- fs[-grep("llks",fs)]
fs <- fs[-grep("counts",fs)]
fnames <- unlist(strsplit(fs,".rdata"))
fitnames <- sapply(fnames,function(f) strsplit(f,paste(results.dir,"/",sep=""))[[1]][2])
fits <- lapply(fs,function(f) {
  load(f)
  return(res)
})
names(fits) <- fitnames

cat("Plotting log posterior during MCMC.\n")
llks <- lapply(1:length(fits),function(i) {
  data.frame(model=names(fits)[i],llk=fits[[i]]$llks,iter=1:fits[[i]]$niter)
})
llks <- do.call(rbind,llks)

llks <- subset(llks,iter<100)
if (opts$dataset == "synthetic") {
  load("data/synthetic.rdata")
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model))+ geom_abline(intercept=true.lpost,slope=0) + labs(x="iteration",y="log posterior",colour="model") + theme_bw() + limits(c(min(llks$llk),0),"y")
} else {
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
}

library(coda)
load(paste("results/",opts$dataset,"/full.2.rdata",sep=""))
r <- melt(res$param)
r <- subset(r,X1 < 10)
q1a <- qplot(X1,value,data=r, colour=factor(X2),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw() + facet_grid(X3~X4,scales="free")


#qplot(iter,log(-llk + 1),data=subset(llks,iter < 20 & iter > 15),geom="line",colour=factor(model)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()

cat("Plotting progress of z at each MCMC iteration.\n")
zs <- lapply(fits,function(f) do.call(rbind,f$zs))
names(zs) <- names(fits)
zs <- melt(zs)
q2 <- qplot(X1,X2,data=zs,geom="tile",fill=factor(value)) + facet_wrap(~L1) + labs(x="iteration",y="node")

cat("Compute distribution of activity using results from the full model.\n")
tb <- table(factor(c(train[,2],train[,3]),1:N))
fx <- grep("full.2",names(fits))
z <- fits[[fx]]$z
df <- data.frame(group=z,count=tb)
q3a <- qplot(log(count.Freq),data=df,geom="histogram")+facet_grid(group~.,scales="free")+labs(y="number of users",x="log(number of events)") + theme_bw()

cat("Compute dyad counts for each pshift using results from the full model.\n")
fx <- grep("full.2",fnames)
if (length(fx) > 0) {
  z <- fits[[fx]]$z
  df <- dyad.ps(train,N)
  df <- melt(df)
  df$i <- z[df$X1]
  df$j <- z[df$X2]
  q3 <- qplot(X3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="shift type",y="count for a given dyad") + opts(axis.text.x=theme_text(angle=90))
} else {
  q3 <- qplot(0,0,label="not available",geom="text")
}

cat("Recall experiment on training data.\n")

load(paste("data/",opts$dataset,".rdata",sep=""))
folder <- paste("results/",opts$dataset,"/ranks/",sep="")
rks <- lapply(dir(folder,full.names=TRUE),function(f) {
  load(f)
  return(list(train=rk.train,test=rk.test))
})
fnames <- unlist(strsplit(dir(folder),".rdata"))
names(rks) <- fnames

ds <- lapply(rks,function(r) {
  d <- list(train=list(recall.200 = recall(r$train,top=1:100),
                       recall.all = recall(r$train,top=seq(1,N^2,length.out=100))),
            test =list(recall.200 = recall(r$test,top=1:100),
                       recall.all = recall(r$test,top=seq(1,N^2,length.out=100))))
  return(melt(d,id.vars="k",measure.vars="recall"))
})
for (i in 1:length(ds)) ds[[i]]$model <- fnames[i]
ds <- do.call(rbind,ds)
rownames(ds) <- c()

# Recall at k
q4 <- qplot(k,value,data=subset(ds,L2=="recall.200"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q4)
#ggsave(paste("figs/",opts$dataset,"/recall.200.pdf",sep=""),width=10,height=8)

q5 <- qplot(k,value,data=subset(ds,L2=="recall.all"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q5)
#ggsave(paste("figs/",opts$dataset,"/recall.all.pdf",sep=""),width=10,height=8)

# # Experimental code for looking at the relation between the rank and the observed activity
# tbtrain <- table(factor(c(train[,2],train[,3]),1:N))
# tbtest <- table(factor(c(test[,2],test[,3]),1:N))
# rks.all <- lapply(1:length(rks),function(i) 
#   rbind(
#   data.frame(model=fnames[i],
#              iter=1:length(rks[[i]]$train),
#              rk=rks[[i]]$train,type="train",
#              tc=tbtrain[train[,2]] + tbtrain[train[,3]]),
#   data.frame(model=fnames[i],
#              iter=1:length(rks[[i]]$test),
#              rk=rks[[i]]$test,type="test",
#              tc=tbtest[test[,2]] + tbtest[test[,3]])
#   )
# )
# rks.all <- do.call(rbind,rks.all)
# qplot(iter,rk,data=subset(rks.all,type=="train"),colour=factor(model),alpha=.2)+facet_grid(model~.)+theme_bw()
# ggsave(paste("figs/",opts$dataset,"/rank.train.pdf",sep=""),width=20,height=12)
# qplot(iter,rk,data=subset(rks.all,type=="test"),colour=factor(model),alpha=.2)+facet_grid(model~.)+theme_bw()
# ggsave(paste("figs/",opts$dataset,"/rank.test.pdf",sep=""),width=20,height=12)
# qplot(tc,rk,data=rks.all,colour=factor(model),alpha=.1)+facet_grid(type ~ model)+theme_bw() + labs(x="sum of sender and recipient degree",y="rank of observed event")
# ggsave(paste("figs/",opts$dataset,"/rank.v.count.pdf",sep=""),width=20,height=12)

cat("Plotting parameter estimates.\n")
betas <- lapply(fits,function(f) f$beta)
names(betas) <- names(fits)
b <- melt(betas)
q6 <- qplot(X1,value,data=b,geom="point",colour=factor(L1)) + facet_grid(X2~X3) + theme_bw() #+ labs(colour="model")
#q6 <- qplot(X1,value,data=b,geom="point") + facet_grid(X2~X3) + theme_bw() #+ labs(colour="model")

# b <- melt(fits[['results/twitter-small/full.3']]$beta)
# b$block <- paste(b$X2,b$X3)
# qplot(factor(X1),value,data=b,geom="point",colour=factor(block))+theme_bw()+labs(x="parameter")

if (opts$predictions) {
  # # Compute out of sample log posterior
  load(paste("data/",opts$dataset,".rdata",sep=""))
  folder <- paste("results/",opts$dataset,"/llks/",sep="")
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
  names(llks.test) <- names(llks.train) <- names(mllks.test) <- names(mllks.train) <- strsplit(dir(folder),".rdata")
  
  llks.train <- melt(llks.train)
  llks.train$event <- 1:nrow(train)
  llks.test <- melt(llks.test)
  llks.test$event <- 1:nrow(test)
  mllks.train <- melt(mllks.train)
  mllks.train$event <- 1:nrow(train)
  mllks.test <- melt(mllks.test)
  mllks.test$event <- 1:nrow(test)
  
  # Examine log likelihood of observations
  qplot(event,value,data=llks.train,geom="point",colour=factor(L1))
  qplot(event,value,data=llks.test,geom="point",colour=factor(L1))
  qplot(event,value,data=mllks.train,geom="point",colour=factor(L1))
  qplot(event,value,data=mllks.test,geom="point",colour=factor(L1)) 
  
  qplot(event,value,data=subset(llks.test,L1 %in% c("full.2","online")),geom="point",colour=factor(L1))
  
  tmp <- subset(llks.test,L1 %in% c("full.2","online"))
  
  # window.size <- 200
  # llks.train <- subset(llks.train,event>window.size)
  # llks.test <- subset(llks.test,event>window.size)
  # mllks.train <- subset(mllks.train,event>window.size)
  # mllks.test <- subset(mllks.test,event>window.size)
  
  r <- cbind(ddply(llks.train,.(L1),summarise,llk=mean(value)),
             ddply(llks.test,.(L1),summarise,llk=mean(value)),
             ddply(mllks.train,.(L1),summarise,llk=mean(value)),
             ddply(mllks.test,.(L1),summarise,llk=mean(value)))
  r <- r[,c(1,2,4,6,8)]
  colnames(r) <- c("method","brem.train","brem.test","multin.train","multin.test")
  for (i in 2:5) r[,i] <- round(r[,i],2)
  r
  print(xtable(r),include.rownames=FALSE,file=paste("figs/",opts$dataset,".llk.tex",sep=""))
}
# 
# load("results/eckmann-small/llks/online.rdata")
# load("results/eckmann-small/counts.rdata")
# llk.test <- llk.test
# ctest <- ctest
# par(mfrow=c(3,1))
# plot(llk.test,xlab="event")
# plot(ctest,xlab="event",ylab="dyad count for observed event")
# plot(ctest,llk.test,xlab="dyad count for observed event")

# trace plots

cat("Creating dashboard.\n")
library(gridExtra)
pdf(paste("figs/",opts$dataset,"/dashboard.pdf",sep=""),width=25,height=10)
blankPanel <- grid.rect(gp=gpar(col="white"))
grid.arrange(q1, q1a, q2, q3, q3a, q4, q5, q6,ncol=4)
dev.off()

#save.image(paste("figs/",opts$dataset,"/dashboard.rdata",sep=""))

# cat("Complete.\n")
# 
# if (opts$save.figs) {
#   print(q1)
#   ggsave(paste("figs/",opts$dataset,"/lpost.pdf",sep=""),height=4,width=5)
#   print(q2)
#   ggsave(paste("figs/",opts$dataset,"/zs.pdf",sep=""),height=4,width=5)
#   print(q3)
#   ggsave(paste("figs/",opts$dataset,"/counts.pdf",sep=""),width=6,height=4)
#   print(q4)
#   ggsave(paste("figs/",opts$dataset,"/recall.pdf",sep=""),width=5,height=4)
# #   print(q5)
# #   ggsave(paste("figs/",opts$dataset,"/test-recall.pdf",sep=""),width=5,height=4)
# #   print(q4a)
# #   ggsave(paste("figs/",opts$dataset,"/train-recall-zoom.pdf",sep=""),width=5,height=4)
# #   print(q5a)
# #   ggsave(paste("figs/",opts$dataset,"/train-recall.pdf",sep=""),width=5,height=4)
#   print(q6)
#   ggsave(paste("figs/",opts$dataset,"/parameters.pdf",sep=""),width=5,height=4)
#   print(q7)
#   ggsave(paste("figs/",opts$dataset,"/traceplots.pdf",sep=""),width=5,height=4)
# }
# 
