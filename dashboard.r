#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
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
fits <- lapply(dir(results.dir,full.names=TRUE),function(f) {
  load(f)
  return(res)
})
fnames <- unlist(strsplit(dir(results.dir),".rdata"))
names(fits) <- fnames

cat("Plotting log posterior during MCMC.\n")
llks <- lapply(1:length(fits),function(i) {
  data.frame(model=fnames[i],llk=fits[[i]]$llks,iter=1:fits[[i]]$niter)
})
llks <- do.call(rbind,llks)
llks <- subset(llks,iter>10)

if (opts$dataset == "synthetic") {
  load("data/synthetic.rdata")
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model))+ geom_abline(intercept=true.lpost,slope=0) + labs(x="iteration",y="log posterior",colour="model") + theme_bw() + limits(c(min(llks$llk),0),"y")
} else {
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
}

cat("Plotting progress of z at each MCMC iteration.\n")
zs <- lapply(fits,function(f) do.call(rbind,f$zs))
names(zs) <- fnames
zs <- melt(zs)
q2 <- qplot(X1,X2,data=zs,geom="tile",fill=factor(value)) + facet_wrap(~L1) + labs(x="iteration",y="node")


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
if (opts$dataset == "synthetic") {
  load("data/synthetic.rdata")
  fits$truth <- list(z=z,beta=beta)
}

# Temp
P <- 13
K <- 1
tmp <- array(0,c(P,K,K))
tmp[12,,] <- 1
fits$counts.only <- list(z=rep(1,N),beta=tmp)

cat("precomputing\n")
strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
strain$precompute()
stest <- new(RemStat,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test),P)
stest$precompute()

options(cores=5)
fnames <- c(names(fits),c("uniform","online","margins"))
rks <- list()
for (i in 1:length(fnames)) {
  get.lrm <- switch(fnames[i],
                    "uniform" = function(x) array(1,c(nrow(x),N,N)),
                    "online"  = function(x) ratemat.online(x,N),
                    "margins" = function(x) ratemat.from.marginals(train,x,N),
                    function(x) brem.lrm.fast(nrow(x),s,fits[[i]]$z,fits[[i]]$beta))
                    #function(x) brem.lrm(x,N,fits[[i]]$z,fits[[i]]$beta)),
  
  cat("train:",fnames[i])
  s <- strain
  lrm <- get.lrm(train)
  cat("ranking\n")
  rk1 <- ranks(train,-lrm,ties.method="random")
  cat("test:",fnames[i])
  s <- stest
  lrm <- get.lrm(test)
  cat("ranking\n")
  rk2 <- ranks(test,-lrm,ties.method="random")
  rm(lrm)
  gc()
  rks[[i]] <- list(model=fnames[i],train=rk1,test=rk2)
}
for (i in 1:length(fnames)) {
  rank.k[[i]]$model <- fnames[i]
}
d <- list(train=list(rank.100 = recall(rk1,top=1:100),
                     rank.all = recall(rk1,top=seq(1,N^2,length.out=100))),
          test =list(rank.100 = recall(rk2,top=1:100),
                     rank.all = recall(rk2,top=seq(1,N^2,length.out=100))))
rank.k[[i]] <- melt(d,id.vars="k",measure.vars="recall")
rank.dyad <- data.frame(model=fnames[i],train=mea)

rank.k <- do.call(rbind,rank.k)
rownames(rks) <- c()

# Recall at k
q4 <- qplot(k,value,data=rks,geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2) +theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q4)
ggsave(paste("figs/",opts$dataset,"/recall.pdf",sep=""),width=10,height=8)

q5 <- qplot

cat("Plotting parameter estimates.\n")
betas <- lapply(fits,function(f) f$beta)
names(betas) <- fnames
b <- melt(betas)
q6 <- qplot(X1,value,data=b,geom="point",colour=factor(L1)) + facet_grid(X2~X3) + theme_bw() #+ labs(colour="model")


# # Compute out of sample log posterior
lposts <- lapply(fits,function(f) {
  brem.lpost(test,N,opts$K,f$z,f$beta,priors=list(beta=list(mu=0,sigma=1)))
})
# lposts <- list(true = brem.lpost(test$A,N,K,z,beta),
#                base = sbm.lpost(test$A,N,K,fit0$z,fit0$beta),
#                diag = brem.lpost(test$A,N,K,fit1$z,fit1$beta),
#                full = brem.lpost(test$A,N,K,fit2$z,fit2$beta),
#                sing = brem.lpost(test$A,N,K,fit3$z,fit3$beta))
# unlist(lposts)


# trace plots
library(coda)
load(paste("results/",opts$dataset,"/full.1.rdata",sep=""))
r <- melt(res$param)
q7 <- qplot(X1,value,data=r, colour=factor(X4),geom="line") + labs(colour="parameters for\n 1x1 block",x="iteration") + theme_bw() + facet_grid(X2~X3)

cat("Creating dashboard.\n")
library(gridExtra)
pdf(paste("figs/",opts$dataset,"/dashboard.pdf",sep=""),width=25,height=10)
blankPanel <- grid.rect(gp=gpar(col="white"))
grid.arrange(q1, q2, q3, q6,q7,q4,q6, ncol=5)
dev.off()

save.image(paste("figs/",opts$dataset,"/dashboard.rdata",sep=""))

cat("Complete.\n")

if (opts$save.figs) {
  print(q1)
  ggsave(paste("figs/",opts$dataset,"/lpost.pdf",sep=""),height=4,width=5)
  print(q2)
  ggsave(paste("figs/",opts$dataset,"/zs.pdf",sep=""),height=4,width=5)
  print(q3)
  ggsave(paste("figs/",opts$dataset,"/counts.pdf",sep=""),width=6,height=4)
  print(q4)
  ggsave(paste("figs/",opts$dataset,"/recall.pdf",sep=""),width=5,height=4)
#   print(q5)
#   ggsave(paste("figs/",opts$dataset,"/test-recall.pdf",sep=""),width=5,height=4)
#   print(q4a)
#   ggsave(paste("figs/",opts$dataset,"/train-recall-zoom.pdf",sep=""),width=5,height=4)
#   print(q5a)
#   ggsave(paste("figs/",opts$dataset,"/train-recall.pdf",sep=""),width=5,height=4)
  print(q6)
  ggsave(paste("figs/",opts$dataset,"/parameters.pdf",sep=""),width=5,height=4)
  print(q7)
  ggsave(paste("figs/",opts$dataset,"/traceplots.pdf",sep=""),width=5,height=4)
}

