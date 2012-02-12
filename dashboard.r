#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata and results at /results/[dataset]/."),
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


cat("Plotting log posterior during MCMC.\n")
niter <- 500
llks <- lapply(fits,function(f) f$llks)#cbind(llk=f$llks,iter=1:f$niter))
names(llks) <- fnames
llks <- melt(llks)
llks$iter <- 1:niter

if (opts$dataset == "synthetic") {
  q1 <- qplot(iter,value,data=subset(llks,iter>10),geom="line",colour=factor(L1)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw() + geom_abline(intercept=true.lpost,slope=0)
} else {
  q1 <- qplot(iter,value,data=subset(llks,iter>10),geom="line",colour=factor(L1)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
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


cat("Recall experiment on test data.\n")
lrms <- lapply(fits,function(f) {
  brem.lrm(test,N,f$z,f$beta)
})
names(lrms) <- fnames
lrms$unif   <- array(1,c(nrow(test),N,N))
lrms$online <- ratemat.online(test,N)
lrms$marg   <- ratemat.from.marginals(train,test,N)
if (opts$dataset == "synthetic") {
  load("data/synthetic.rdata")
  lrms$true <- brem.lrm(test,N,z,beta)
}

# TODO: add SBM baseline
# TODO: add 

ps <- lapply(lrms,function(lrm) {
  recall(ranks(test,-lrm,ties.method="random"),top=1:100)
})
res <- melt(ps,id.vars=c("k"),measure.vars="recall")
q4 <- qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")

chosen <- seq(1,N^2,length.out=100)
ps <- lapply(lrms,function(lrm) {
  recall(ranks(test,-lrm,ties.method="random"),top=chosen)
})
res <- melt(ps,id.vars=c("k"),measure.vars="recall")
q5 <- qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")

cat("Plotting parameter estimates.\n")
betas <- lapply(fits,function(f) f$beta)
names(betas) <- fnames
b <- melt(betas)
q6 <- qplot(X1,value,data=b,geom="point",colour=factor(L1)) + facet_grid(X2~X3) + theme_bw() + labs(colour="model")


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


cat("Creating dashboard.\n")
library(gridExtra)
pdf(paste("figs/",opts$dataset,"/dashboard.pdf",sep=""),width=20,height=10)
blankPanel <- grid.rect(gp=gpar(col="white"))
grid.arrange(q1, q2, q3, q4, q5, q6, ncol=3)
dev.off()

#ggsave(paste("figs/",opts$dataset,"/dashboard.pdf",sep=""),width=20,height=10)

cat("Complete.\n")

if (opts$save.figs) {
  print(q1)
  ggsave(paste("figs/",opts$dataset,"/lpost.pdf",sep=""),height=4,width=5)
  print(q2)
  ggsave(paste("figs/",opts$dataset,"/zs.pdf",sep=""),height=4,width=5)
  print(q3)
  ggsave(paste("figs/",opts$dataset,"/counts.pdf",sep=""),width=6,height=4)
  print(q4)
  ggsave(paste("figs/",opts$dataset,"/test-recall.pdf",sep=""),width=5,height=4)
}

