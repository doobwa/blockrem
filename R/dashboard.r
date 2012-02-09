#!Rscript

dataset <- "eckmann-small"

library(ggplot2)

#' Perform analysis for a given dataset given the the current results of model fits.  
#' @dataset name of dataset.  Searches /results/[dataset]/ for results, and saves figures to /figs/[dataset]/
#' Each fit object should have:
#' @z latent class assignments
#' @beta values of parameters in the last iteration
#' @llks log posterior at each MCMC iteration
#' @param array of parameter values at each MCMC iteration
#' @zs list of latent class assignment vectors at each MCMC iteration
#' The synthetic dataset should also have true.lpost variable available.
results.dir <- paste("results/",dataset,sep="")

# Pull data from each saved file and grab the name of the fit
fits <- lapply(dir(results.dir,full.names=TRUE),function(f) {
  load(f)
  return(res)
})
fnames <- unlist(strsplit(dir(results.dir),".rdata"))

# Get likelihood for each fit.
niter <- 2000
llks <- lapply(fits,function(f) f$llks)#cbind(llk=f$llks,iter=1:f$niter))
names(llks) <- fnames
llks <- melt(llks)
llks$iter <- 1:niter
q <- qplot(iter,value,data=subset(llks,iter>10),geom="line",colour=factor(L1)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
if (dataset == "synthetic") {
  q + geom_abline(intercept=true.lpost,slope=0)
} else {
  q
}
ggsave(paste("figs/",dataset,"/lpost.pdf",sep=""),height=4,width=5)


# Plot progress of z at each MCMC iteration
zs <- lapply(fits,function(f) do.call(rbind,f$zs))
names(zs) <- fnames
zs <- melt(zs)
qplot(X1,X2,data=zs,geom="tile",fill=factor(value)) + facet_wrap(~L1)
ggsave(paste("figs/",dataset,"/zs.pdf",sep=""),height=4,width=5)


# Compute dyad counts for each pshift
source("R/utils.r")
df <- dyad.ps(train,N)
df <- melt(df)
df$i <- z[df$X1]
df$j <- z[df$X2]
qplot(X3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="shift type",y="count for a given dyad") + opts(axis.text.x=theme_text(angle=90))
ggsave(paste("figs/",dataset,"/counts.pdf",sep=""),width=6,height=4)


# Prediction experiment on test data: precision
table(test[,2],test[,3])
lrms <- list(unif = array(1,c(M,N,N)),
             true = brem.lrm(test$A,N,z,beta),
             base = sbm.lrm(test$A,N,fit0$z,fit0$beta),
             diag = brem.lrm(test$A,N,fit1$z,fit1$beta),
             full = brem.lrm(test$A,N,fit2$z,fit2$beta),
             sing = brem.lrm(test$A,N,fit3$z,fit3$beta))
ps <- lapply(lrms,function(lrm) {
  recall(ranks(test$A,-lrm,ties.method="random"),top=1:100)
})
res <- melt(ps,id.vars=c("k"),measure.vars="recall")
qplot(k,value,data=res,geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
qplot(k,value,data=subset(res,k<50),geom="line",colour=factor(L1),group=factor(L1))+theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
ggsave("figs/syn/test-recall.pdf",width=5,height=4)


# Compute out of sample log posterior
lposts <- list(true = brem.lpost(test$A,N,K,z,beta),
               base = sbm.lpost(test$A,N,K,fit0$z,fit0$beta),
               diag = brem.lpost(test$A,N,K,fit1$z,fit1$beta),
               full = brem.lpost(test$A,N,K,fit2$z,fit2$beta),
               sing = brem.lpost(test$A,N,K,fit3$z,fit3$beta))
unlist(lposts)
save(lposts,file="data/syn/lpost.rdata")
