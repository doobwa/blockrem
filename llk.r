datasets <- c("synthetic-1", "eckmann-small", "classroom-16", "classroom-17", "classroom-27", #"classroom-29", "classroom-31",
              "realitymining-small", "twitter-small", "enron-small")#, "irvine")
source("utils.r")

library(brem)
library(multicore)
library(reshape2)
options(cores=8)
res <- mclapply(datasets,function(dataset) {
  cat(dataset,"\n")
  load(paste("data/",dataset,".rdata",sep=""))
  M <- nrow(A)
  P <- 13
  ego <- 1
  train.ix <- 1:nrow(train)
  test.ix  <- (1:nrow(A))[-train.ix]
  s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,M,P,ego)
  s$precompute()
  s$transform()
  folder <- paste("results/",dataset,"/fits/",sep="")
  models <- modelnames(folder)
  llks <- lapply(models,function(model) {
    load(paste(folder,model,".rdata",sep=""))
    K <- dim(fit$params$beta)[2]
    list(K = sapply(fit$samples,function(s) max(s$z)),
         iter = fit$iter,
         train = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, train.ix - 1)[train.ix],
         test  = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, test.ix  - 1)[test.ix])
  })
  names(llks) <- models
  return(llks)
})
names(res) <- datasets

r <- res
bad <- which(sapply(res,length) == 0)
if (length(bad) > 0) r <- res[-bad]
r <- melt(r)
colnames(r)[2:4] <- c("type","model","dataset")
llks <- subset(r,! type %in% c("K","iter"))
ks <-  subset(r,type == "K")
iters <- subset(r,type=="iter")

ma <- lapply(unique(llks$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(llks,ma,by="model")
ks <- merge(ks,ma,by="model")
iters <- merge(iters,ma,by="model")

save(llks,ks,iters,x,ma,file="results/llk.rdata")
load("results/llk.rdata")

ktb <- dcast(ks,dataset + xsigalpha + xsigbeta ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
ktb
itb <- dcast(iters,dataset + xsigalpha + xsigbeta ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
itb

tmp <- dcast(x,xsigalpha +xsigbeta+type+ dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[1:2] <- c("alpha","beta")

tmp[,-c(1:4)] <- round(tmp[,-c(1:4)],3)
chosen <- c(1:4,8,12,15,16,23)
subset(tmp,alpha %in% c(1,5))[,chosen]

res <- subset(tmp,alpha %in% c(5) & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]
colnames(res) <- c("Dataset","SBM","K=1","K=2","K=3","K=10")

library(xtable)
xr <- xtable(res,caption="Comparing mean loglikelihood under different methods for each event in a given test set.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-llk2.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

