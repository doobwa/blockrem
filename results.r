
datasets <- c("synthetic-1","eckmann-small", "realitymining-small", "classroom-16", "classroom-17", "classroom-27", "classroom-29", "classroom-31", "enron-small")
res <- lapply(datasets,function(x) {
  f <-paste("results/",x,"/final/results.rdata",sep="")
  if (file.exists(f)) load(f)
  return(df)
})
names(res) <- datasets
lapply(res,dim)

res <- do.call(rbind,res)

res$dataset <- factor(res$dataset,datasets)
res <- res[,-c(6,7)]
library(reshape2)
res <- res[,c(1:3,5,4)]

r <- subset(res,type=="test" & L1 %in% c("kinit10.sm0.nb1.deg0","marg","online","uniform"))
r <- dcast(r,dataset + metric ~ L1,fun.aggregate=sum)

r <- subset(r, !dataset %in% c("enron-small","realitymining-small"))
colnames(r)[3] <- c("brem")
r$dataset    <- as.character(r$dataset)
r$metric <- as.character(r$metric)
r <- r[,c(1,2,6,4,5,3)]
r$dataset[-seq(1,nrow(r),by=2)] <- ""

a <- subset(r, metric %in% c("rem","mult"))
#r$likelihood[-seq(1,nrow(r),by=2)] <- ""
library(xtable)
xr <- xtable(a,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-llk.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

b <- subset(r, metric %in% c("rank5","rank20"))
b$metric <- c("@  5", "@ 20")
xr <- xtable(r,caption="Predictive accuracy on a recall task.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-recall.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")


# Clean up
colnames(r)[7:9] <- paste("K=",1:3,sep="")
r$dataset    <- as.character(r$dataset)
r$likelihood <- as.character(r$likelihood)
r$dataset[-seq(1,nrow(r),by=4)] <- ""
r$likelihood[-seq(1,nrow(r),by=2)] <- ""
#for (i in 4:ncol(r)) r[,i] <- round(r[,i],3)

 r <- xtable(r,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details. [TODO: More datasets as well as results for the DP version.  Could also remove training set results.]",label="tab:results",digits=3)
print(r,include.rownames=FALSE,file=paste("figs/results.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")
