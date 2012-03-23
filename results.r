
res <- lapply(c("synthetic","eckmann-small"),function(x) {
  load(paste("results/",x,"/final/results.rdata",sep=""))
  return(df)
})
res <- do.call(rbind,res)
res$dataset <- factor(res$dataset,c("synthetic","eckmann-small"))

r <- cast(res,dataset + likelihood + type ~ L1,add.missing=TRUE)

# Clean up
colnames(r)[7:9] <- paste("K=",1:3,sep="")
r$dataset    <- as.character(r$dataset)
r$likelihood <- as.character(r$likelihood)
r$dataset[-seq(1,nrow(r),by=4)] <- ""
r$likelihood[-seq(1,nrow(r),by=2)] <- ""
#for (i in 4:ncol(r)) r[,i] <- round(r[,i],3)

 r <- xtable(r,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details. [TODO: More datasets as well as results for the DP version.  Could also remove training set results.]",label="tab:results",digits=3)
print(r,include.rownames=FALSE,file=paste("figs/results.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")
