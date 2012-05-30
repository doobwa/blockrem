
datasets <- c("synthetic-1","eckmann-small","realitymining-small", "classroom-16", "classroom-17", "classroom-27", "enron-small")
res <- lapply(datasets,function(x) {
  f <-paste("results/",x,"/final/results.rdata",sep="")
  if (file.exists(f)) {
    load(f)
    return(df[,1:5])
  } else {
    return(data.frame())
  }
})
names(res) <- datasets

bad <- which(sapply(res,dim)[1,] == 0)
if (length(bad) > 0) res <- res[-bad]
t(sapply(res,dim))

res <- do.call(rbind,res)

res$dataset <- factor(res$dataset,datasets)
#res <- res[,-c(6,7)]
library(reshape2)
res <- res[,c(1:3,5,4)]

r <- subset(res,type=="test")
#chosen.models <- c("kinit1.sm0.nb0.pshift1.deg1.trans1.collapse1",
#                   "kinit10.sm0.nb0.pshift0.deg0.trans1.collapse1",
#                   "kinit10.sm0.nb0.pshift1.deg1.trans1.collapse1")
chosen.models <- c("kinit1.sm0.nb0.pshift1.deg1.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift0.deg0.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift1.deg0.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift1.deg1.trans1.collapse0")
r1 <- subset(res,L1 %in% chosen.models & type=="test")
r2 <- subset(res,L1 %in% c("marg","online","uniform") & type=="test")
r <- dcast(rbind(r1,r2),dataset + metric ~ L1,fun.aggregate=sum)
colnames(r)[3:6] <- c("BREM.0","BREM.01","BREM.012","K=1")

#r <- subset(r, !dataset %in% c("realitymining-small"))
r$dataset    <- as.character(r$dataset)
r$metric <- as.character(r$metric)
r <- r[c("dataset","metric","uniform","marg","online","K=1","BREM.0","BREM.01","BREM.012")]
r$dataset[-seq(1,nrow(r),by=2)] <- ""

# Add truth for synthetic data
r$truth <- NA
r$truth[1:4] <- subset(res,L1 == "truth" & type=="test")$value

a <- subset(r, metric %in% c("rem","mult"))
b <- subset(r, metric %in% c("rank5","rank20"))
b$metric <- c("@  5", "@ 20")

#r$likelihood[-seq(1,nrow(r),by=2)] <- ""
library(xtable)
xr <- xtable(a,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-llk.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

xr <- xtable(b,caption="Predictive accuracy on a recall task.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-recall.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

