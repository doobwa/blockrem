# Go through available pred objects.

datasets <- c("synthetic-1","eckmann-small","classroom-16", "classroom-17", "classroom-27","realitymining-small","enron-small","twitter")

res <- lapply(datasets,function(dataset) {
  folder <- paste("results/",dataset,"/preds/",sep="")
  models <- modelnames(folder)
  r <- lapply(models,function(model) {
    load(paste(folder,model,".rdata",sep=""))
    return(pred)
  })
  names(r) <- models
  return(r)
})
names(res) <- datasets

r <- res
ix <- which(sapply(r,length) == 0)
if (length(ix) > 0) r <- r[-ix]
for (i in 1:length(r)) {
  for (j in 1:length(r[[i]])) {
    for (k in 1:length(r[[i]][[j]])) {
      r[[i]][[j]][[k]] <- list(train=mean(r[[i]][[j]][[k]]$train),
                               test =mean(r[[i]][[j]][[k]]$test))
    } 
  }
}
# recall(r$test,top=1:30)
r <- melt(r)
colnames(r)[2:5] <- c("type","metric","model","dataset")
ma <- lapply(unique(r$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(r,ma,by="model")
tmp <- dcast(x,xsigalpha +xsigbeta+metric+type+ dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[1:2] <- c("alpha","beta")

subset(tmp,metric=="llk" & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]

#datasets <- c("synthetic-1","eckmann-small","realitymining-small",
 datasets <- c("synthetic-1","eckmann-small","classroom-16", "classroom-17", "classroom-27","realitymining-small","enron-small","twitter")
res <- lapply(datasets,function(x) {
  f <- paste("results/",x,"/final/results.rdata",sep="")
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

#res$dataset <- factor(res$dataset,datasets)
#res <- res[,-c(6,7)]
library(reshape2)
res <- res[,c(1:3,5,4)]

r <- res
chosen.models <- c("kinit1.sm0.nb0.pshift1.deg1.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift0.deg0.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift1.deg0.trans1.collapse0",
                   "kinit10.sm0.nb0.pshift1.deg1.trans1.collapse0",
                   "kinit1.sm0.nb0.pshift1.deg1.trans1.collapse1",
                   "kinit10.sm0.nb0.pshift1.deg0.trans1.collapse1",
                  "kinit10.sm0.nb0.pshift0.deg0.trans1.collapse1",
                  "kinit10.sm0.nb0.pshift1.deg1.trans1.collapse1")
r1 <- subset(res,!L1 %in% c("marg","online","uniform"))#chosen.models)
r2 <- subset(res,L1 %in% c("marg","online","uniform"))
r <- dcast(rbind(r1,r2),dataset + metric + type ~ L1,fun.aggregate=sum)
                   
colnames(r)[4:15] <- c("100","101","110","111","300","301","310","311","500","501","510","511")#c("BREM.0","BREM.0.c","BREM.01","BREM.01.c","BREM.012","BREM.012.c","K=1")#,"K=1.c")

#r <- subset(r, !dataset %in% c("realitymining-small"))
r$dataset    <- as.character(r$dataset)
r$metric <- as.character(r$metric)
r$type <- as.character(r$type)
#r <- r[c("type","dataset","metric","uniform","marg","online","K=1","BREM.0","BREM.01","BREM.012","BREM.0.c","BREM.01.c","BREM.012.c")]
#r$dataset[-seq(1,nrow(r),by=2)] <- ""

# Add truth for synthetic data
#r$truth <- NA
#r$truth[1:4] <- subset(res,L1 == "truth" & type=="test")$value

subset(r, metric %in% c("rem","mult") & type=="train")[,-c(4,5,6,19,18,16)]
subset(r, metric %in% c("rem") & type=="train")[,-c(4,5,6,19,18,16)]
subset(r, metric %in% c("rem") & type=="test")[,-c(19,18,16)]
subset(r, metric %in% c("rank5","rank20") & type=="train")
subset(r, metric %in% c("rank5","rank20") & type=="test")

a <- subset(r, metric %in% c("rem","mult") & type=="test")
b <- subset(r, metric %in% c("rank5","rank20") & type=="test")
b$metric <- c("@  5", "@ 20")

#r$likelihood[-seq(1,nrow(r),by=2)] <- ""
library(xtable)
xr <- xtable(a,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-llk.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

xr <- xtable(b,caption="Predictive accuracy on a recall task.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-recall.tex",sep=""),NA.string="",table.placement="t",size="footnotesize")

