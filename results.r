# Go through available pred objects.
source("utils.r")
datasets <- c("synthetic-1","classroom-16", "classroom-17", "classroom-27","eckmann-small","realitymining-small","enron-small","twitter-small")

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
allpred <- res

save(allpred,file="results/allpred.rdata")

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
library(reshape2)
r <- melt(r)
colnames(r)[2:5] <- c("type","metric","model","dataset")
ma <- lapply(unique(r$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(r,ma,by="model")
tmp <- dcast(x,xsigalpha +xsigbeta+metric+type+ dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[1:2] <- c("alpha","beta")
subset(tmp,metric=="mllk" & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]


a <- subset(tmp,metric=="llk" & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]
b <- subset(x,model %in% c("online","marg","uniform"))
b <- dcast(b,metric + type + dataset ~ model)
b <- b[1:length(datasets),c(3,6,4,5)]
a <- a[,-1]
final <- cbind(b,a)[c(7,3,4,5,6,8),]
colnames(final)[1:9] <- c("Dataset","\\texttt{unif}","\\texttt{marg}","\\texttt{online}","\\texttt{BM}","$K=1$","$K=2$","$K=3$","$K=10$")
final[,1] <- c("Synthetic","Classroom","University Email","Enron Email","Mobile Phone Calls","Twitter Dir. Messages")
library(xtable)
xr <- xtable(final,caption="Comparing mean loglikelihood for each event across methods for each dataset.  Larger values are better.  See text for details.",label="tab:results",digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-llk4.tex",sep=""),NA.string="",table.placement="t",size="footnotesize",sanitize.text.function=identity)


r <- allpred
ix <- which(sapply(r,length) == 0)
if (length(ix) > 0) r <- r[-ix]
rks <- list()
for (i in 1:length(r)) {
  rks[[i]] <- list()
  for (j in 1:length(r[[i]])) {
    tr <- recall(r[[i]][[j]]$rks$train)$recall#[c(5,20)]
    te <- recall(r[[i]][[j]]$rks$test)$recall#[c(5,20)]
    rks[[i]][[j]] <- list(train=list("5"=tr[1],"20"=tr[2]),
                          test =list("5"=te[1],"20"=te[2]))
  }
  names(rks[[i]]) <- names(r[[i]])
}
names(rks) <- names(r)

# 5 and 20 in same table
r <- melt(rks)
colnames(r)[2:5] <- c("metric","type","model","dataset")
ma <- lapply(unique(r$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(r,ma,by="model")
tmp <- dcast(x,dataset+xsigalpha + xsigbeta+type+ metric ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[2:3] <- c("alpha","beta")
subset(tmp,metric%in%c("5","20") & alpha==5 & type=="test")[,c("dataset","metric","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]



r <- melt(rks)
colnames(r)[2:5] <- c("metric","type","model","dataset")
ma <- lapply(unique(r$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(r,ma,by="model")
tmp <- dcast(x,xsigalpha +xsigbeta+metric+type+ dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[1:2] <- c("alpha","beta")
subset(tmp,metric=="5" & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]
subset(tmp,metric=="20" & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]

cutoff <- "20"
a <- subset(tmp,metric==cutoff & alpha==5 & type=="test")[,c("dataset","2_10_0_0","2_1_1_1","2_2_1_1","2_3_1_1","2_10_1_1")]
b <- subset(x,model %in% c("online","marg","uniform"))
b <- dcast(b,metric + type + dataset ~ model)
b <- subset(b,metric==cutoff & type=="test")
b <- b[,c(3,6,4,5)]
a <- a[,-1]
final <- cbind(b,a)[c(7,3,4,5,6),]
colnames(final)[1:9] <- c("Dataset","\\texttt{unif}","\\texttt{marg}","\\texttt{online}","\\texttt{BM}","$K=1$","$K=2$","$K=3$","$K=10$")
final[,1] <- c("Synthetic","Classroom","University Email","Enron Email","Mobile Phone Calls")#,"Twitter Direct Messages")
library(xtable)
xr <- xtable(final,caption=paste("Comparing recall at cutoff",cutoff,"across methods for each dataset.  Larger values are better.  See text for details."),label=paste("tab:recall",cutoff,sep=""),digits=3)
print(xr,include.rownames=FALSE,file=paste("figs/results-recall-",cutoff,".tex",sep=""),NA.string="",table.placement="t",size="footnotesize",sanitize.text.function=identity)


###################
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
