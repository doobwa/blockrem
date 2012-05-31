datasets <- c("synthetic-1", "eckmann-small", "classroom-16", "classroom-17", "classroom-27", #"classroom-29", "classroom-31",
              "realitymining-small", "twitter-small", "enron-small")#, "irvine")
modelatts <- function(model) {
  atts <- strsplit(model,"\\.")[[1]]
  g <- function(x,y) as.numeric(strsplit(x,y)[[1]][2])
  xs <- c("kinit","kmax","sm","nb","pshift","deg","trans","collapse","xsigalpha","xsigbeta")
  a <- lapply(1:length(xs),function(i) g(atts[i],xs[i]))
  names(a) <- xs
  return(a)
}
modelnames <- function(folder) {
  as.vector(sapply(dir(folder),function(x) strsplit(x,".rdata")[[1]][1]))
}

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

ktb <- dcast(ks,type + dataset + xsigalpha + xsigbeta ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
ktb
itb <- dcast(iters,type + dataset + xsigalpha + xsigbeta ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
itb

tmp <- dcast(x,xsigalpha +xsigbeta+type+ dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp$xsigalpha <- tmp$xsigalpha/1000
tmp$xsigbeta  <- tmp$xsigbeta/1000
colnames(tmp)[1:2] <- c("alpha","beta")

tmp[,-c(1:4)] <- round(tmp[,-c(1:4)],3)
chosen <- c(1:4,8,12,15,16,23)
subset(tmp,alpha %in% c(1,5))[,chosen]
z=subset(tmp,type=="test")[,chosen]
z=subset(tmp,type=="test" & xsigalpha==1)[,c(1:4,8,12,15,16,23)]

write.table(z,file="results-5.dat",sep="\t")
write(print(z)
