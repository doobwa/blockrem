datasets <- c("synthetic-1", "eckmann-small", "classroom-16", "classroom-17", "classroom-27", "classroom-29", "classroom-31", "realitymining-small", "twitter-small", "enron-small")#, "irvine")
modelatts <- function(model) {
  atts <- strsplit(model,"\\.")[[1]]
  g <- function(x,y) as.numeric(strsplit(x,y)[[1]][2])
  xs <- c("kinit","kmax","sm","nb","pshift","deg","trans","collapse")
  a <- lapply(1:length(xs),function(i) g(atts[i],xs[i]))
  names(a) <- xs
  return(a)
}
modelnames <- function(folder) {
  as.vector(sapply(dir(folder),function(x) strsplit(x,".rdata")[[1]][1]))
}


res <- lapply(datasets,function(dataset) {
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
    list(train = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, train.ix - 1)[train.ix],
         test  = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, test.ix  - 1)[test.ix])
  })
  names(llks) <- models
  return(llks)
})
names(res) <- datasets

bad <- which(sapply(res,length) == 0)
if (length(bad) > 0) r <- res[-bad]
r <- melt(r)
colnames(r)[2:4] <- c("type","model","dataset")

ma <- lapply(unique(r$model),function(m) data.frame(model=m,modelatts(m)))
ma <- do.call(rbind,ma)
x <- merge(r,ma,by="model")

library(reshape2)
tmp <- dcast(x,type + dataset ~ kinit + kmax + pshift + deg, fun.aggregate=mean)
tmp
