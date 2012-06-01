
opts <- list(dataset="classroom-16",predictions=TRUE)
dataset <- opts$dataset

# Pull data from each saved file and grab the name of the fit
folder <- paste("results/",opts$dataset,"/fits",sep="")
load(paste("data/",opts$dataset,".rdata",sep=""))
modelnames <- function(folder) {
  as.vector(sapply(dir(folder),function(x) strsplit(x,".rdata")[[1]][1]))
}
modelatts <- function(model) {
  atts <- strsplit(model,"\\.")[[1]]
  g <- function(x,y) as.numeric(strsplit(x,y)[[1]][2])
  xs <- c("kinit","sm","nb","pshift","deg","trans","collapse")
  a <- lapply(1:length(xs),function(i) g(atts[i],xs[i]))
  names(a) <- xs
  return(a)
}

fs <- list.files(folder,full.names=TRUE)
models <- modelnames(folder)
atts <- lapply(models,function(m) as.data.frame(modelatts(m)))
atts <- do.call(rbind,atts)
atts$model <- models

fits <- lapply(fs,function(f) {
  load(f)
  return(fit)
})
names(fits) <- models

f <- fits[[2]]  # K=1
g <- fits[[8]]

pred.f <- evaluate(A,N,train.ix,test.ix,f,niters=2)
pred.g <- evaluate(A,N,train.ix,test.ix,g,niters=2)

load("results/classroom-16/final/results.rdata")

fllks <- lposterior(f$params,priors)$y
gllks <- lposterior(g$params,priors)$y

tail(f$llks)
tail(g$llks[1:g$iter])

mean(f$llks)
mean(g$llks[1:g$iter]
