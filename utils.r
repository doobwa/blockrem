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

progress <- function(dataset,model,action,f="progress.rdata") {
  if (!file.exists(f)) {
    current <- list()
    save(current,file=f)
  } 
  load(f)
  if (length(current[[dataset]]) == 0) {
    current[[dataset]] <- list()
  }
  current[[dataset]][[model]] <- action
  save(current,file=f)
}

get_progress <- function(f="progress.rdata") {
  load(f)
  x <- melt(current)
  colnames(x) <- c("action","model","dataset")
  ma <- lapply(unique(x$model),function(m) data.frame(model=m,modelatts(m)))
  ma <- do.call(rbind,ma)
  merge(x,ma,by="model")
}

library(testthat)
test_that("view progress makes correct list", {
  f <- "tmp.rdata"
  file.remove(f)
  progress("a","b","start",f=f)
  load(f)
  expect_that(current, equals(list(a = list(b = "start"))))
  progress("a","c","start",f=f)
  load(f)
  expect_that(current, equals(list(a = list(b = "start", c = "start"))))
  progress("a","b","end",f=f)
  load(f)
  expect_that(current, equals(list(a = list(b = "end", c = "start"))))  
  file.remove(f)
})
