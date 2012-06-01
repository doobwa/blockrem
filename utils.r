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

# TODO: Just input a tuple...?
progress <- function(f,x,num.match=2) {
  if (!file.exists(f)) {
    current <- matrix(x,nr=1)
  } else {
    load(f)
    ix <- which(current[,1] == x[1] & current[,2] == x[2])
    if (length(ix) == 0) {
      current <- rbind(current,x)
    } else {
      current[ix,] <- x
    }
  }
  save(current,file=f)
}

get_progress <- function(f) {
  if (!file.exists(f)) stop("progress file does not exist")
  load(f)
  x <- data.frame(dataset=current[,1],model=current[,2],action=current[,3],time=as.numeric(current[,4]))
  rownames(x) <- NULL
  x$model <- as.character(x$model)
  ma <- lapply(unique(x$model),function(m) data.frame(model=m,modelatts(m)))
  ma <- do.call(rbind,ma)
  merge(x,ma,by="model")
}

## progress <- function(dataset,model,action,notes=NULL,f="progress.rdata") {
##   if (!file.exists(f)) {
##     current <- list()
##     save(current,file=f)
##   } 
##   load(f)
##   if (length(current[[dataset]]) == 0) {
##     current[[dataset]] <- list()
##   }
##   if (is.null(notes)) {
##     current[[dataset]][[model]] <- action
##   } else {
##     current[[dataset]][[model]] <- list("action"=action,"notes"=notes)
##   }
##   save(current,file=f)
## }

## get_progress <- function(f="progress.rdata") {
##   if (!file.exists(f)) stop("progress file does not exist")
##   load(f)
##   x <- melt(current)
##   x <- dcast(x,L1 + L2 ~ L3)
##   colnames(x) <- c("dataset","model","action","time")
##   ma <- lapply(unique(x$model),function(m) data.frame(model=m,modelatts(m)))
##   ma <- do.call(rbind,ma)
##   merge(x,ma,by="model")
## }

