#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(brem))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-m", "--model"), default=NULL,
              help=""))
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

options(verbose=FALSE)

# Pull data from each saved file and grab the name of the fit
load(paste("data/",opts$dataset,".rdata",sep=""))

test.ix <- (1:nrow(A))[-(1:nrow(train))]

if (!opts$model %in% c("uniform","online","marg")) {
  if (opts$model=="countsonly") {
    beta <- rep(0,13)
    beta[12] <- 1
    beta <- array(beta,c(13,1,1))
    z <- rep(1,N)
    res <- list(z=z,beta=beta)
  } else {
    f <- paste("results/",opts$dataset,"/",opts$model,".rdata",sep="")
    load(f)
  }
  cat("precomputing\n")
  P <- 13
  strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P)
  strain$precompute()
  stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
  stest$precompute()
  cat("lambdas (train)\n")
  lrm.train <- brem.lrm.fast(nrow(train), strain, res$z, res$beta)
  cat("lambdas (test)\n")
  lrm.test  <- brem.lrm.fast2(stest, res$z, res$beta,test.ix)
  # Compute multinomial likelihoods
  m.train <- exp(lrm.train)
  m.test  <- exp(lrm.test)
  for (i in 1:nrow(train)) {
    m.train[i,,] <- m.train[i,,]/sum(m.train[i,,],na.rm=TRUE)
  }
  for (i in 1:length(test.ix)) {
    m.test[i,,] <- m.test[i,,]/sum(m.test[i,,],na.rm=TRUE)
  }
} else {
    get.ms <- switch(opts$model,
     "uniform" = function(x) array(1,c(nrow(x),N,N)),
     "online"  = function(x) ratemat.online(x,N),
     "marg"    = function(x) {
       b <- table(factor(train[, 2], 1:N), factor(train[, 3], 1:N))
       rowrates <- rowSums(b)
       colrates <- colSums(b)
       r <- rowrates %*% t(colrates)
       ma <- array(0, c(nrow(x), N, N))
       for (i in 1:nrow(x)) ma[i, , ] <- r
       return(ma)
    })

    eps <- 1  # smoothing
    cat("lambdas (train)\n")
    m.train <- get.ms(train)
    m.train[which(m.train==-Inf)] <- 0
    for (i in 1:nrow(train)) {
      lam <- m.train[i,,] + eps
      m.train[i,,] <- lam/sum(lam)
    }
    m.test <- get.ms(A)
    m.test <- m.test[test.ix,,]
    m.test[which(m.test==-Inf)] <- 0
    for (i in 1:length(test.ix)) {
      lam <- m.test[i,,] + eps
      m.test[i,,] <- lam/sum(lam)
    }
    # Get lambda estimates using global rate
    lrm.train <- log(m.train * nrow(train)/train[nrow(train),1])
    lrm.test  <- log(m.test  * nrow(train)/train[nrow(train),1])
}

# Compute multinomial socre for each observation
#' Multinomial probability of observed data
#' @m multinomial probability array (M x N x N)
#' @x observed event history
multinomial.score <- function(m,x) {
  M <- nrow(x)
  r <- rep(0, M)
  for (i in 1:M) {
    r[i] <- m[i,x[i,2],x[i,3]]
  }
  return(r)
}
llkm.train <- log(multinomial.score(m.train,train))
llkm.test  <- log(multinomial.score(m.test, test))

# Compute loglikelihood of each observation
llk.train <- loglikelihood_vec_from_lrm(lrm.train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
llk.test <- loglikelihood_vec_from_lrm(lrm.test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

# Compute rank of each observation
cat("ranks (train)\n")
rk.train <- ranks(train,-lrm.train,ties.method="random")
cat("ranks (test)\n")
rk.test  <- ranks(test, -lrm.test, ties.method="random")

# Save
dir.create(paste("results/",opts$dataset,"/llks/",sep=""))
save(llkm.train,llkm.test,llk.train,llk.test,opts,file=paste("results/",opts$dataset,"/llks/",opts$model,".rdata",sep=""))
dir.create(paste("results/",opts$dataset,"/ranks/",sep=""))
save(rk.train,rk.test,opts,file=paste("results/",opts$dataset,"/ranks/",opts$model,".rdata",sep=""))

