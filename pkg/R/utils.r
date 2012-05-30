##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot histograms of the log probability of data and model parameters
##' @param fit brem object
##' @param priors 
##' @return 
##' @author chris
plot.posterior <- function(train,N,fit) {
  
  M <- nrow(train)
  P <- dim(fit$params$beta)[1]
  K <- dim(fit$params$beta)[2]
  ego <- fit$ego*1  # RemStat doesn't want boolean
  s <- new(RemStat,
           train[,1],
           as.integer(train[,2])-1,
           as.integer(train[,3])-1,
           N,M,P,ego)
  s$precompute()  
  s$transform()
  sigma <- fit$params$sigma
  mu <- fit$params$mu
  beta  <- fit$params$beta
  z <- fit$params$z
  if (length(sigma) == 1) {
    mu <- rep(mu,P)
    sigma <- rep(sigma,P)
  }

 # par(mfrow=c(K,K),mar=c(2,2,1,1))
  lbetas <- ws <- ys <- list()
  for (k1 in 1:K) {
    lbetas[[k1]] <- ys[[k1]] <- ws[[k1]] <- list()
    for (k2 in 1:K) {
      lbetas[[k1]][[k2]] <- dnorm(beta[,k1,k2],mu,sigma,log=TRUE)
      k1nodes <- which(z==k1)
      k2nodes <- which(z==k2)
      ys[[k1]][[k2]] <- RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,beta,z-1,s$ptr(),K)
      ws[[k1]][[k2]] <- ys[[k1]][[k2]][-length(ys[[k1]][[k2]])]
     # hist(ys[[k1]][[k2]])
    }
  }

  lsigmas <- dgamma(sigma,fit$priors$sigma$alpha,fit$priors$sigma$beta,log=TRUE)

  tb <- table(factor(z,1:K))
  lzs <- log(tb[z] + fit$priors$alpha)

  par(mfrow=c(3,2))
  hist(unlist(lsigmas),main="Parameters sigma_p",xlab="logprob")
  hist(unlist(lbetas),main="Parameters beta_p,k1,k2",xlab="logprob")
  hist(unlist(ys),main="Observations",xlab="logprob")
  hist(unlist(ws),main="Observations",xlab="logprob")
  hist(unlist(lzs),main="Parameters z")
}

plot.posterior.pc <- function(lp) {
  par(mfrow=c(1,3))
  hist(lp$sigma,main="Parameters sigma_p",xlab="logprob")
  hist(lp$beta,main="Parameters beta_p,k1,k2",xlab="logprob")
  hist(lp$y,main="Observations",xlab="logprob")
}

##' @title 
##' @param A 
##' @param N 
##' @param fit 
##' @param i 
##' @param s 
##' @return 
##' @author chris
eval.llk <- function(edgelist,N,fit,i,s) {
  if (i > nrow(edgelist)) stop("index out of bounds")
  if (i==1) {
    lrm <- brem.lrm.fast.subset(s, fit$z, fit$beta, i)
    llk <- lrm[edgelist[i,2],edgelist[i,3]]
  } else {
    lrm <- brem.lrm.fast.subset(s, fit$z, fit$beta, (i-1):i)
  # Currently written to account for how the following function
  # deals with the first event
    llks <- RemLogLikelihoodVecFromArray(lrm,edgelist[(i-1):i,1],as.integer(edgelist[(i-1):i,2])-1,as.integer(edgelist[(i-1):i,3])-1,N,2)
    llk <- llks[2]
  }
  return(llk)
}

##' @title 
##' @param edgelist 
##' @param lrm 
##' @param i 
##' @return 
##' @author chris
eval.brem <- function(edgelist,lrm,i) {
  N <- dim(lrm)[2]
  sen <- as.integer(edgelist[i,2])
  rec <- as.integer(edgelist[i,3])
  if (i == 1) {
    llk <- lrm[sen,rec]
  } else {
    delta <- c(0,edgelist[i,1] - edgelist[i-1,1])
    a <- sample((1:N)[-sen],1)
    b <- sample((1:N)[-rec],1)
    sen <- c(a,sen)  # arbitrary first event
    rec <- c(b,rec)
    lam <- array(1,dim=c(2,dim(lrm)))
    lam[2,,] <- lrm
    llk <- RemLogLikelihoodVecFromArray(lam,delta,sen-1,rec-1,N,2)
    llk <- llk[2]
  }
  return(llk)
}

##' @title Compute the multinomial log likelihood of a given event
##' @param edgelist Mx3 matrix of (time, sender, receiver)
##' @param lrm current log rate matrix
##' @param i index of event to evaluate
##' @return
##' @export
##' @author chris
eval.mult <- function(edgelist,lrm,i) {
  m <- exp(lrm)
  m <- m/sum(m)
  log(m[edgelist[i,2],edgelist[i,3]])
}

##' @title Compute the rank of a given event
##' @param edgelist Mx3 matrix of (time, sender, receiver)
##' @param lrm current log rate matrix
##' @param i index of event to evaluate
##' @param ... parameters to pass to rank()
##' @return
##' @export
##' @author chris
eval.rank <- function(edgelist,lrm,i,...) {
  N <- dim(lrm)[1]
  rk <- matrix(rank(-lrm,...),N,N)
  rk[edgelist[i,2],edgelist[i,3]]
}

eval.mrank <- function(edgelist,lrm,i,nranks=10) {
  x <- sapply(1:nranks,function(ix) eval.rank(train,lrm,i,ties.method="random"))
  mean(x)
}


##' Compute multinomial log likelihood and ranks of train and test events using a fitted BREM model in an online fashion.
##' @param edgelist Mx3 matrix of (time, sender, receiver)
##' @param N number of nodes
##' @param train.ix indices of edgelist to use as training data
##' @param test.ix indices of edgelist to use as test data
##' @param fit BREM object having elements beta, z, and ego
##' @param ... options to pass to rank()
##' @return list of mllk, rks
##' @export
evaluate <- function(edgelist,N,train.ix,test.ix,fit,niters=NULL,nranks=10) {
  train <- edgelist[train.ix,]
  test  <- edgelist[test.ix,]
  P <- dim(fit$params$beta)[1]
  strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P,as.integer(fit$ego))
  strain$precompute()
  if (fit$transform) strain$transform()
  stest <- new(RemStat,edgelist[,1],as.integer(edgelist[,2])-1,as.integer(edgelist[,3])-1,N,nrow(edgelist),P,as.integer(fit$ego))
  stest$precompute()
  if (fit$transform) stest$transform()

  # If not averaging across samples, use latest sample
  last <- fit$iter #length(fit$samples)
  if (is.null(niters)) {
    iters <- last
  } else {
    iters <- (last-niters):last
  }
  
  # Compute loglikelihoods
  llk <- list(train = rep(0,nrow(train)),
              test  = rep(0,length(test.ix)))
  rks <- mllk <- llk
  for (m in 1:nrow(train)) {
    lrm <- posterior.mean.lrm(strain,fit$samples,m,iters)
#    lrm <- brem.lrm.fast.subset(strain, fit$params$z, fit$params$beta, m)
    lrm <- lrm[1,,]    
    llk$train[m]  <- eval.brem(train,lrm,m)
    mllk$train[m] <- eval.mult(train,lrm,m)
    rks$train[m]  <- eval.rank(train,lrm,m,ties.method="random")    
#    rks$train[m]  <- eval.mrank(train,lrm,m,nranks)
  }
  for (m in 1:length(test.ix)) {
    lrm <- posterior.mean.lrm(stest,fit$samples,test.ix[m],iters)
#    lrm <- brem.lrm.fast.subset(stest, fit$params$z, fit$params$beta, test.ix[m])
    lrm <- lrm[1,,]
    llk$test[m]  <- eval.brem(edgelist,lrm,test.ix[m])
    mllk$test[m] <- eval.mult(edgelist,lrm,test.ix[m])
    rks$test[m] <- eval.rank(edgelist,lrm,test.ix[m],ties.method="random")    
#    rks$test[m]  <- eval.mrank(edgelist,lrm,test.ix[m],nranks)
  }
  return(list(llk=llk,mllk=mllk,rks=rks))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Compute the posterior mean of a log intensity array lambda_ij(t_m)
##' @param strain RemStat pointer 
##' @param samples samples from the posterior
##' @param ix event indices to computel oglikelihood for
##' @param iters indices of the samples list to use
##' @return 
##' @author chris
posterior.mean.lrm <- function(strain,samples,ix,iters) {
  lrms <- lapply(samples[iters],function(s) {
    brem.lrm.fast.subset(strain, s$z, s$beta, ix)
  })
  if (length(iters) > 1) {
    x <- lrms[[1]]
    for (i in 2:length(iters)) x <- x + lrms[[i]]
    lrms <- list(x / length(iters))
  } 
  return(lrms[[1]])
}
posterior.mean.lrm.s <- function(strain,samples,ix,iters) {
  lrms <- lapply(samples[iters],function(s) {
    brem.lrm.fast.subset(strain, s$z, s$beta, ix)
  })
  lrm <- abind(lrms,along=4)
  apply(lrm,1:3,mean)
}

##' @title Same as eval.online, but for a particular baseline method
##' @param edgelist Mx3 matrix of (time, sender, receiver)
##' @param N number of nodes
##' @param train.ix indices of edgelist to use as training data
##' @param test.ix indices of edgelist to use as test data
##' @param model type of baseline (online, marg, or uniform)
##' @param ... 
##' @return 
##' @author chris
evaluate.baseline <- function(edgelist,N,train.ix,test.ix,model="online",nranks=10) {
  if (! model %in% c("online","marg","uniform")) stop("unrecognized baseline")

  train <- edgelist[train.ix,]
  test  <- edgelist[test.ix,]
  eps <- 1
  
  r <- list()
  r$uniform <- matrix(1,N,N)
  
  ## Compute marginals
  b <- table(factor(train[, 2], 1:N), factor(train[, 3], 1:N))
  r$marg <- rowSums(b) %*% t(colSums(b))

  ## Initialize online counts
  r$online <- matrix(0,N,N)

  llk <- list(train = rep(0,nrow(train)),
              test  = rep(0,length(test.ix)))
  rks <- mllk <- llk
  for (m in 1:nrow(train)) {
    if (m > 1) {
      i <- train[m-1,2]
      j <- train[m-1,3]
      r$online[i,j] <- r$online[i,j] + 1
    }    
    lam <- r[[model]] + eps
    lam <- lam/sum(lam)
    lrm <- log(lam * nrow(train)/train[nrow(train),1])
    diag(lrm) <- -Inf
    llk$train[m]  <- eval.brem(train,lrm,m)
    mllk$train[m] <- eval.mult(train,lrm,m)
    rks$train[m]  <- eval.rank(train,lrm,m,ties.method="random")
    #rks$train[m]  <- eval.mrank(train,lrm,m,nranks=nranks)
  }
  for (m in 1:length(test.ix)) {
    i <- edgelist[test.ix[m]-1,2]
    j <- edgelist[test.ix[m]-1,3]
    r$online[i,j] <- r$online[i,j] + 1
    lam <- r[[model]] + eps
    lam <- lam/sum(lam)
    ## use expected rate from training set
    lrm <- log(lam * nrow(train)/train[nrow(train),1])  
    diag(lrm) <- -Inf
    llk$test[m]  <- eval.brem(edgelist,lrm,test.ix[m])
    mllk$test[m] <- eval.mult(edgelist,lrm,test.ix[m])
    rks$test[m]  <- eval.rank(edgelist,lrm,test.ix[m],ties.method="random")
    #rks$test[m]  <- eval.mrank(edgelist,lrm,test.ix[m],nranks=nranks)
  }

  return(list(llk=llk,mllk=mllk,rks=rks))
}
  
##   # Compute multinomial likelihood

## # Compute loglikelihood of each observation
## llk.train <- RemLogLikelihoodVecFromArray(pred$lrm$train,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train))
## llk.test <- RemLogLikelihoodVecFromArray(pred$lrm$test,test[,1],as.integer(test[,2])-1,as.integer(test[,3])-1,N,nrow(test))

##' Compute log intensity arrays for a particular baseline for the training set and test set
##' @param train training event history
##' @param A entire event history
##' @param test.ix indices of A that represent the test set
##' @param fit object as returned by brem.mcmc
##' @export
get.pred.baseline <- function(train,A,N,test.ix,model="online") {
  get.ms <- switch(model,
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
  m.train <- get.ms(train)
  m.train[which(m.train==-Inf)] <- 0
  for (i in 1:nrow(train)) {
    lam <- m.train[i,,] + eps
    diag(lam) <- 0
    m.train[i,,] <- lam/sum(lam)
  }
  m.test <- get.ms(A)
  m.test <- m.test[test.ix,,]
  m.test[which(m.test==-Inf)] <- 0
  for (i in 1:length(test.ix)) {
    lam <- m.test[i,,] + eps
    diag(lam) <- 0
    m.test[i,,] <- lam/sum(lam)
  }
  # Set rates of diagonals to 0
  for (i in 1:nrow(train)) {
    diag(m.train[i,,]) <- 0
  }
  for (i in 1:length(test.ix)) {
    diag(m.test[i,,]) <- 0
  }
  # Get lambda estimates using global rate
  lrm.train <- log(m.train * nrow(train)/train[nrow(train),1])
  lrm.test  <- log(m.test  * nrow(train)/train[nrow(train),1])
  return(list(lrm=list(train=lrm.train,test=lrm.test),m=list(train=m.train,test=m.test)))
}

##' Compute log intensity arrays for the given REM fit for the training set and test set
##' @param train training event history
##' @param A entire event history
##' @param test.ix indices of A that represent the test set
##' @param fit object as returned by brem.mcmc
##' @export
get.pred <- function(train,A,N,test.ix,fit) {
  P <- dim(fit$beta)[1]
  ego <- 0
  strain <- new(RemStat,train[,1],as.integer(train[,2])-1,as.integer(train[,3])-1,N,nrow(train),P,fit$ego)
  strain$precompute()
  stest <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P,ego)
  stest$precompute()
  lrm <- list()
  lrm$train <- brem.lrm.fast(strain, fit$z, fit$beta)
  lrm$test  <- brem.lrm.fast(stest,  fit$z, fit$beta)
  lrm$test  <- lrm$test[test.ix,,]
  
  # Compute multinomial likelihoods
  m.train <- exp(lrm$train)
  m.test <- exp(lrm$test)
  for (i in 1:nrow(train)) {
    m.train[i,,] <- m.train[i,,]/sum(m.train[i,,])
  }
  for (i in 1:length(test.ix)) {
    m.test[i,,] <- m.test[i,,]/sum(m.test[i,,],na.rm=TRUE)
  }
  list(lrm=lrm,m=list(train=m.train,test=m.test))
}

##' Compute the multinomial likelihood for the given array and event history
##' @param p M x N x N array where each slice m is a matrix of probabilities of each dyad
##' @param x M x 3 event history
multinomial.score <- function(p,x) {
  M <- nrow(x)
  r <- rep(0, M)
  for (i in 1:M) {
    r[i] <- p[i,x[i,2],x[i,3]]
  }
  return(r)
}


#' Returns the number of occurrences prior to each event for the observed dyad
#' @x event history
#' @N number of nodes
counts.of.observed <- function(x,N) {
  counts <- ratemat.online(x,N)
  M <- nrow(x)
  r <- rep(0, M)
  for (i in 1:M) {
    r[i] <- counts[i,x[i,2],x[i,3]]
  }
  return(r)
}

plotmat <- function(mat,labels=c("Sender","Receiver"),limits=c(0,max(mat[,3])),color="black") {
  require(ggplot2)
  colnames(mat)[3] <- "value"
  qq <- ggplot(mat,aes(x=X2,y=X1)) +
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low="white",high=color,limits=limits)+
    theme_bw() + labs(x=labels[2],y=labels[1],fill="Probability")+
    coord_equal(ratio=1) + 
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) +
    opts(legend.position = "none",
         panel.grid.minor=theme_blank(),
         panel.grid.major=theme_blank())
  qq
}

get.indices <- function(A,N) {
  # Fast enough for our scale of problems
  # Returns 0-based index
  # Includes M-1 in all lists
  M <- nrow(A)
  xs <- vector("list",N)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    xs[[i]] <- c(xs[[i]],m - 1)
    xs[[j]] <- c(xs[[j]],m - 1)
  }
  for (i in 1:N) 
    xs[[i]] <- c(xs[[i]],M - 1)
  return(xs)
}

#' Compute the rank of the each event in the observed edgelist given an array of intensities
#' @edgelist Mx3 matrix of event times, senders and recipients
#' @ratemats MxNxN array of intensities.  
ranks.fast <- function(edgelist,ratemats,...) {
  M <- nrow(edgelist)
  N <- dim(ratemats)[2]
  r <- rep(0,M)
  for (m in 1:M) {
    i <- edgelist[m,2]
    j <- edgelist[m,3]
    k <- N*(j-1) + i
    r[m] <- rank(ratemats[i,,],...)[k]
  }
  return(r)
}

ranks <- function(edgelist,ratemats,...) {
  M <- nrow(edgelist)
  n <- dim(ratemats)[2]
  r <- rep(0,M)
  for (i in 1:M) {
    if (length(dim(ratemats))==2) rmat <- ratemats
    else rmat <- ratemats[i,,]
    rk <- matrix(rank(rmat,...),n,n)
    r[i] <- rk[edgelist[i,2],edgelist[i,3]]
  }
  return(r)
}

#' Compute the recall metric using a vector of ranks
#' @rs vector of observed ranks
#' @top vector of cutpoints to compute
recall <- function(rs,top=c(1:20)) {
  data.frame(k=top,recall=sapply(top,function(k) mean(rs <= k)))
}

dyad.ps <- function(A,N,ego=TRUE) {
  x <- array(0,c(N,N,7))
  mp <- matrix(0,N,N)
  M <- nrow(A)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    # Find last time lambda_ij changed
    if (mp[i,j]!=0) {
      a <- A[mp[i,j],2]
      b <- A[mp[i,j],3]
      x[i,j,] <- x[i,j,] + pshift(i,j,a,b)
    }
    mp[i,] <- m
    mp[j,] <- m
    if (ego == 0) {
      mp[,i] <- m
      mp[,j] <- m
    }
  }
  dimnames(x)[[3]] <- c("AB-BA","AB-BY","AB-XA","AB-XB","AB-AY","AB-AB","Total")
  return(x)
}

node.ps <- function(A,N) {
  require(relevent)
  x <- matrix(0,N,2)
  p <- matrix(0,N,5)
  M <- nrow(A)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    a <- x[i,1]
    b <- x[i,2]
    p[i,] <- p[i,] + pshift(i,j,a,b)
    a <- x[j,1]
    b <- x[j,2]
    p[j,] <- p[j,] + pshift(i,j,a,b)
    
    # Update last event for i and j
    x[i,] <- x[j,] <- c(i,j)
  }
  colnames(p) <- c("AB-BA","AB-BY","AB-XA","AB-XB","AB-AY")
  return(p)
}

RLogLambda <- function(i,j,a,b,beta,px) {
  lam <- 0
  if (px[1]) {
    lam <- lam + beta[1]
  }
  lam <- lam + sum(px[-1] * beta[-1] * pshift(i,j,a,b))
  return(lam)
}
pshift <- function(i,j,a,b) {
  x <- rep(0,7)
  if (i!=a & i==b & j==a & j!=b) { # ab-ba
    x[1] <- 1
  }
  if (i!=a & i==b & j!=a & j!=b) { # ab-by
    x[2] <- 1
  }
  if (i!=a & i!=b & j==a & j!=b) { # ab-xa
    x[3] <- 1
  }
  if (i!=a & i!=b & j!=a & j==b) { # ab-xb
    x[4] <- 1
  }
  if (i==a & i!=b & j!=a & j!=b) { # ab-ay
    x[5] <- 1
  }
  if (i==a & i!=b & j!=a & j==b) { # ab-ab
    x[6] <- 1
  }
  x[7] <- 1
  return(x)
}
block.ps <- function(A,z) {
  require(relevent)
  s <- A[,2]
  r <- A[,3]
  zs <- sort(unique(z))
  setup <- expand.grid(i=zs,j=zs)
  ds <- lapply(1:nrow(setup),function(k) {
    ix <- which(z[s] == setup$i[k] & z[r]==setup$j[k])
    pshifts <- accum.ps(A[ix,])[length(ix),c(1,3,8,9,10,13)]
    data.frame(i=setup$i[k],j=setup$j[k],pshift=names(pshifts),value=pshifts)
  })
  ds <- do.call(rbind,ds)
  rownames(ds) <- c()
  return(ds)
}
rcatlog <- function(ps) {
  ps <- ps - max(ps)
  ps <- exp(ps)
  ps <- ps/sum(ps)
  ps <- cumsum(ps)
  u <- runif(1)
  for (i in 1:length(ps)) {
    if (u < ps[i]) return(i)
  }
}
plot.blockmodel <- function(B,z) {
  nmap <- order(z)
  s <- match(B[,2],nmap)
  r <- match(B[,3],nmap)
  plot(s,r,pch=".",xlab="sender",ylab="recipient")
  cpoints <- sapply(1:K,function(k) min(which(sort(z)==k)))
  abline(v=cpoints[-1])
  abline(h=cpoints[-1])
}

#' Compute an (M x N x N) array of previous events for each dyad.  
#' @edgelist event history matrix
#' @n total number of nodes
ratemat.online <- function(edgelist,n) {
  M <- nrow(edgelist)
  rms <- array(0,c(M,n,n))
  for (i in 2:M) {
    rms[i:M,edgelist[i-1,2],edgelist[i-1,3]] <- rms[i:M,edgelist[i-1,2],edgelist[i-1,3]] + 1
    diag(rms[i,,]) <- -Inf
  }
  return(rms)
}

#' Compute an (M.test x N x N) array of intensities for each dyad using the marginals from a training set.
#' @train event history
#' @test  event history to compute the array for
#' @N     number of nodes
ratemat.from.marginals <- function(train,test,N) {
  x <- table(factor(train[,2],1:N),factor(train[,3],1:N))
  rowrates <- rowSums(x)
  colrates <- colSums(x)
  r <- rowrates %*% t(colrates)
  r <- r/sum(r)
  diag(r) <- -Inf
  M <- nrow(test)
  lrm <- array(0,c(M,N,N))
  for (i in 1:M) lrm[i,,] <- r
  return(lrm)
}


#' Full lograte array using computeLambdaFast.  lrm[m,i,j] is lambda_ij prior m'th event.  
#' All lambda=0 for m=0.
lrm_slow <- function(beta,z,s,M,N,K,P) {
  lrms <- array(0,c(M,N,N))
  for (m in 2:M) {
    for (i in 1:N) {
      for (j in 1:N) {
        if (i!=j) {
          x <- s$get_s(m-2,i-1,j-1)  # why m-2?  because this should be for the rateprior to current event m-1
          lrms[m,i,j] <- compute_lambda_fast(i,j,z[i],z[j],x,beta,N,K,P)
        } 
      }
    }
  }
  return(lrms)
}

llk_indiv <- function(a,lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 0:(M-1)) {
    i <- sen[m+1]
    j <- rec[m+1]
    if (i==a | j==a | m==(M-1)) {
      llks[m+1] <- lrm[m+1,i+1,j+1]
      for (r in 0:(N-1)) {
        if (r != i & r!=j) {
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * sum(exp(lrm[m+1,i+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * sum(exp(lrm[m+1,r+1,i+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * sum(exp(lrm[m+1,j+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * sum(exp(lrm[m+1,r+1,j+1]))
        }
      }
      llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * sum(exp(lrm[m+1,i+1,r+1]))
      llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,r+1,i+1]))
    }
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  return(llks)
}



test_taus_from_s <- function(times,sen,rec,N,M,ego) {
  s <- new(RemStat,times,sen,rec,N,M,ego)
  s$precompute()
  taus <- array(0,c(M,N,N))
  for (m in 0:(M-1)) {
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i != j) {
          taus[m+1,i+1,j+1] <- s$get_tau(m,i,j)
        }
      }
    }
  }
  return(taus)
}

test_taus <- function(lrm,times,sen,rec,M,N) {
  mp <- matrix(0,N,N)
  llks <- matrix(0,N,N)
  M <- length(times)
  taus <- array(0,c(M,N,N))
  for (m in 0:(M-2)) {
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i!=j) {
          taus[m+1,i+1,j+1] <- times[mp[i+1,j+1]+1]
        }
      }
    }
    i <- sen[m+1]
    j <- rec[m+1]
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
    
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        taus[m+1,i+1,j+1] <- times[mp[i+1,j+1]+1]
      }
    }
  }
  return(taus)
}
llk_fast_last <- function(lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- matrix(0,N,N)
  taus <- matrix(0,N,N)
  for (m in 1:(M-2)) {
    i <- sen[m]
    j <- rec[m]
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        taus[i+1,j+1] <- times[mp[i+1,j+1]+1]
        llks[i+1,j+1] <- (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  return(list(taus=taus,llks=llks,mp=mp))
}


# 0 based indexing on sen and rec
RRemLogLikelihoodFromArraySlow <- function(lrm,times,sen,rec,N,M) {
  sen <- sen+1
  rec <- rec+1
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 1:M) diag(lrm[m,,]) <- -Inf
  llks[1] = lrm[1,sen[1],rec[1]]
  for (m in 2:M) {
    llks[m] <- lrm[m,sen[m],rec[m]] - (times[m]-times[m-1]) * sum(exp(lrm[m,,]))
  }
  return(llks)
}

RRemLogLikelihoodFromArrayMult <- function(lrm,times,sen,rec,M,N) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  llks[1] <- lrm[1,sen[1]+1,rec[1]+1]
  for (m in 1:(M-2)) {
    i <- sen[m+1]
    j <- rec[m+1]
    llks[m+1] <- lrm[m+1,i+1,j+1]
    for (r in 0:(N-1)) {
      if (r != i & r != j) {
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * exp(lrm[m+1,i+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * exp(lrm[m+1,r+1,i+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * exp(lrm[m+1,j+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * exp(lrm[m+1,r+1,j+1])
      }
    }
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * exp(lrm[m+1,i+1,j+1])
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * exp(lrm[m+1,j+1,i+1])
    
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        #cat(i," ",j," ",times[mp[i+1,j+1]+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  return(llks)
}

RRemLogLikelihoodFromArrayAlt <- function(lrm,times,sen,rec,N,M) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  llks[1] <- lrm[1,sen[1]+1,rec[1]+1]
  for (m in 0:(M-2)) {
    i <- sen[m+1]
    j <- rec[m+1]
    llks[m+1] <- lrm[m+1,i+1,j+1]
    for (r in 0:(N-1)) {
      if (r != i & r != j) {
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * exp(lrm[m+1,i+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * exp(lrm[m+1,r+1,i+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * exp(lrm[m+1,j+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * exp(lrm[m+1,r+1,j+1])
      }
    }
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * exp(lrm[m+1,i+1,j+1])
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * exp(lrm[m+1,j+1,i+1])
    
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  llks[m+1] <- lrm[m+1,sen[m+1]+1,rec[m+1]+1]
  other <- 0
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        other <- other - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  llks[m+1] <- llks[m+1] + other
  return(llks)
}

precomputeTau <- function(A,N) {
  M <- nrow(A)
  sen <- A[,2]
  rec <- A[,3]
  tau <- array(0,c(M,N,N))
  for (m in 1:(M-1)) {
    i <- sen[m]
    j <- rec[m]
    tau[(m+1):M,,i] <- m-1
    tau[(m+1):M,,j] <- m-1
    tau[(m+1):M,i,] <- m-1
    tau[(m+1):M,j,] <- m-1
  }
  return(tau)
}

