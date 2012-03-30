#' General purpose Metropolis-Hastings sampler
#' @param current
#' @param likelihood
#' @param olp
#' @param sd

#' @return ...
mh <- function(current,lposterior,olp=NULL,sd=.1) {
  if (is.null(olp)) {
    olp <- attr(current,"log.density") <- lposterior(current)
  }
  cand <- current + rnorm(length(current),0,sd)
  clp <- lposterior(cand)
  if (clp - olp > log(runif(1))) {
    current <- cand
    attr(current,"log.density") <- clp
  }
  return(current)
}

slice <- function(current,lposterior,olp=NULL,m=20) {
  if (is.null(olp)) {
    olp <- attr(current,"log.density") <- lposterior(current)
  }
  lpost <- function(val) {
    current[p] <- val
    lposterior(current)
  }
  for (p in 1:length(current)) {
    val <- uni.slice.alt(current[p],lpost,gx0=olp,m=m)
    current[p] <- val
    olp <- attr(val,"log.density")
  }
  attr(current,"log.density") <- olp
  return(current)
}

# Changed U and grad_U to not be negative llk, etc.
# Added log.density attributes
hmc <- function (current_q, lposterior, epsilon, L) 
{
  U      <- function(x) -lposterior(x)
  grad_U <- function(x) -lposterior(x,lp=FALSE)
    q = current_q
    p = rnorm(length(q), 0, 1)
    current_p = p
    p = p - epsilon * grad_U(q)/2
    for (i in 1:L) {
        q = q + epsilon * p
        if (i != L) 
            p = p - epsilon * grad_U(q)
    }
    p = p - epsilon * grad_U(q)/2
    p = -p
    current_U = U(current_q)
    current_K = sum(current_p^2)/2
    proposed_U = U(q)
    proposed_K = sum(p^2)/2
    attr(current_q,"log.density") <- -current_U
    attr(q,"log.density") <- -proposed_U
    if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
        return(q)
    }
    else {
        return(current_q)
    }  
}

lposterior <- function(A,s,beta,z,px,k1,k2,priors=list(beta=list(mu=0,sigma=1)),lp=TRUE,lgrad=FALSE) { 
  if (lp) {
    res <- block(A,s,beta,z,k1,k2,px,priors)
    if (lgrad) attr(res,"lgrad") <- block.grad(A,s,beta,z,k1,k2,px)
  } else {
    res <- block.grad(A,s,beta,z,k1,k2,px)
  }
  return(res)
}

gibbs.collapsed <- function(ix,beta,z,s,nextra=0,
                            priors=list(alpha=1,beta=list(mu=0,sigma=1))) {

  # count number of unused clusters
  P <- dim(beta)[1]
  K <- dim(beta)[2]
  counts <- table(factor(z,1:K))

  # decide which clusters will be resampled (at least empty clusters, at most nextra)
  tosample <- which(counts==0)
  nextra <- nextra - length(tosample)
  if (nextra > 0) {
    tosample <- c(tosample, (K+1):(K+nextra))
  }
  if (length(tosample) == 0) stop("No new clusters proposed")
  K <- K + nextra
  counts <- table(factor(z,1:K))
  
  # augment beta with draws from prior
  b <- array(0,c(P,K,K))
  b[,-tosample,-tosample] <- beta
  for (k in tosample) {
    for (l in 1:K) {
      b[which(px==1),k,l] <- rnorm(sum(px),priors$beta$mu,priors$beta$sigma)
      b[which(px==1),l,k] <- rnorm(sum(px),priors$beta$mu,priors$beta$sigma)
    }
  }
  
  ys <- matrix(0,N,K)
  for (a in 1:N) {
    counts[z[a]] <- counts[z[a]] - 1
    for (k in 1:K) {
      z[a] <- k
      ys[a,k] <- RemLogLikelihoodActorPc(a-1,b,z-1,s,K) + log(counts[k] + priors$alpha)
    }
    z[a] <- rcatlog(ys[a,])
    counts[z[a]] <- counts[z[a]] + 1
  }

  newmembers <- sum(table(factor(z,1:K))[tosample])
  list(beta=b,z=z,probs=ys,newmembers=newmembers)
}

mcmc <- function(A,N,K,px,method,niter=100,beta=NULL,z=NULL,priors=list(beta=list(mu=0,sigma=1)),verbose=TRUE,...) {
  
  if (sum(px) == 0)
    stop("No parameters selected.")

  ## precompute statistics
  s <- new(brem:::RemStat,A[,1],A[,2]-1,A[,3]-1,N,nrow(A),length(px))
  s$precompute()
  s$transform()
  
  lp <- rep(0,niter)
  param <- array(0,c(niter,dim(beta)))
  for (iter in 1:niter) {
    ## loop through blocks
    for (k1 in 1:K) {
      for (k2 in 1:K) {
        ## sample new parameter values for this block
        b <- beta[which(px==1),k1,k2]
        beta[which(px==1),k1,k2] <- method(b,lposterior,...)
      }
    }
    lp[iter] <- brem.lpost.fast(A,N,K,z,s,beta,priors)
    param[iter,,,] <- beta
    if (verbose) {
      cat(lp[iter],"\n")
    }
  }
  return(list(param=param,lp=lp))
}

block <- function(A,s,beta,z,k1,k2,px,priors=list(beta=list(mu=0,sigma=1))) { 
  beta[which(px==0),,] <- 0
  M <- nrow(A)
  K <- dim(beta)[2]
  zs <- z[A[,2]]
  zr <- z[A[,3]]
  ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
  if (length(ix) == 0) return(0)
  llk <- sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ix))
  if (!is.null(priors)) {
    llk <- llk + sum(dnorm(unlist(beta[which(px==1),k1,k2]),priors$beta$mu,priors$beta$sigma,log=TRUE))
  }
  return(llk)
}

block.grad <- function(A,s,beta,z,k1,k2,px,priors=list(beta=list(mu=0,sigma=1))) {   
  M <- nrow(A)
  zs <- z[A[,2]]
  zr <- z[A[,3]]
  ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
  if (length(ix) == 0) return(0)
  lg <- RemGradientPcSubset(beta,z-1,s$ptr(),K,ix,which(px==1)-1)
  lg <- lg[which(px==1)]
  if (!is.null(priors)) {
    lg <- lg - 2 * (beta[which(px==1),k1,k2] - priors$beta$mu)/priors$beta$sigma
  }
  return(lg)
}
## RRemGradient <- function(lrm,times,sen,rec,N,M,P) {
  
##   s <- InitializeStatisticsArray(N,P);
##   for (m in 2:M) {
    
##     i <- A[m-1,2]
##     j <- A[m-1,3]
    
##     s <- UpdateStatisticsArray(s,m-1,i-1,j-1,N,P)
##   mp <- matrix(0,N,N)
##   llks <- rep(0,M)
##   llks[1] <- lrm[1,sen[1]+1,rec[1]+1]
##   for (m in 0:(M-2)) {
##     i <- sen[m+1]
##     j <- rec[m+1]
##     llks[m+1] <- lrm[m+1,i+1,j+1]
##     for (r in 0:(N-1)) {
##       if (r != i & r != j) {
##         llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * exp(lrm[m+1,i+1,r+1])
##         llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * exp(lrm[m+1,r+1,i+1])
##         llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * exp(lrm[m+1,j+1,r+1])
##         llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * exp(lrm[m+1,r+1,j+1])
##       }
##     }
##     llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * exp(lrm[m+1,i+1,j+1])
##     llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * exp(lrm[m+1,j+1,i+1])
    
##     mp[i+1,] <- m
##     mp[,i+1] <- m
##     mp[j+1,] <- m
##     mp[,j+1] <- m
##   }
##   m <- M-1
##   llks[m+1] <- lrm[m+1,sen[m+1]+1,rec[m+1]+1]
##   other <- 0
##   for (i in 0:(N-1)) {
##     for (j in 0:(N-1)) {
##       if (i!=j) {
##         other <- other - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
##       }
##     }
##   }
##   llks[m+1] <- llks[m+1] + other
##   return(llks)
## }
