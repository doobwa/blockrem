

##' Split phi for cluster k into k and l
split_phi <- function(phi.merge,k,l,priors) {
  P <- dim(phi.merge)[1]
  K <- dim(phi.merge)[2]
  phi.split <- phi.merge
  for (r in 1:K) {
    phi.split[,l,r] <- rnorm(P,phi.merge[,k,r],priors$sigma)
    phi.split[,r,l] <- rnorm(P,phi.merge[,r,k],priors$sigma)
  }
  phi.split[,l,l] <- rnorm(P,phi.merge[,k,k],priors$sigma)
  phi.split[,k,k] <- rnorm(P,phi.merge[,k,k],priors$sigma)
  return(phi.split)
}

##' Merge cluster l into k
merge_phi <- function(phi,k,l,priors) {
  P <- dim(phi)[1]
  K <- dim(phi)[2]
  phi[,k,k] <- rnorm(P,.5 * (phi[,k,k] + phi[,l,l]),priors$sigma)
  if (K > 2) {
    rs <- (1:K)[-c(k,l)]
    for (r in rs) {
      phi[,k,r] <- rnorm(P,.5 * (phi[,k,r] + phi[,l,r]),priors$sigma)
      phi[,r,k] <- rnorm(P,.5 * (phi[,r,k] + phi[,r,l]),priors$sigma)
    }
  }
  phi <- sample_cluster_from_prior(phi,l,priors)
  return(phi)
}

##' Draw cluster l from the prior
add_cluster <- function(phi) {
  P <- dim(phi)[1]
  K <- dim(phi)[2]
  b <- array(0,c(P,K+1,K+1))
  b[,1:K,1:K] <- phi
  return(b)
}

remove_cluster <- function(phi,l) {
  if (any(l > dim(phi)[2] | l==0)) stop("this cluster id does not exist")
  return(phi[,-l,-l,drop=FALSE])
}

sample_cluster_from_prior <- function(phi,l,priors) {
  P <- dim(phi)[1]
  K <- dim(phi)[2]
  if (l > K | l==0) stop("this cluster id does not exist")
  for (r in 1:K) {
    phi[,l,r] <- rnorm(P,priors$phi$mu,priors$phi$sigma)
    phi[,r,l] <- rnorm(P,priors$phi$mu,priors$phi$sigma)
  }
  phi[,l,l] <- rnorm(P,priors$phi$mu,priors$phi$sigma)
  return(phi)
}

##' Return normalized vector of probabilities of log values
##' @param ps log values
##' @return normalized probabilities
lnormalize <- function(ps) {
  ps <- ps - max(ps)
  ps <- exp(ps)
  return(ps/sum(ps))
}

gibbs <- function(phi,z,S,llk_node,N,priors) {

  # Randomly assign nodes that are currently either in k or l
  K <- dim(phi)[2]
  z.init <- z
  counts <- table(factor(z,1:K))
  ys <- matrix(0,N,K)
  for (a in S) {
    counts[z[a]] <- counts[z[a]] - 1
    for (j in 1:K) {
      z[a] <- j
      ys[a,j] <- llk_node(a,phi,z,priors) + log(counts[j] + priors$alpha)
    }
    z[a] <- sample(1:K,1,prob=lnormalize(ys[a,])) #rcatlog(ys[a,])
    counts[z[a]] <- counts[z[a]] + 1
  }
  ys[S,] <- t(apply(ys[S,,drop=FALSE],1,lnormalize))
  q <- sum(log(ys[cbind(S,z[S])]))
  list(z=z,transition=q,ys=ys,z.init=z.init)
}


##' If prob.only, computes the probability of a restricted Gibbs scan
##' p(z|phi) = prod_{a in S} p(z_a|z_prev, phi).
##' Randomly initializes z_a's.  If prob. only, set z_a to provided z.
##' Otherwise, Gibbs sample p(z|phi). Returns z and the transition prob.
gibbs_restricted <- function(phi,z,S,k,l,llk_node,priors,prob.only=FALSE) {

  # Randomly assign nodes that are currently either in k or l
  K <- dim(phi)[2]
  clusters <- c(k,l)
#  S <- which(z %in% clusters)
  z.init <- z
  z[S] <- sample(clusters,length(S),replace=TRUE)
  counts <- table(factor(z,1:K))
  ys <- matrix(0,N,K)
  for (a in S) {
    counts[z[a]] <- counts[z[a]] - 1
    for (j in clusters) {
      z[a] <- j
      ys[a,j] <- llk_node(a,phi,z,priors) + log(counts[j] + priors$alpha)
    }
    if (!prob.only) {
      z[a] <- sample(clusters,1,prob=lnormalize(ys[a,clusters])) #rcatlog(ys[a,])
    } else {
      z[a] <- z.init[a]
    }
    counts[z[a]] <- counts[z[a]] + 1
  }
#  browser()
  ys[S,clusters] <- t(apply(ys[S,clusters],1,lnormalize))
  q <- sum(log(ys[cbind(S,z[S])]))
  list(z=z,transition=q,ys=ys,z.init=z.init)
}


##' Computes the ratio q(phi^m -> phi^s)/p(phi^s)
pm2ps <- function(phi.merge,phi.split,k,l,K,priors) {
  q <-
    dnorm(phi.split[,k,k],phi.merge[,k,k],priors$sigma,log=TRUE) +
    dnorm(phi.split[,l,l],phi.merge[,k,k],priors$sigma,log=TRUE) +
    dnorm(phi.split[,k,l],phi.merge[,k,k],priors$sigma,log=TRUE) +
    dnorm(phi.split[,l,k],phi.merge[,k,k],priors$sigma,log=TRUE) -
    dnorm(phi.split[,k,k],priors$phi$mu,priors$phi$sigma,log=TRUE) -
    dnorm(phi.split[,l,l],priors$phi$mu,priors$phi$sigma,log=TRUE) -
    dnorm(phi.split[,k,l],priors$phi$mu,priors$phi$sigma,log=TRUE) -
    dnorm(phi.split[,l,k],priors$phi$mu,priors$phi$sigma,log=TRUE)

  if (K > 2) {
    rs <- (1:K)[-c(k,l)]
    for (r in rs) {
      q <- q +
        dnorm(phi.split[,r,l],phi.merge[,r,k],priors$sigma,log=TRUE) +
        dnorm(phi.split[,l,r],phi.merge[,k,r],priors$sigma,log=TRUE) -
        dnorm(phi.split[,r,l],priors$phi$mu,priors$phi$sigma,log=TRUE) -
        dnorm(phi.split[,l,r],priors$phi$mu,priors$phi$sigma,log=TRUE)
    }
  }
  return(sum(q))
}

##' Computes the ratio q(phi^s -> phi^m)/p(phi^m)
ps2pm <- function(phi.split,phi.merge,k,l,K,priors) {
  q <-
    dnorm(phi.merge[,k,k],.5*(phi.split[,k,k] + phi.split[,l,l]),priors$sigma,log=TRUE) -
  dnorm(phi.merge[,k,k],priors$phi$mu,priors$phi$sigma,log=TRUE)
  return(sum(q))
}

splitmerge <- function(phi,z,lposterior,llk_node,priors,sigma=.1,verbose=TRUE) {
  N <- length(z)
  K <- dim(phi)[2]
  init <- split <- merge <- list(z=z,phi=phi)
  ij <- sample(1:N,2,replace=FALSE)
  i <- ij[1]; j <- ij[2]
  k <- init$z[i]
  l <- init$z[j]
  S <- which(z %in% z[ij])
  q <- list()
  if (k == l) {
#    browser()
    l <- K+1
    split$phi <- add_cluster(split$phi)
    split$phi <- split_phi(split$phi,k,l,priors)
    h <- gibbs_restricted(split$phi,split$z,S,k,l,llk_node,priors,prob.only=FALSE)
    split$z <- h$z    
  } else {
    merge$z[S] <- k
    merge$phi <- merge_phi(split$phi,k,l,priors)
    h <- gibbs_restricted(merge$phi,split$z,S,k,l,llk_node,priors,prob.only=TRUE)
  }
  lp <- list()
  lp$split <- lposterior(split$phi,split$z,priors)
  lp$merge <- lposterior(merge$phi,merge$z,priors)
  lp$zm2zs <- h$transition
#  if (lp$zm2zs < -200) browser()
  lp$pm2ps <- pm2ps(merge$phi,split$phi,k,l,K,priors)
  lp$ps2pm <- ps2pm(split$phi,merge$phi,k,l,K,priors)
  lp$pzs2pzm <- log(factorial(sum(split$z[S] == k) ) * # TODO: Fix -1 issue
    factorial(sum(split$z[S] == l) ) / factorial(length(S) ))
  alpha <- lp$split - lp$merge + 0 - lp$zm2zs + lp$pzs2pzm + lp$ps2pm - lp$pm2ps
  
  initial <- list(type="nomove",phi=init$phi,z=init$z)
  if (init$z[i] == init$z[j]) {
    proposed <- list(type="split",phi=split$phi,z=split$z,lprob=alpha)
  } else {
    proposed <- list(type="merge",phi=merge$phi,z=merge$z,lprob=-alpha)
  }
  
  if (runif(1) < exp(proposed$lprob)) {
    final <- proposed
    final$accepted <- TRUE
    cat("move accepted\n")
  } else {
    final <- initial
    final$accepted <- FALSE
    cat("move rejected\n")
  }

  if (verbose) {
    type <- ifelse(init$z[i] == init$z[j],"split","merge") 
    cat("move:",type," clusters:",init$z[ij],"\n")
    cat("initial z:",init$z,"\n")
    cat("merge z:",merge$z,"\n")
    cat("split z:",split$z,"\n")
    cat("lprob:",proposed$lprob,"\n")
    #print(lapply(phi[c("init","merge","split")],function(p) p[1,,]))
    ## cat("split",lp$split,"\n")
    ## cat("merge",lp$merge,"\n")
    print(unlist(lp))
  }
  list(phi=phi,chosen=c(i,j),alpha=alpha,z=z,q=q,final=final)
}



## Sketch out mcmc routine
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param lposterior function taking phi, z, and priors
##' @param llk_node function taking a, phi, z, and priors
##' @param priors list of priors for phi, sigma, alpha
##' @param N number of actors
##' @param P number of parameters
##' @param K number of clusters
##' @param niter 
##' @param do.sm perform split-merge moves
##' @param num.extra propose extra clusters at each iteration
##' @param verbose 
##' @param sigma proposal distribution variance parameter
##' @export
##' @return 
##' @author chris
mcmc.blockmodel <- function(lposterior,llk_node,priors,N,P,K,niter=20,do.sm=TRUE,num.extra=0,verbose=FALSE,sigma=.1) {
  priors$sigma <- sigma
  phi <- array(0,c(P,K,K))
  phi[1,,] <- rnorm(K^2)
  z <- sample(1:K,N,rep=T)

  samples <- list()
  lps <- acc <- rep(0,niter)
  for (iter in 1:niter) {

    ## Split merge move
    if (do.sm) {
      sm <- splitmerge(phi,z,lposterior,llk_node,priors,verbose=verbose)
      phi <- sm$final$phi
      z <- sm$final$z
      acc[iter] <- sm$final$accepted
    }

    ## Sample phi
    phi <- sample_phi(phi,z,lposterior,priors)$phi

    ## Add clusters from prior (Neal 2000)
    if (num.extra > 0) {
      for (j in 1:num.extra) {
        phi <- add_cluster(phi)
        phi <- sample_cluster_from_prior(phi,dim(phi)[2],priors)
      }
    }

    ## Gibbs sample assignments
    h <- gibbs(phi,z,1:N,llk_node,N,priors)
    z <- h$z

    ## Remove empty clusters
    iz <- which(table(factor(z,1:dim(phi)[2]))==0)
    for (l in rev(iz)) {  
      z[which(z > l)] <- z[which(z > l)] - 1
      phi <- remove_cluster(phi,l)
    }
    K <- max(z)

    ## Save progress
    lps[iter] <- lposterior(phi,z,priors)
#    browser()
    cat(iter,":",lps[iter],"\n")
    if (verbose) cat(z,"\n")
    samples[[iter]] <- list(phi=phi,z=z)
  }
  return(list(lps=lps,samples=samples,acc=acc,phi=phi,z=z))
}
