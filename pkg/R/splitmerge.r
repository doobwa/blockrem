

##' Split phi for cluster k into k and l
split_phi <- function(phi.merge,k,l,sigma=.1) {
  P <- dim(phi.merge)[1]
  K <- dim(phi.merge)[2]
  phi.split <- phi.merge
  for (r in 1:K) {
    phi.split[,l,r] <- rnorm(P,phi.merge[,k,r],sigma)
    phi.split[,r,l] <- rnorm(P,phi.merge[,r,k],sigma)
  }
  phi.split[,l,l] <- rnorm(P,phi.merge[,k,k],sigma)
  return(phi.split)
}

##' Merge cluster l into k
merge_phi <- function(phi,k,l) {
  P <- dim(phi)[1]
  K <- dim(phi)[2]
  rs <- (1:K)[-c(k,l)]
  for (p in 1:P) {
    for (r in rs) {
      phi[p,r,k] <- sample(phi[p,r,c(k,l)],1)
      phi[p,k,r] <- sample(phi[p,c(k,l),r],1)
    }
    phi[p,k,k] <- sample(c(phi[p,k,k],phi[p,l,l]),1)
  }
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

##' @return new z and the transition probabilities

lnormalize <- function(ps) {
  ps <- ps - max(ps)
  ps <- exp(ps)
  return(ps/sum(ps))
}


gibbs <- function(phi,z,S,llk_node,prob.only=FALSE) {

  # Randomly assign nodes that are currently either in k or l
  K <- dim(phi)[2]
  z.init <- z
  counts <- table(factor(z,1:K))
  ys <- matrix(0,N,K)
  for (a in S) {
    counts[z[a]] <- counts[z[a]] - 1
    for (j in 1:K) {
      z[a] <- j
      ys[a,j] <- llk_node(a,phi,z) + log(counts[j] + priors$alpha)
    }
    z[a] <- sample(1:K,1,prob=lnormalize(ys[a,])) #rcatlog(ys[a,])
    counts[z[a]] <- counts[z[a]] + 1
  }
  ys[S,] <- t(apply(ys[S,],1,lnormalize))
  q <- sum(log(ys[cbind(S,z[S])]))
  list(z=z,transition=q,ys=ys,z.init=z.init)
}

gibbs_restricted <- function(phi,z,S,k,l,llk_node,prob.only=FALSE) {

  # Randomly assign nodes that are currently either in k or l
  K <- dim(phi)[2]
  clusters <- c(k,l)
#  S <- which(z %in% clusters)
  z[S] <- sample(clusters,length(S),replace=TRUE)
  z.init <- z
  counts <- table(factor(z,1:K))
  ys <- matrix(0,N,K)
  for (a in S) {
    counts[z[a]] <- counts[z[a]] - 1
    for (j in clusters) {
      z[a] <- j
      ys[a,j] <- llk_node(a,phi,z) + log(counts[j] + priors$alpha)
    }
    if (!prob.only) {
      z[a] <- sample(clusters,1,prob=lnormalize(ys[a,clusters])) #rcatlog(ys[a,])
    }
    counts[z[a]] <- counts[z[a]] + 1
  }
#  browser()
  ys[S,clusters] <- t(apply(ys[S,clusters],1,lnormalize))
  q <- sum(log(ys[cbind(S,z[S])]))
  list(z=z,transition=q,ys=ys,z.init=z.init)
}

## ps2pm <- function(phi.split,phi.merge) {
##   P <- 1#dim(phi.merge)[1]
##   K <- dim(phi.merge)[2]
##   log(.5^(P*(2*K-3)))
## }
## pm2ps <- function(phi.merge,phi.split,k,l,sigma) {
##   P <- 1
##   K <- dim(phi.split)[2]
##   rs <- (1:K)[-c(k,l)]
##   ps <- lapply(rs,function(r) {
##     c(dnorm(phi.split[1:P,l,r],phi.merge[1:P,k,r],sigma,log=TRUE),
##       dnorm(phi.split[1:P,r,l],phi.merge[1:P,r,k],sigma,log=TRUE))
##   })
##   ps <- c(ps,list(dnorm(phi.split[1:P,l,l],phi.merge[1:P,k,k],sigma,log=TRUE)))
##   ps <- unlist(ps)
##   return(sum(ps))
## }


ps2pm <- function(phi.split,phi.merge) {
  P <- 1#dim(phi.merge)[1]
  K <- dim(phi.merge)[2]
  log(.5^(P*(2*K-3)))
}
pm2ps <- function(phi.merge,phi.split,k,l,sigma) {
  P <- 1
  K <- dim(phi.split)[2]
  rs <- (1:K)[-c(k,l)]
  ps <- lapply(rs,function(r) {
    c(dnorm(phi.split[1:P,l,r],phi.merge[1:P,k,r],sigma,log=TRUE),
      dnorm(phi.split[1:P,r,l],phi.merge[1:P,r,k],sigma,log=TRUE))
  })
  ps <- c(ps,list(dnorm(phi.split[1:P,l,l],phi.merge[1:P,k,k],sigma,log=TRUE)))
  ps <- unlist(ps)
  return(sum(ps))
}

splitmerge <- function(phi,z,lposterior,llk_node,priors,verbose=TRUE) {
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
    K <- max(split$z)
    split$z[S] <- sample(c(k,l),length(S),replace=TRUE)
    split$phi <- split_phi(split$phi,k,l)
    for (iter in 1:3) {
      split$phi <- sample_phi(split$phi,split$z,lposterior,kx=1:K)$phi #fix kx
      h <- gibbs_restricted(split$phi,split$z,S,k,l,llk_node,prob.only=FALSE)
      split$z <- h$z
    }
  } else {
    merge$z[S] <- k
    merge$phi <- remove_cluster(merge$phi,l)
    merge$z[which(merge$z > l)] <- merge$z[which(merge$z > l)] - 1
    K <- max(merge$z)
    for (iter in 1:3) {
      merge$phi <- sample_phi(merge$phi,merge$z,lposterior,kx=1:K)$phi #fix kx
    }
  }
  lp.split <- lposterior(split$phi,split$z,priors)
  lp.merge <- lposterior(merge$phi,merge$z,priors)
  zm2zs <- .5*length(S)
  pzs2pzm <- factorial(sum(split$z[S] == k) ) *
    factorial(sum(split$z[S] == l) ) / factorial(length(S) )
  alpha <- lp.split - lp.merge - log(zm2zs) + log(pzs2pzm)
  
  if (any(is.nan(unlist(q)))) browser()

  initial <- list(type="nomove",phi=init$phi,z=init$z)
  if (init$z[i] == init$z[j]) {
    proposed <- list(type="split",phi=split$phi,z=split$z,prob=exp(alpha))
  } else {
    proposed <- list(type="merge",phi=merge$phi,z=merge$z,prob=1/exp(alpha))
  }
  
  if (runif(1) < proposed$prob) {
    final <- proposed
    cat("move accepted\n")
  } else {
    final <- initial
    cat("move rejected\n")
  }

  if (verbose) {
    type <- ifelse(init$z[i] == init$z[j],"split","merge") 
    cat("move:",type," clusters:",init$z[ij],"\n")
    cat("initial z:",init$z,", K:",dim(init$phi)[2],"\n")
    cat("merge z:",merge$z,", K:",dim(merge$phi)[2],"\n")
    cat("split z:",split$z,", K:",dim(split$phi)[2],"\n")
    cat("prob:",final$prob,"\n")
    #print(lapply(phi[c("init","merge","split")],function(p) p[1,,]))
    cat("split",lp.split,"\n")
    cat("merge",lp.merge,"\n")
    #print(unlist(q))
  }
  list(phi=phi,chosen=c(i,j),alpha=alpha,z=z,q=q,final=final)
}


splitmerge_old <- function(phi,z,lposterior,llk_node,priors,verbose=TRUE) {
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
    l <- K+1
    split$phi <- add_cluster(split$phi)
    split$phi <- split_phi(split$phi,k,l,priors$sigma)
    h <- gibbs_restricted(split$phi,merge$z,S,k,l,llk_node,prob.only=FALSE)
  } else {
    merge$z[S] <- k
    merge$phi <- merge_phi(split$phi,k,l)
    h <- gibbs_restricted(split$phi,merge$z,S,k,l,llk_node,prob.only=TRUE)
  }
  split$z <- h$z
  q$zm2zs <- h$transition
  q$pm2ps <- pm2ps(merge$phi,split$phi,k,l,priors$sigma)
  q$ps2pm <- ps2pm(split$phi,merge$phi)
  q$zs2zm <- 0
  lp.split <- lposterior(split$phi,split$z,priors)
  lp.merge <- lposterior(merge$phi,merge$z,priors)
  alpha <- lp.split - lp.merge + q$zs2zm - q$zm2zs + q$ps2pm - q$pm2ps
  
  if (any(is.nan(unlist(q)))) browser()
  if (init$z[i] == init$z[j]) {
    final <- list(type="split",phi=split$phi,z=split$z,prob=exp(alpha))
  } else {
    final <- list(type="merge",phi=merge$phi,z=merge$z,prob=1/exp(alpha))
  }
  if (verbose) {
    type <- ifelse(init$z[i] == init$z[j],"split","merge") 
    cat("move:",type," clusters:",init$z[ij],"\n")
    cat("initial z:",init$z,", K:",dim(init$phi)[2],"\n")
    cat("merge z:",merge$z,", K:",dim(merge$phi)[2],"\n")
    cat("split z:",split$z,", K:",dim(split$phi)[2],"\n")
    cat("prob:",final$prob,"\n")
    #print(lapply(phi[c("init","merge","split")],function(p) p[1,,]))
    cat("split",lp.split,"\n")
    cat("merge",lp.merge,"\n")
    print(unlist(q))
  }
  list(phi=phi,chosen=c(i,j),alpha=alpha,z=z,q=q,final=final,ys=h$ys)
}
