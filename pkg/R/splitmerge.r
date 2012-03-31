

##' Split phi for cluster k into k and l
split_phi <- function(phi.merge,k,l,sigma=.1) {
  P <- dim(phi.merge)[1]
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
  return(phi[,-l,-l])
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


gibbs_restricted <- function(phi,z,S,k,l,llk_node,ztrue=NULL) {

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
    z[a] <- sample(clusters,1,prob=lnormalize(ys[a,clusters])) #rcatlog(ys[a,])
    if (!is.null(ztrue)) z[a] <- ztrue[a]
    counts[z[a]] <- counts[z[a]] + 1
  }
#  browser()
  ys[S,clusters] <- t(apply(ys[S,clusters],1,lnormalize))
  q <- sum(log(ys[cbind(S,z[S])]))
  list(z=z,transition=q,ys=ys,z.init=z.init)
}

ps2pm <- function(phi.split,phi.merge) {
  P <- dim(phi.merge)[1]
  K <- dim(phi.merge)[2]
  log(.5^(P*(2*K-3)))
}
pm2ps <- function(phi.merge,phi.split,k,l,sigma) {
  K <- dim(phi.split)[2]
  rs <- (1:K)[-c(k,l)]
  ps <- lapply(rs,function(r) {
    c(dnorm(phi.split[,l,r],phi.merge[,k,r],sigma,log=TRUE),
      dnorm(phi.split[,r,l],phi.merge[,r,k],sigma,log=TRUE))
  })
  ps <- c(ps,list(dnorm(phi.split[,l,l],phi.merge[,k,k],sigma,log=TRUE)))
  ps <- unlist(ps)
  return(sum(ps))
}

splitmerge <- function(phi,z,lposterior,llk_node,priors,verbose=TRUE) {
  N <- length(z)
  K <- dim(phi)[2]
  ij <- sample(1:N,2,replace=FALSE)
  i <- ij[1]; j <- ij[2]
  S <- which(z %in% z[ij])
  z <- list(init=z)
  phi <- list(init=phi)
  q <- list()
  if (z$init[i] == z$init[j]) {
    phi$merge <- add_cluster(phi$init)
    phi$merge <- sample_cluster_from_prior(phi$merge,K+1,priors)
    z$merge <- z$init
    phi$split <- split_phi(phi$merge,z$init[i],K+1,priors$sigma)
    h <- gibbs_restricted(phi$split,z$merge,S,z$init[i],K+1,llk_node)
    z$split <- h$z
    q$zm2zs <- h$transition
    q$pm2ps <- pm2ps(phi$merge,phi$split,z$init[i],K+1,priors$sigma)
  } else {
    phi$split <- phi$init
    z$split <- z$merge <- z$init
    z$merge[S] <- z$init[i]
    phi$merge <- merge_phi(phi$split,z$init[i],z$init[j])
    #phi$merge <- sample_cluster_from_prior(phi$merge,z$init[j],priors)
    #phi$split <- sample_cluster_from_prior(phi$split,z$init[j],priors)
    g <- gibbs_restricted(phi$split,z$merge,S,z$init[i],z$init[j],llk_node,ztrue=z$split)
    q$zm2zs <- g$transition
    if (q$zm2zs == -Inf) browser()
    q$pm2ps <- pm2ps(phi$merge,phi$split,z$init[i],z$init[j],priors$sigma)
  }
  if (any(is.nan(unlist(q)))) browser()
  q$ps2pm <- ps2pm(phi$split,phi$merge)
  q$zs2zm <- 0
  lp.split <- lposterior(phi$split,z$split,priors)
  lp.merge <- lposterior(phi$merge,z$merge,priors)
  alpha <- lp.split - lp.merge + q$zs2zm - q$zm2zs + q$ps2pm - q$pm2ps
  
  if (z$init[i] == z$init[j]) {
    final <- list(type="split",phi=phi$split,z=z$split,prob=exp(alpha))
  } else {
    final <- list(type="merge",phi=phi$merge,z=z$merge,prob=1/exp(alpha))
  }
  if (verbose) {
    type <- ifelse(z$init[i] == z$init[j],"split","merge") 
    cat("move:",type," clusters:",z$init[ij],"\n")
    cat("initial z:",z$init,", K:",dim(phi$init)[2],"\n")
    cat("merge z:",z$merge,", K:",dim(phi$merge)[2],"\n")
    cat("split z:",z$split,", K:",dim(phi$split)[2],"\n")
    cat("prob:",final$prob,"\n")
    print(lapply(phi[c("init","merge","split")],function(p) p[1,,]))
    print(unlist(q))
  }
  list(phi=phi,chosen=c(i,j),alpha=alpha,z=z,q=q,final=final,ys=g$ys)
}
