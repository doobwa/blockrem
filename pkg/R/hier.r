brem <- function(train,N,K=2,effects=c("intercept","abba","abby","abay"),ego=TRUE,do.sm=FALSE,num.extra=2,niter=20,verbose=TRUE) {
  M <- nrow(train)
  P <- 13
  ego <- ego*1  # RemStat doesn't want boolean
  s <- new(RemStat,
           train[,1],
           as.integer(train[,2])-1,
           as.integer(train[,3])-1,
           N,M,P,ego)
  s$precompute()
  priors$px <- which(effects %in%
                     c("intercept","abba","abby","abxa","abxb","abay","abab",
                       "sen_outdeg","rec_outdeg","sen_indeg","rec_indeg",
                       "dyad_count","changepoint_count"))

  # Initialize parameters
  phi <- array(0,c(P,K,K))
  phi[1,,] <- rnorm(K^2)
  
  z     <- sample(1:K,N,rep=T)
  mu    <- rep(0,P)
  sigma <- rgamma(P,priors$sigma$alpha,priors$sigma$beta)

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
    phi <- sample_phi(phi,z,mu,sigma,priors)$phi

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

    ## Fit hierarchical portion
    ## TODO: Sample these
    mu <- 0
    sigma <- 1

    ## Save progress
    params <- list(phi=phi,z=z,mu=mu,sigma=sigma)
    lps[iter] <- lposterior(params,priors)
    cat(iter,":",lps[iter],"\n")
    if (verbose) cat(z,"\n")
    samples[[iter]] <- params
#    browser()

  }

  # TODO: Save to outfile.  Pull it out of mcmc.blockmodel.

  fit <- list(params=params,samples=samples,ego=ego,priors=priors,s=s)
  class(fit) <- "brem"
  return(fit)
}

##' RemStat object s required in environment
##' @title 
##' @param params 
##' @param priors 
##' @param collapse.sigma 
##' @return 
##' @author chris
lposterior <- function(params,priors,collapse.sigma=TRUE) {
  K <- dim(params$phi)[2]
  P <- dim(params$phi)[1]
  if (s$get_P() != P) stop("Mismatch in parameter vector lengths")

  # Likelihood
  pr.y <- sum(RemLogLikelihoodPc(params$phi,params$z-1,s$ptr(),K))

  # Prior on phi
  if (collapse.sigma) {
    pr.phi <- lprior.nosigma.hier.gaussian(params$phi,params$mu,priors,grad=FALSE)
  } else {
    pr.phi <- lprior.hier.gaussian(params$phi,params$mu,params$sigma,priors,grad=FALSE)
  }

  # Prior on z
  tb <- table(factor(params$z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  
  return(pr.y + pr.phi + pr.z)
}

llk_node <- function(a,phi,z,priors) {
  K <- dim(phi)[2]
  sum(RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K))
#  sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
}

##' Requires s to be in environment
##' @title Sample block level parameters phi
##' @param phi 
##' @param z 
##' @param mu 
##' @param sigma 
##' @param priors 
##' @param kx 
##' @return 
##' @author chris
sample_phi <- function(phi,z,mu,sigma,priors,kx=NULL) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  for (k1 in kx) {
    for (k2 in kx) {
      cat(".")
      olp <- NULL
      for (p in priors$px) {
        brem.lpost.pk1k2 <- function(x) {
          k1nodes <- which(z==k1)
          k2nodes <- which(z==k2)
          phi[p,k1,k2] <- x
          pr.y <- sum(RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,phi,z-1,s$ptr(),K))
          pr.phi <- sum(dnorm(phi[,k1,k2],mu,sigma,log=TRUE)) 
          return(pr.y + pr.phi)
        }
        res <- slice(phi[p,k1,k2],brem.lpost.pk1k2,m=20,olp=olp)
        olp <- attr(res,"log.density")
        phi[p,k1,k2] <- res
      }
    }
  }
  return(list(phi=phi,z=z))
}
