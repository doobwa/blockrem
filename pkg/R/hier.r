
##' RemStat object s required in environment
##' @title 
##' @param params 
##' @param priors 
##' @param collapse.sigma 
##' @return 
##' @author chris
lposterior <- function(params,priors,collapse.sigma=TRUE) {
  K <- dim(params$beta)[2]
  P <- dim(params$beta)[1]
  if (s$get_P() != P) stop("Mismatch in parameter vector lengths")

  # Likelihood
  pr.y <- sum(RemLogLikelihoodPc(params$beta,params$z-1,s$ptr(),K))

  # Prior on beta
  if (collapse.sigma) {
    pr.beta <- lprior.nosigma.hier.gaussian(params$beta,params$mu,priors,grad=FALSE)
  } else {
    pr.beta <- lprior.hier.gaussian(params$beta,params$mu,params$sigma,priors,grad=FALSE)
  }

  # Prior on z
  tb <- table(factor(params$z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  
  return(pr.y + pr.beta + pr.z)
}

llk_node <- function(a,beta,z,priors) {
  K <- dim(beta)[2]
  sum(RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K))
#  sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
}

##' Requires s to be in environment
##' @title Sample block level parameters beta
##' @param beta 
##' @param z 
##' @param mu 
##' @param sigma 
##' @param priors 
##' @param kx 
##' @return 
##' @author chris
sample_beta <- function(beta,z,mu,sigma,priors,kx=NULL,collapse.sigma=TRUE) {
  K <- dim(beta)[2]
  if (is.null(kx)) kx <- 1:K
  for (k1 in kx) {
    for (k2 in kx) {
      cat(".")
      olp <- NULL
      for (p in priors$px) {
        brem.lpost.pk1k2 <- function(x) {
          k1nodes <- which(z==k1)
          k2nodes <- which(z==k2)
          beta[p,k1,k2] <- x
          pr.y <- sum(RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,beta,z-1,s$ptr(),K))
#          pr.beta <- sum(dnorm(beta[,k1,k2],mu,sigma,log=TRUE))
          if (collapse.sigma) {
            pr.beta <- lprior.nosigma.hier.gaussian(beta[,k1,k2],mu,priors,grad=FALSE)
          } else {
            pr.beta <- lprior.hier.gaussian(beta[,k1,k2],mu,sigma,priors,grad=FALSE)
          }
          return(pr.y + pr.beta)
        }
        res <- slice(beta[p,k1,k2],brem.lpost.pk1k2,m=20,olp=olp)
        olp <- attr(res,"log.density")
        beta[p,k1,k2] <- res
      }
    }
  }
  return(list(beta=beta,z=z))
}


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
  beta <- array(0,c(P,K,K))
  beta[1,,] <- rnorm(K^2)
  
  z     <- sample(1:K,N,rep=T)
  mu    <- rep(0,P)
  sigma <- rgamma(P,priors$sigma$alpha,priors$sigma$beta)

  samples <- list()
  lps <- acc <- rep(0,niter)
  for (iter in 1:niter) {

    ## Split merge move
    if (do.sm) {
      sm <- splitmerge(beta,z,lposterior,llk_node,priors,verbose=verbose)
      beta <- sm$final$phi
      z <- sm$final$z
      acc[iter] <- sm$final$accepted
    }

    ## Sample phi
    beta <- sample_beta(beta,z,mu,sigma,priors)$beta

    ## Add clusters from prior (Neal 2000)
    if (num.extra > 0) {
      for (j in 1:num.extra) {
        beta <- add_cluster(beta)
        beta <- sample_cluster_from_prior(beta,dim(beta)[2],priors)
      }
    }

    ## Gibbs sample assignments
    h <- gibbs(beta,z,1:N,llk_node,N,priors)
    z <- h$z

    ## Remove empty clusters
    iz <- which(table(factor(z,1:dim(beta)[2]))==0)
    for (l in rev(iz)) {  
      z[which(z > l)] <- z[which(z > l)] - 1
      beta <- remove_cluster(beta,l)
    }
    K <- max(z)

    ## Fit hierarchical portion
    theta <- t(array(beta,dim=c(P,K^2)))
    pr <- list(theta=theta,mu=mu,sigma=sigma)
    mu <- gibbs.mu.hier.gaussian(pr,priors)$mu
    sigma <- gibbs.sigma.hier.gaussian(pr,priors)$sigma
    mu[-priors$px] <- sigma[-priors$px] <- 0  # not included in likelihood, so set to 0 for display purposes

    ## Save progress
    params <- list(beta=beta,z=z,mu=mu,sigma=sigma)
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
