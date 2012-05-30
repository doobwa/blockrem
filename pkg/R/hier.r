lpz <- function(count,priors) {
  if (!is.null(priors$nb)) {
    dnbinom(count,priors$nb$shape,mu=priors$nb$mu,log=TRUE)
  } else {
    log(count + priors$alpha)
  }
}

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
  y <- sum(RemLogLikelihoodPc(params$beta,params$z-1,s$ptr(),K))

  # Prior on beta
  if (collapse.sigma) {
    beta <- lprior.nosigma.hier.gaussian(params$beta,params$mu,priors,grad=FALSE)
  } else {
    beta <- lprior.hier.gaussian(params$beta,params$mu,params$sigma,priors,grad=FALSE)
  }

  # Prior on z
  tb <- table(factor(params$z,1:K))
  z <- sapply(tb,function(count) count * lpz(count,priors))

  # Prior on sigma
  sigma <- dgamma(params$sigma,priors$sigma$alpha,priors$sigma$beta,log=TRUE)

  r <- list(y=y,beta=beta,z=z,sigma=sigma)
  r$all <- sum(r$sigma[priors$px]) + y + beta + sum(z)
  return(r)
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param train edgelist
##' @param N maximum id of actors
##' @param priors alpha, sigma.proposal, phi=list(mu=0,sigma=1), mu=list(mu=0,sigma=1), sigma=list(alpha=2,beta=2))
##' @param K initial number of clusters
##' @param effects character vector that is a subset of "intercept","abba","abby","abxa","abxb","abay","abab","sen_outdeg","rec_outdeg","sen_indeg","rec_indeg","dyad_count","changepoint_count"
##' @param ego 
##' @param transform 
##' @param do.sm perform split merge moves
##' @param num.extra number of extra clusters to propose from the prior each iteration
##' @param niter number of iterations
##' @param collapse.sigma 
##' @param verbose 
##' @param outfile save progress to this file
##' @return 
##' @author chris
brem <- function(train,N,priors,K=2,effects=c("intercept","abba","abby","abay"),ego=TRUE,transform=FALSE,do.sm=FALSE,num.extra=2,niter=20,collapse.sigma=TRUE,verbose=TRUE,outfile=NULL) {
  M <- nrow(train)
  P <- 13
  ego <- ego*1  # RemStat doesn't want boolean
  s <<- new(RemStat,
           train[,1],
           as.integer(train[,2])-1,
           as.integer(train[,3])-1,
           N,M,P,ego)
  s$precompute()
  if (transform) s$transform()
  # TODO: Transform

  enam <- c("intercept","abba","abby","abxa","abxb","abay","abab",
            "sen_outdeg","rec_outdeg","sen_indeg","rec_indeg",
            "dyad_count","changepoint_count")
  priors$px <- match(effects,enam)

  # Initialize parameters
  beta <- array(0,c(P,K,K))
  beta[1,,] <- rnorm(K^2)
  
  z     <- sample(1:K,N,rep=T)
  mu    <- rep(0,P)
  sigma <- rgamma(P,priors$sigma$alpha,priors$sigma$beta)
  beta  <- sample_beta(beta,z,mu,sigma,priors,collapse.sigma=collapse.sigma)$beta

  samples <- list()
  lps <- llks <- acc <- rep(0,niter)
  for (iter in 1:niter) {

    ## Add clusters from prior (Neal 2000)
    if (num.extra > 0 & K!=1) {
      prs <- priors
      prs$phi <- list(mu=mu,sigma=sigma)
      for (j in 1:num.extra) {
        beta <- add_cluster(beta)
        beta <- sample_cluster_from_prior(beta,dim(beta)[2],prs)
      }
    }

    ## Gibbs sample assignments
    h <- gibbs(beta,z,1:N,llk_node,lpz,N,priors)
    z <- h$z

    ## Remove empty clusters
    iz <- which(table(factor(z,1:dim(beta)[2]))==0)
    for (l in rev(iz)) {  
      z[which(z > l)] <- z[which(z > l)] - 1
      beta <- remove_cluster(beta,l)
    }
    K <- max(z)

    ## Split merge move
    if (do.sm & K != 1) {
      lpost <- function(phi,z,priors) {
        lposterior(list(beta=phi,z=z,mu=mu,sigma=sigma),priors)$all
      }
      prs <- priors
      prs$phi <- list(mu=mu,sigma=sigma)
      prs$tau <- sigma/2
      sm <- splitmerge(beta,z,lpost,llk_node,lpz,prs,verbose=verbose)
      beta <- sm$final$phi
      z <- sm$final$z
      acc[iter] <- sm$final$accepted
    }

    ## Sample phi
    beta <- sample_beta(beta,z,mu,sigma,priors)$beta
    
    ## Fit hierarchical portion
    if (K != 1) {
      theta <- t(array(beta,dim=c(P,K^2)))
      pr <- list(theta=theta,mu=mu,sigma=sigma)
      mu <- gibbs.mu.hier.gaussian(pr,priors)$mu
      sigma <- gibbs.sigma.hier.gaussian(pr,priors)$sigma
    }
    
    ## Save progress
    # not included in likelihood, so set to 0 for display purposes
    mu[-priors$px] <- sigma[-priors$px] <- 0
    params <- list(beta=beta,z=z,mu=mu,sigma=sigma,ego=ego)
    lp <- lposterior(params,priors)
    lps[iter] <- lp$all
    llks[iter] <- lp$y
    cat("\n",iter,": llk",llks[iter]," lp",lps[iter],"\n")
    if (verbose) cat(z,"\n")
    samples[[iter]] <- params
    fit <- list(params=params,samples=samples,ego=ego,transform=transform,priors=priors,lps=lps,llks=llks,effects=effects,lp=lp,iter=iter,niter=niter)
    class(fit) <- "brem"
    if (!is.null(outfile)) save(fit,file=outfile)
  }

  return(fit)
}
