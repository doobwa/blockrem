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

  # Do an iteration for the lowerlevel blockmodel
  fit <- mcmc.blockmodel(lp,llk_node,priors,N,P,K,do.sm=do.sm,num.extra=num.extra,niter=niter,verbose=TRUE)

  # Fit hierarchical portion
  mu <- 0
  sigma <- 1

  # TODO: Save to outfile.  Pull it out of mcmc.blockmodel.

  fit$ego <- ego
  fit$mu <- mu
  fit$sigma <- sigma
  fit$beta <- fit$samples[[niter]]$phi
  fit$z <- fit$samples[[niter]]$z
  fit$priors <- priors
  fit$s <- s
  class(fit) <- "brem"
  return(fit)
}

lp <- function(phi,z,priors) {
  K <- dim(phi)[2]
  P <- dim(phi)[1]
  if (s$get_P() != P) stop("Mismatch in parameter vector lengths")
  pr.y <- sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
  pr.phi <- sum(dnorm(unlist(phi),priors$phi$mu,priors$phi$sigma,log=TRUE))
  tb <- table(factor(z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  return(pr.y + pr.phi + pr.z)
}

llk_node <- function(a,phi,z,priors) {
  K <- dim(phi)[2]
  sum(RemLogLikelihoodActorPc(a-1,phi,z-1,s$ptr(),K))
#  sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
}

# Requires s, train, and px to be in environment
sample_phi <- function(phi,z,lpost,priors,kx=NULL) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  for (k1 in kx) {
    for (k2 in kx) {
      cat(".")
      olp <- NULL
      for (p in priors$px) {
        brem.lpost.pk1k2 <- function(x,lp=TRUE) {
          k1nodes <- which(z==k1)
          k2nodes <- which(z==k2)
          phi[p,k1,k2] <- x
          pr.y <- sum(RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,phi,z-1,s$ptr(),K))
#          pr.y <- sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
          pr.phi <- sum(dnorm(unlist(phi),priors$phi$mu,priors$phi$sigma,log=TRUE))     
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


## lp <- function(phi,z,priors) {
##   K <- dim(phi)[2]
##   pr.y <- sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
##   pr.phi <- sum(dnorm(unlist(phi),priors$phi$mu,priors$phi$sigma,log=TRUE))
##   tb <- table(factor(z,1:K))
##   tb <- tb[which(tb>0)]
##   pr.z <- sum(log(sapply(tb - 1,factorial)))
##   return(pr.y + pr.phi + pr.z)
## }

## llk_node <- function(a,phi,z,priors) {
##   K <- dim(phi)[2]
##   sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
## }

## sample_phi <- function(phi,z,lpost,priors,kx=NULL,px=1:13) {
##   K <- dim(phi)[2]
##   if (is.null(kx)) kx <- 1:K
##   for (p in px) {
##     lpost <- function(x,lp=TRUE,lgrad=FALSE) {
##       phi[p,k1,k2] <- x
##       lp(phi,z,priors)
##     }
##     for (k1 in kx) {
##       for (k2 in kx) {
##         phi[p,k1,k2] <- slice(phi[p,k1,k2],lpost,m=20)
##       }
##     }
##   }
##   return(list(phi=phi,z=z))
## }
