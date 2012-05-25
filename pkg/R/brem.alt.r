lp <- function(phi,z,priors) {
  K <- dim(phi)[2]
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
  olp <- NULL
  if (is.null(kx)) kx <- 1:K
  for (k1 in kx) {
    for (k2 in kx) {
      cat(".")
      for (p in px) {
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
#        olp <- attr(res,"log.density")
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
