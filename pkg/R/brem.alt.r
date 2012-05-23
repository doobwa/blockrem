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

# Requires s, train to be in environment
sample_phi <- function(phi,z,lpost,priors,kx=NULL,px=1:13) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  for (p in px) {
    lpost <- function(x,lp=TRUE) {
      M <- nrow(A)
      zs <- z[train[,2]]
      zr <- z[train[,3]]
      ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
      ## Skip first and last events. TODO: Why?
      if (length(ix) > 0) {
        ix <- setdiff(ix,c(0,M-1))
      }
      
      phi[p,k1,k2] <- x
      pr.y <- sum(RemLogLikelihoodPcSubset(phi,z-1,s$ptr(),K,ix))
      if (0 %in% ix | M %in% ix) stop("out of bounds")
      pr.phi <- sum(dnorm(unlist(phi),priors$phi$mu,priors$phi$sigma,log=TRUE))     
      return(pr.y + pr.phi)
    }
    for (k1 in kx) {
      for (k2 in kx) {
        phi[p,k1,k2] <- slice(phi[p,k1,k2],lpost,m=20)
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
