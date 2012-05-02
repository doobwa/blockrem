
lp <- function(phi,z,priors) {
  K <- dim(phi)[2]
  pr.y <- sum(RemLogLikelihoodPc(phi,z-1,s$ptr(),K))
  pr.phi <- sum(dnorm(unlist(phi),priors$phi$mu,priors$phi$sigma,log=TRUE))
  tb <- table(factor(z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  return(pr.y + pr.phi + pr.z)
}
lp_node <- function(a,phi,z,priors) {
  lp(phi,z,priors)
}

# Sample the first dimension of phi
sample_phi <- function(phi,z,lpost,priors,kx=NULL,px=1:13) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  for (p in px) {
    lpost <- function(x,lp=TRUE,lgrad=FALSE) {
      phi[p,k1,k2] <- x
      lp(phi,z,priors)
    }
    for (k1 in kx) {
      for (k2 in kx) {
#        olp[k1,k2] <- lpost(phi[p,k1,k2])
        phi[p,k1,k2] <- slice(phi[p,k1,k2],lpost,m=20)
#        nlp[k1,k2] <- lpost(phi[p,k1,k2])
      }
    }
  }
  return(list(phi=phi,z=z))
}
