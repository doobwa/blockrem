llk_node_test <- function(a,phi,z) {
  tb <- table(edgelist[,2],edgelist[,3])
  tm <- edgelist[nrow(edgelist),1]
  llk <- 0
  for (i in 1:N) {
    for (j in (1:N)[-i]) {
      lam <- exp(phi[1,z[i],z[j]]) * tm
      llk <- llk + dpois(tb[i,j], lam, log=TRUE)
    }
  }
  return(llk)
}
llk_node <- llk_node_test


lposterior <- function(phi,z,priors) {
  K <- dim(phi)[2]
  phi[-1,,] <- 0
  llk <- llk_node_test(1,phi,z)
  pr.phi <- sum(dnorm(phi[which(phi != 0)],priors$phi$mu,priors$phi$sigma,log=TRUE))
  tb <- table(factor(z,1:K))
  tb <- tb[which(tb>0)]
  pr.z <- sum(log(sapply(tb - 1,factorial)))
  return(llk + pr.phi + pr.z)
}
lposterior_test <- lposterior

# Sample the first dimension of phi
sample_phi <- function(phi,z,lpost,kx=NULL) {
  K <- dim(phi)[2]
  if (is.null(kx)) kx <- 1:K
  lpost <- function(x,lp=TRUE,lgrad=FALSE) {
    phi[1,k1,k2] <- x
    lposterior_test(phi,z,priors)
  }
  olp <- nlp <- matrix(0,K,K)
  ratios <- c(); r <- 0
  for (k1 in kx) {
    for (k2 in kx) {
      olp[k1,k2] <- lpost(phi[1,k1,k2])
      phi[1,k1,k2] <- slice(phi[1,k1,k2],lpost,m=20)
      nlp[k1,k2] <- lpost(phi[1,k1,k2])
    }
  }
  return(list(phi=phi,olp=olp,nlp=nlp,ratios=ratios))
}


