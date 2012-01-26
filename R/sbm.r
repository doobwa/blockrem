sbm.llk <- function(A,N,K,z,beta) {
#   s <- factor(A[,2],1:N)
#   r <- factor(A[,3],1:N)
  zs <- factor(z[A[,2]],1:K)
  zr <- factor(z[A[,3]],1:K)
  mat <- table(zs,zr)
  t_M <- A[nrow(A),1]
  R <- table(factor(z,1:K))
  R <- outer(R,R)
  mat * log(beta) - t_M * R * beta
  #dmultinom(mat,size=sum(mat),prob=eta,log=TRUE)  
}
sbm.lpost <- function(A,N,K,z,beta,alpha=1) {
  beta <- exp(beta)
  sum(sbm.llk(A,N,K,z,beta)) + sum(dgamma(beta,1,1,log=TRUE))
}
sbm.lrm <- function(A,N,z,beta) {
  a <- array(0,c(nrow(A),N,N))
  for (i in 1:nrow(A)) {
    a[i,,] <- beta[z,z]
  }
  return(a)
}
sbm.mcmc <- function(A,N,K,niter=100,z=NULL,mcmc.sd=.1) {
  if (is.null(z)) z <- sample(1:K,N,replace=TRUE)
  beta <- matrix(rnorm(K^2),K,K)
  llks <- rep(0,niter)
  current <- beta
  for (iter in 1:niter) {
    # For each effect sample via MH
    olp <- sbm.lpost(A,N,K,z,current)
    for (p in 1:5) {
      cand <- current
      cand  <- cand + rnorm(K^2,0,mcmc.sd)
      clp <- sbm.lpost(A,N,K,z,cand)
      if (clp - olp > log(runif(1))) {
        current <- cand
        olp <- clp
      }
    }
    
#     for (i in 1:N) {
#       ps <- rep(0,K)
#       for (k in 1:K) {
#         z[i] <- k
#         ps[k] <- sbm.lpost(A,N,K,z,current)
#       }
#       ps <- exp(ps - max(ps))
#       z[i] <- sample(1:K,size=1,prob=ps)
#     }

    llks[iter] <- olp
    cat("iter",iter,":",llks[iter],"z:",table(z),"\n")
  }
  return(list(beta=current,z=z))
}
mult.dir <- function(dataset,dims,prior=1/prod(dims)) {
  if (ncol(dataset)==2)
    e <- paste(dataset[,1],dataset[,2])
  if (ncol(dataset)==3)
    e <- paste(dataset[,1],dataset[,2],dataset[,3])
  tb <- table(e)
  prob <- (tb + prior)/(sum(tb + prior) + prod(dims)*prior)
  b <- strsplit(names(tb)," ")
  b <- lapply(b,as.numeric)
  b <- do.call(rbind,b)
  b <- cbind(b,tb,prob)
  rownames(b) <- c()
  colnames(b) <- c(paste("d",1:length(dims),sep=""),"count","prob")
  b <- as.data.frame(b)
  return(b)
}