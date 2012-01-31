
simulate.brem <- function(M,N,z,beta) {

  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  lambda <- beta$intercept[z,z]#matrix(beta[1],N,N)
  lrm <- array(0,c(M,N,N))
  s <- list("abba" = matrix(0,N,N),
            "abby" = matrix(0,N,N),
            "abxa" = matrix(0,N,N),
            "abxb" = matrix(0,N,N),
            "abay" = matrix(0,N,N),
            "abab" = matrix(0,N,N),
            "sod" = matrix(0,N,N),
            "rod" = matrix(0,N,N),
            "sid" = matrix(0,N,N),
            "rid" = matrix(0,N,N))
  z <- z-1
  for (m in 2:M) {
    
    i <- A[m-1,2]
    j <- A[m-1,3]
    
    s <- brem$updateStatistics(s,i-1,j-1,N)
    # Compute changes to lambda
    for (r in 1:N) {
      lambda[r,j] <- brem$computeLambda(r-1,j-1,z[r],z[j],s,beta)
      lambda[j,r] <- brem$computeLambda(j-1,r-1,z[j],z[r],s,beta)
      lambda[i,r] <- brem$computeLambda(i-1,r-1,z[i],z[r],s,beta)
      lambda[r,i] <- brem$computeLambda(r-1,i-1,z[r],z[i],s,beta)
    }
    diag(lambda) <- -Inf
    
    # Draw the next event
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    i <- cells[drawcell,1]
    j <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    A <- rbind(A,c(time,i,j))
    
    lrm[m,,] <- lambda
  }
  return(list(A=A,lrm=lrm))
}

brem.lrm <- function(A,N,z,beta) {
  M <- nrow(A)
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  z <- z - 1
  brem$lrm(beta,times,sen,rec,z,N,M)
}
brem.llk <- function(A,N,z,beta) {
  M <- nrow(A)
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  z <- z-1
  lrm <- brem$lrm(beta,times,sen,rec,z,N,M)
  brem$llk2(lrm,times,sen,rec,N,M)
}
brem.mle <- function(A,N,z,beta) {
  fn <- function(par) {
    p <- array(0,c(1,1,length(par)))
    p[1,1,] <- par
    -brem.llk(A,N,z,par)
  }
  optim(as.vector(beta),fn)$par
}
brem.lpost <- function(A,N,K,z,beta) {
  llks <- brem.llk(A,N,z,beta)
  lprior <- sum(dnorm(unlist(beta),0,1,log=TRUE)) + N * log(1/K)
  sum(llks)+lprior
}

brem.mcmc <- function(A,N,K,niter=5,model.type="full",mcmc.sd=.1,beta=NULL,z=NULL,gibbs=TRUE,mh=TRUE) {
  llks <- rep(0,niter)
  M <- nrow(A)
  P <- length(beta)
  param <- array(0,c(niter,K,K,P))
  mu <- log(M/A[M,1]/N/N)  # MLE for beta for uniform rate across all events
  sigma <- .5
  current <- lapply(1:P,function(i) matrix(rnorm(K^2,mu,sigma),K,K))
  if (is.null(z))    z <- sample(1:K,N,replace=TRUE)
  if (!is.null(beta)) current <- beta
  
  for (iter in 1:niter) {
    
    # For each effect sample via MH
    first.iter <- (iter==1)
    current <- brem.mh(A,N,K,P,z,current,model.type,mcmc.sd,first.iter)
    
    # Sample empty cluster parameters from prior
#     for (k in 1:K) {
#       if (length(which(z==k))==0) {
#         current[k,,] <- rnorm(K*P,0,5)
#         current[,k,] <- rnorm(K*P,0,5)
#       }
#     }
    
    if (gibbs) {
      # Gibbs sample assignments
      for (i in 1:N) {
        ps <- rep(0,K)
        for (k in 1:K) {
          z[i] <- k
          ps[k] <- brem.lpost(A,N,K,z,current)
        }
        ps <- exp(ps - max(ps))
        z[i] <- sample(1:K,size=1,prob=ps)
      }
    }
    
    param[iter,,,] <- current
    llks[iter] <- brem.lpost(A,N,K,z,current)
    
    cat("iter",iter,":",llks[iter],"z:",table(z),"\n")
    
  }
  return(list(z=z,llks=llks,param=param,beta=current))
}
brem.mh <- function(A,N,K,P,z,current,model.type="baserates",mcmc.sd=.1,first.iter=FALSE) {
  if (first.iter) olp <- -Inf
  else olp <- brem.lpost(A,N,K,z,current)
  cand <- current
  if (model.type=="baserates") {
    cand[["intercept"]] <- cand[["intercept"]] + rnorm(K^2,0,mcmc.sd)
    #cand[,,-1] <- 0
    clp <- brem.lpost(A,N,K,z,cand)
    if (clp - olp > log(runif(1))) {
      current <- cand
      olp <- clp
    }
  }
#   if (model.type=="diag.rem") {
#     a <- cand[1,1,] + rnorm(P,0,mcmc.sd)
#     b <- cand[1,2,] + rnorm(P,0,mcmc.sd)
#     for (k1 in 1:K) {
#       for (k2 in 1:K) {
#         if (k1==k2) {
#           cand[k1,k2,] <- a
#         } else {
#           cand[k1,k2,] <- b
#         }
#       }
#     }
#     clp <- brem.lpost(A,N,K,z,cand)
#     if (clp - olp > log(runif(1))) {
#       current <- cand
#       olp <- clp
#     }
#   }
  if (model.type=="full") {
    for (p in 1:length(beta)) {
      cand <- current
      cand[[p]]  <- cand[[p]] + rnorm(K^2,0,mcmc.sd)
      clp <- brem.lpost(A,N,K,z,cand)
      if (clp - olp > log(runif(1))) {
        current <- cand
        olp <- clp
      }
    }
  }
  return(current)
}