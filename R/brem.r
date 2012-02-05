
simulate.brem <- function(M,N,z,beta) {

  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  lambda <- beta[1,z,z]
  lrm <- array(0,c(M,N,N))
  s <- brem$initializeStatistics(N,P);
  z <- z-1
  for (m in 2:M) {
    
    i <- A[m-1,2]
    j <- A[m-1,3]
    
    s <- brem$updateStatistics(s,i-1,j-1,N,P)
    # Compute changes to lambda
    for (r in 1:N) {
      lambda[r,j] <- brem$computeLambda(r-1,j-1,z[r],z[j],s,beta,N,K,P)
      lambda[j,r] <- brem$computeLambda(j-1,r-1,z[j],z[r],s,beta,N,K,P)
      lambda[i,r] <- brem$computeLambda(i-1,r-1,z[i],z[r],s,beta,N,K,P)
      lambda[r,i] <- brem$computeLambda(r-1,i-1,z[r],z[i],s,beta,N,K,P)
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
  P <- dim(beta)[1]
  K <- dim(beta)[2]
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  z <- z - 1
  brem$lrm(beta,times,sen,rec,z,N,M,K,P)
}
brem.llk <- function(A,N,z,beta,use.lrm=FALSE) {
  M <- nrow(A)
  P <- dim(beta)[1]
  K <- dim(beta)[2]
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  z <- z-1
  if (use.lrm) {
    lrm <- brem$lrm(beta,times,sen,rec,z,N,M,K,P)
    return(brem$llk2(lrm,times,sen,rec,N,M))
  } else {
    return(brem$llk(beta,times,sen,rec,z,N,M,K,P))
  }
}
brem.mle <- function(A,N,K,P,z,beta=NULL) {
  fn <- function(par) {
    beta <- array(par,c(P,K,K))
    beta[7:11,,] <- 0
    -brem.llk(A,N,z,beta)
  }
  if (is.null(beta)) beta <- array(0,c(P,K,K))
  beta.hat <- optim(as.vector(beta),fn)$par
  array(beta.hat,c(P,K,K))
}
brem.lpost <- function(A,N,K,z,beta) {
  llks <- brem.llk(A,N,z,beta)
  lprior <- sum(dnorm(unlist(beta),0,1,log=TRUE)) + N * log(1/K)
  sum(llks)+lprior
}
brem.lpost.fast <- function(A,N,K,z,sptr,beta) {
  llks <- brem$llkfast(beta,z-1,sptr,K)
  lprior <- sum(dnorm(unlist(beta),0,1,log=TRUE)) + N * log(1/K)
  sum(llks)+lprior
}

brem.mcmc <- function(A,N,K,s,niter=5,model.type="full",mcmc.sd=.1,beta=NULL,z=NULL,gibbs=TRUE,mh=TRUE) {
  llks <- rep(0,niter)
  M <- nrow(A)
  P <- 11
  param <- array(0,c(niter,K,K,P))
  px <- c(1,1,1,1,1,1,0,0,0,0)  # skip degree effects
  mu <- 0#log(M/A[M,1]/N/N)  # MLE for beta for uniform rate across all events
  sigma <- .5
  current <- array(rnorm(P*K^2,mu,sigma),c(P,K,K))
  current[7:11,,] <- 0
  if (is.null(z))    z <- sample(1:K,N,replace=TRUE)
  if (!is.null(beta)) current <- beta
  
  for (iter in 1:niter) {
    
    # For each effect sample via MH
    first.iter <- (iter==1)
    current <- brem.mh(A,N,K,P,z,current,model.type,mcmc.sd,first.iter,px)
    
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
          ps[k] <- brem.lpost.fast(A,N,K,z,s,current)
        }
        ps <- exp(ps - max(ps))
        z[i] <- sample(1:K,size=1,prob=ps)
      }
    }
    
    param[iter,,,] <- current
    llks[iter] <- brem.lpost.fast(A,N,K,z,s,current)
    
    cat("iter",iter,":",llks[iter],"z:",table(z),"\n")
    
  }
  return(list(z=z,llks=llks,param=param,beta=current))
}
brem.mh <- function(A,N,K,P,z,current,model.type="baserates",mcmc.sd=.1,first.iter=FALSE,px=c(1,1,1,1,1,1,1,0,0,0,0)) {
  if (first.iter) {
    olp <- -Inf
  }  else {
    olp <- brem.lpost.fast(A,N,K,z,s,current)
  }
  
  cand <- current
  if (model.type=="baserates") {
    cand[1,,] <- cand[1,,] + rnorm(K^2,0,mcmc.sd)
    cand[-1,,] <- 0
    clp <- brem.lpost.fast(A,N,K,z,s,cand)
    if (clp - olp > log(runif(1))) {
      current <- cand
      olp <- clp
    }
  }
  if (model.type=="diag.rem") {
    for (p in which(px==1)) {
      
      # Sample new values for diagonal element (a) and all off diagonal (b)
      cand <- current
      a <- cand[p,1,1] + rnorm(1,0,mcmc.sd)  
      b <- cand[p,1,2] + rnorm(1,0,mcmc.sd) 
      
      # Place these values in the parameter array
      for (k1 in 1:K) {
        for (k2 in 1:K) {
          if (k1==k2) {
            cand[p,k1,k2] <- a
          } else {
            cand[p,k1,k2] <- b
          }
        }
      }
      
      # MH sampler
      clp <- brem.lpost.fast(A,N,K,z,s,cand)
      if (clp - olp > log(runif(1))) {
        current <- cand
        olp <- clp
      }
    }
  }
  if (model.type=="full") {
    cand <- current
    for (p in which(px==1)) {
      cand[p,,]  <- cand[p,,] + rnorm(K^2,0,mcmc.sd)
      clp <- brem.lpost.fast(A,N,K,z,s,cand)
      if (clp - olp > log(runif(1))) {
        current <- cand
        olp <- clp
      }
    }
  }
  return(current)
}
brem.llk.slow <- function(lrm,times,sen,rec,z,N,M) {
  mp <- matrix(1,N,N)
  llk <- 0
  llks <- rep(0,4)
  for (m in 2:(M-1)) {
    i = sen[m];
    j = rec[m];
    zi = z[i];
    zj = z[j];
    llk = llk + lrm[m,i,j]
    
    for (r in 1:N) {
      zr = z[r];
      if (r != i) {
        lam  = lrm[m,i,r]
        llk  = llk - (times[m] - times[mp[i,r]]) * exp(lam);
        lam  = lrm[m,r,i]
        llk  = llk - (times[m] - times[mp[r,i]]) * exp(lam);
        mp[i,r] = m;
        mp[r,i] = m;
      }
      if (r != j) {
        lam  = lrm[m,j,r]
        llk  = llk - (times[m] - times[mp[j,r]]) * exp(lam);
        lam  = lrm[m,r,j]
        llk  = llk - (times[m] - times[mp[r,j]]) * exp(lam);
        mp[j,r] = m;
        mp[r,j] = m;
      }
    }
    llks[m] <- llk
  }
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        lam  = lrm[M,i,j]
        llk  = llk - (times[M] - times[mp[i,j]]) * exp(lam);
      }
    }
  }
  llks[M] <- llk
  return(llks)
}