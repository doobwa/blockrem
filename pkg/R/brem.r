#' Simulate data from the block relational event model
#' @M number of events
#' @N number of actors
#' @z latent class assignments
#' @beta P x K x K array of parameters
simulate.brem <- function(M,N,z,beta) {
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
  lrm <- log_intensity_array(beta,times,sen,rec,z,N,M,K,P)
  for (i in 1:M) diag(lrm[i,,]) <- -Inf
  return(lrm)
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
    tmp <- log_intensity_array(beta,times,sen,rec,z,N,M,K,P)
    return(loglikelihood_from_lrm(tmp,times,sen,rec,N,M))
  } else {
    return(loglikelihood(beta,times,sen,rec,z,N,M,K,P))
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
brem.lpost <- function(A,N,K,z,beta,priors) {
  llks <- brem.llk(A,N,z,beta)
  lprior <- sum(dnorm(unlist(beta),priors$beta$mu,priors$beta$sigma,log=TRUE)) + N * log(1/K)
  sum(llks)+lprior
}
brem.lpost.fast <- function(A,N,K,z,s,beta,priors=list(beta=list(mu=0,sigma=1))) {   
  sum(loglikelihood_fast(beta,z-1,s$ptr(),K))+
  sum(dnorm(unlist(beta),priors$beta$mu,priors$beta$sigma,log=TRUE)) + 
  N * log(1/K)
}

brem.mcmc <- function(A,N,K,s,niter=5,model.type="full",mcmc.sd=.1,beta=NULL,z=NULL,gibbs="fast",mh=TRUE,outdir=getwd(),priors=list(beta=list(mu=0,sigma=1))) {
  llks <- rep(0,niter)
  M <- nrow(A)
  P <- 11
  param <- array(0,c(niter,K,K,P))
  zs <- NULL

  current <- array(rnorm(P*K^2,priors$beta$mu,priors$beta$sigma),c(P,K,K))
  current[7:11,,] <- 0
  
  if (!is.null(beta)) current <- beta
  if (is.null(z))    z <- sample(1:K,N,replace=TRUE)
  
  for (iter in 1:niter) {
    
#     for (k in 1:K) {
#       if (length(which(z==k))==0) {
#         other <- sample((1:K)[-k],1)
#         cat("attempting split",other,"into",k,"\n")
#         z.new <- z
#         z.new[chosen] <- k
#         current.new <- brem.mh(A,N,K,P,z,s,current,model.type,priors,mcmc.sd)
#       }
#     }
    
    olp <- brem.lpost.fast(A,N,K,z,s,current)
    
    z.new <- sample(1:K,N,replace=TRUE)
    cand <- current
    for (i in 1:5) {
      cand <- brem.mh(A,N,K,P,z.new,s,cand,model.type,priors,mcmc.sd)
    }
    clp <- brem.lpost.fast(A,N,K,z.new,s,cand)
    if (clp - olp > log(runif(1))) {
      z <- z.new
      current <- cand
      olp <- clp
      cat("Node ordering accepted:",z,"\n")
    }
    
    # For each effect sample via MH
    for (i in 1:2) {
      current <- brem.mh(A,N,K,P,z,s,current,model.type,priors,mcmc.sd)
    }
    
    # Sample empty cluster parameters from prior.  Only sample baserates.
    for (k in 1:K) {
      if (length(which(z==k))==0) {
        current[,k,] <- current[,,k] <- 0
        current[1,k,] <- current[1,,k] <- rnorm(K,priors$beta$mu,priors$beta$sigma)
        cat("sampling empty cluster params\n",current[1,,],"\n")
      }
    }
    
    if (gibbs=="slow") {
      # Gibbs sample assignments
      for (i in 1:N) {
        ps <- rep(0,K)
        for (k in 1:K) {
          z[i] <- k
          ps[k] <- brem.lpost.fast(A,N,K,z,s,current,priors)
        }
        ps <- exp(ps - max(ps))
        z[i] <- sample(1:K,size=1,prob=ps)
      }
    }
    if (gibbs=="fast") {
      z <- gibbs(current,z-1,s$ptr(),K)$z + 1
    }
    
    zs[[iter]] <- z
    param[iter,,,] <- current
    llks[iter] <- brem.lpost.fast(A,N,K,z,s,current,priors)
    
    cat("iter",iter,":",llks[iter],"z:",z,"\n")
    
    res <- list(z=z,beta=current,llks=llks,param=param,zs=zs,niter=niter)
    outfile <- paste(outdir,"/",model.type,".",K,".rdata",sep="")
    save(res,file=outfile)
  }
  return(res)
}
brem.mh <- function(A,N,K,P,z,s,current,model.type="baserates",priors,mcmc.sd=.1) {
  px <- c(1,1,1,1,1,1,1,0,0,0,0)  # TODO: Eventually allow MH updates on degree effects
  olp <- brem.lpost.fast(A,N,K,z,s,current,priors)
  
  cand <- current
  if (model.type=="baserates") {
    cand[1,,] <- cand[1,,] + rnorm(K^2,0,mcmc.sd)
    cand[-1,,] <- 0
    cand[1,1,1] <- 0  # identifiability?
    clp <- brem.lpost.fast(A,N,K,z,s,cand)
    if (clp - olp > log(runif(1))) {
      current <- cand
      olp <- clp
    }
  }
  if (model.type=="shared") {
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
      
      cand[1,1,1] <- 0  # identifiability?
      
      # MH sampler
      clp <- brem.lpost.fast(A,N,K,z,s,cand,priors)
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
      cand[1,1,1] <- 0  # identifiability?
      clp <- brem.lpost.fast(A,N,K,z,s,cand,priors)
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