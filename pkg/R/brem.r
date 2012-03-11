#' Simulate data from the block relational event model
#' @M number of events
#' @N number of actors
#' @z latent class assignments
#' @beta P x K x K array of parameters
simulate.brem <- function(M,N,z,beta) {
  
  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  lambda <- beta[1,z,z]
  lrm <- array(0,c(M,N,N))
  s <- InitializeStatisticsArray(N,P);
  z <- z-1
  for (m in 2:M) {
    
    i <- A[m-1,2]
    j <- A[m-1,3]
    
    s <- UpdateStatisticsArray(s,m-1,i-1,j-1,N,P)
    # Compute changes to lambda
    for (r in 1:N) {
      lambda[r,j] <- LogLambda(r-1,j-1,z[r],z[j],s,beta,N,K,P)
      lambda[j,r] <- LogLambda(j-1,r-1,z[j],z[r],s,beta,N,K,P)
      lambda[i,r] <- LogLambda(i-1,r-1,z[i],z[r],s,beta,N,K,P)
      lambda[r,i] <- LogLambda(r-1,i-1,z[r],z[i],s,beta,N,K,P)
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
  lrm <- RemLogIntensityArray(beta,times,sen,rec,z,N,M,K,P)
  for (i in 1:M) diag(lrm[i,,]) <- -Inf
  return(lrm)
}
brem.lrm.fast <- function(M,s,z,beta) {
  K <- dim(beta)[2]
  lrm <- LogIntensityArrayPc(beta,z-1,s$ptr(),K)
  for (i in 1:(dim(lrm)[1])) diag(lrm[i,,]) <- -Inf
  return(lrm)
}
brem.lrm.fast.subset <- function(s,z,beta,ix) {
  K <- dim(beta)[2]
  ix <- as.integer(ix)
  lrm <- LogIntensityArrayPcSubset(beta,z-1,s$ptr(),K,ix-1)
  for (i in 1:(dim(lrm)[1])) diag(lrm[i,,]) <- -Inf
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
    tmp <- LogIntensityArray(beta,times,sen,rec,z,N,M,K,P)
    return(RemLogLikelihoodFromArray(tmp,times,sen,rec,N,M))
  } else {
    return(RemLogLikelihood(beta,times,sen,rec,z,N,M,K,P))
  }
}
brem.mle <- function(A,N,K,P,z,beta=NULL) {
  fn <- function(par) {
    beta <- array(par,c(P,K,K))
    beta[7:13,,] <- 0
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
  sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K)) +
  sum(dnorm(unlist(beta),priors$beta$mu,priors$beta$sigma,log=TRUE)) + 
  N * log(1/K)
}
#' @A event history matrix
#' @z cluster assignments (1-based)
brem.lpost.fast.block <- function(A,N,K,z,s,beta,k1,k2,priors=list(beta=list(mu=0,sigma=1))) {   
  M <- nrow(A)
  zs <- z[A[,2]]
  zr <- z[A[,3]]
  ix <- which(zs==k1 | zr==k2) - 1  # 0 based indexing for c++
  if (length(ix) > 0) {
    ix <- setdiff(ix,c(0,M-1))
  }
  llk <- sum(dnorm(unlist(beta),priors$beta$mu,priors$beta$sigma,log=TRUE)) + N * log(1/K)
  if (length(ix)==0) {
    return(llk)
  } else {
    return(llk + sum(RemLogLikelihoodPcSubset(beta,z-1,s$ptr(),K,ix))) 
  }    
}

brem.mcmc <- function(A,N,K,s,niter=5,model.type="full",mcmc.sd=.1,m=20,beta=NULL,z=NULL,gibbs="fast",mh=FALSE,outfile=NULL,priors=list(beta=list(mu=0,sigma=1)),verbose=FALSE,px=NULL) {
  llks <- rep(0,niter)
  M <- nrow(A)
  P <- 13
  param <- array(0,c(niter,P,K,K))
  zs <- NULL

  current <- array(rnorm(P*K^2,priors$beta$mu,priors$beta$sigma),c(P,K,K))
  for (p in which(px==0)) {
    current[p,,] <- 0
  }
  current[1,1,1] <- 0  # identifiability?
  
  if (!is.null(beta)) current <- beta
  if (is.null(z))     z <- sample(1:K,N,replace=TRUE)
  if (is.null(px))    px <- rep(1,P-1)  # last feature is event index m
  
  olp <- brem.lpost.fast(A,N,K,z,s,current)
  
  deltas <- matrix(0,niter,5)
  for (iter in 1:niter) {
    stime <- proc.time()
    # For each effect sample via MH
    if (mh) {
      for (i in 1:2) {
        res <- brem.mh(A,N,K,P,z,s,current,px,model.type,priors,mcmc.sd,olp)
        current <- res$current
        olp <- res$olp
      }
    } else {
      res <- brem.slice(A,N,K,P,z,s,current,px,model.type,priors,olp,m=m)
      current <- res$current
      olp <- res$olp
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
      cat("gibbs")
      z <- gibbs(1:N-1,current,z-1,s$ptr(),K)$z + 1
    }
    
    zs[[iter]] <- z
    param[iter,,,] <- current
    llks[iter] <- brem.lpost.fast(A,N,K,z,s,current,priors)
    
    etime <- proc.time()
    cat("\niter",iter,":",llks[iter],"z:",z,"\n")
    cat(current[,1,1],"\n")
    
    deltas[iter,] <- etime-stime
    if (K>1) cat(current[,2,2],"\n")
    cat("time:",deltas[iter],"\n")
    
    res <- list(z=z,beta=current,llks=llks,param=param,zs=zs,niter=niter,deltas=deltas)
    if (!is.null(outfile)) {
      save(res,file=outfile)
    }
  }
  return(res)
}
brem.mh <- function(A,N,K,P,z,s,current,px,model.type="baserates",priors,mcmc.sd=.1,olp=NULL) {
  if (is.null(olp)) {
    olp <- brem.lpost.fast(A,N,K,z,s,current,priors)
  }
  
  cand <- current
  if (model.type=="baserates") {
    cand[1,,] <- cand[1,,] + rnorm(K^2,0,mcmc.sd)
    cand[-1,,] <- 0   # 
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
#     for (p in which(px==1)) {
#       cand[p,,]  <- cand[p,,] + rnorm(K^2,0,mcmc.sd)
#       cand[1,1,1] <- 0  # identifiability?
#       clp <- brem.lpost.fast(A,N,K,z,s,cand,priors)
#       if (clp - olp > log(runif(1))) {
#         current <- cand
#         olp <- clp
#       }
#     }
    for (k1 in 1:K) {
      for (k2 in 1:K) {
        olp <- brem.lpost.fast.block(A,N,K,z,s,current,k1,k2,priors)
        for (p in which(px==1)) {
          if (!(p==1 & k1==1 & k2==1)) {
            cat(".")
            cand[p,k1,k2]  <- cand[p,k1,k2] + rnorm(1,0,mcmc.sd)
            clp <- brem.lpost.fast.block(A,N,K,z,s,cand,k1,k2,priors)
            if (clp - olp > log(runif(1))) {
              current <- cand
              olp <- clp
            }
          }
        }
      }
    }
  }
  return(list(current=current,olp=olp))
}

brem.slice <- function(A,N,K,P,z,s,beta,px,model.type="baserates",priors,olp=NULL,m=20) {
  
  #' Needs p,k1,k2,A,N,K,z,s,beta,priors,model.type
  slicellk <- function(x) {
    beta[p,k1,k2] <- x
    if (model.type=="shared") {
      beta <- use.first.blocks(beta)
    }
    #brem.lpost.fast.block(A,N,K,z,s,beta,k1,k2,priors)
    brem.lpost.fast(A,N,K,z,s,beta,priors)
  }
  
  use.first.blocks <- function(beta) {
    for (p in 1:P) {
    for (j1 in 1:K) {
      for (j2 in 1:K) {
        if (j1==j2) {
          beta[p,j1,j2] <- beta[p,1,1]
        } 
        if (j1!=j2) {
          beta[p,j1,j2] <- beta[p,1,2]
        }
      } 
    }
    }
    return(beta)
  }
  
  olp <- NULL
  if (is.null(olp)) {
    olp <- brem.lpost.fast(A,N,K,z,s,beta,priors)
  }
  
  kx1 <- 1:K
  kx2 <- 1:K
  beta[1,1,1] <- 0
  if (model.type=="baserates") {
    px <- 1
    beta[-1,,] <- 0
  } else if (model.type=="full") {
  } else if (model.type=="shared" & K>1) {
    kx1 <- 1
    kx2 <- 1:2
  }
  
  olps <- olp
  #olp <- brem.lpost.fast(A,N,K,z,s,beta,priors)
  for (k1 in kx1) {
    for (k2 in kx2) {
      #olp <- brem.lpost.fast.block(A,N,K,z,s,beta,k1,k2,priors)
      olp <- brem.lpost.fast(A,N,K,z,s,beta,priors)
      for (p in which(px==1)) {
        cat(".")
        if (!(k1==1 & k2==1 & p==1)) {  # identifiability
          newval <- uni.slice.alt(beta[p,k1,k2],slicellk,gx0=olp,m=m)
          beta[p,k1,k2] <- newval
          if (model.type=="shared") beta <- use.first.blocks(beta)
          olp <- attr(newval,"log.density")
          olps <- c(olps,olp)
        }
      }
    }
  }

  return(list(current=beta,olp=olp,olps=olps))
}

# Sample empty cluster parameters from prior.  Only sample baserates.
#     for (k in 1:K) {
#       if (length(which(z==k))==0) {
#         current[,k,] <- current[,,k] <- 0
#         current[1,k,] <- current[1,,k] <- rnorm(K,priors$beta$mu,priors$beta$sigma)
#         cat("sampling empty cluster params\n",current[1,,],"\n")
#       }
#     }
