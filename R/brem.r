
dyad.ps <- function(A,N) {
  require(relevent)
  x <- array(0,c(N,N,6))
  mp <- matrix(0,N,N)
  M <- nrow(A)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    # Find last time lambda_ij changed
    if (mp[i,j]!=0) {
      a <- A[mp[i,j],2]
      b <- A[mp[i,j],3]
      x[i,j,] <- x[i,j,] + pshift(i,j,a,b)
    }
    mp[i,] <- m
    mp[,i] <- m
    mp[j,] <- m
    mp[,j] <- m
  }
  dimnames(x)[[3]] <- c("AB-BA","AB-BY","AB-XA","AB-XB","AB-AY","AB-AB")
  return(x)
}

node.ps <- function(A,N) {
  require(relevent)
  x <- matrix(0,N,2)
  p <- matrix(0,N,5)
  M <- nrow(A)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    a <- x[i,1]
    b <- x[i,2]
    p[i,] <- p[i,] + pshift(i,j,a,b)
    a <- x[j,1]
    b <- x[j,2]
    p[j,] <- p[j,] + pshift(i,j,a,b)
    
    # Update last event for i and j
    x[i,] <- x[j,] <- c(i,j)
  }
  colnames(p) <- c("AB-BA","AB-BY","AB-XA","AB-XB","AB-AY")
  return(p)
}

computeLambda <- function(i,j,a,b,beta,px) {
  lam <- 0
  if (px[1]) {
    lam <- lam + beta[1]
  }
  lam <- lam + sum(px[-1] * beta[-1] * pshift(i,j,a,b))
  return(lam)
}
pshift <- function(i,j,a,b) {
  x <- rep(0,6)
  if (i!=a & i==b & j==a & j!=b) { # ab-ba
    x[1] <- 1
  }
  if (i!=a & i==b & j!=a & j!=b) { # ab-by
    x[2] <- 1
  }
  if (i!=a & i!=b & j==a & j!=b) { # ab-xa
    x[3] <- 1
  }
  if (i!=a & i!=b & j!=a & j==b) { # ab-xb
    x[4] <- 1
  }
  if (i==a & i!=b & j!=a & j!=b) { # ab-ay
    x[5] <- 1
  }
#   if (i==a & i!=b & j!=a & j==b) { # ab-ab
#     x[6] <- 1
#   }
  return(x)
}
block.ps <- function(A,z) {
  require(relevent)
  s <- A[,2]
  r <- A[,3]
  zs <- sort(unique(z))
  setup <- expand.grid(i=zs,j=zs)
  ds <- lapply(1:nrow(setup),function(k) {
    ix <- which(z[s] == setup$i[k] & z[r]==setup$j[k])
    pshifts <- accum.ps(A[ix,])[length(ix),c(1,3,8,9,10,13)]
    data.frame(i=setup$i[k],j=setup$j[k],pshift=names(pshifts),value=pshifts)
  })
  ds <- do.call(rbind,ds)
  rownames(ds) <- c()
  return(ds)
}
simulate.brem <- function(M,N,z,beta,px) {
  implementedEffects <- c("Intercept","PSAB-BA","PSAB-BY","PSAB-XA","PSAB-XB","PSAB-AY","PSAB-AB")
  if (dim(beta)[3] != length(implementedEffects)) stop("wrong dimensions for parameter array")
  P <- length(implementedEffects)
  
  # Use baserate to initialize lambda

  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  lambda <- beta[z,z,1]#matrix(beta[1],N,N)
  lrm <- array(0,c(M,N,N))
  for (m in 2:M) {
    
    i <- A[m-1,2]
    j <- A[m-1,3]
    # Compute changes to lambda
    for (r in 1:N) {
      lambda[r,j] <- computeLambda(r,j,i,j,beta[z[r],z[j],],px)
    }
    for (r in 1:N) {
      lambda[j,r] <- computeLambda(j,r,i,j,beta[z[j],z[r],],px)
    }
    for (r in 1:N) {
      lambda[i,r] <- computeLambda(i,r,i,j,beta[z[i],z[r],],px)
    }
    for (r in 1:N) {
      lambda[r,i] <- computeLambda(r,i,i,j,beta[z[r],z[i],],px)
    }
    
    # Draw the next event
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    i <- cells[drawcell,1]
    j <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    A <- rbind(A,c(time,i,j))
    
    diag(lambda) <- -Inf
    lrm[m,,] <- lambda
  }
  return(list(A=A,lrm=lrm))
}

drem.lrm <- function(A,N,beta,ix,jx,px) {
  M <- nrow(A)
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  ix <- ix-1
  jx <- jx-1
  drem$lrm(beta,times,sen,rec,ix,jx,px,N,M)
}
# ix: unique senders
# jx: unique receivers
drem.llk <- function(A,N,beta,ix,jx,px) {
  M <- nrow(A)
  times <- A[,1]
  sen <- A[,2]-1
  rec <- A[,3]-1
  ix <- ix-1
  jx <- jx-1
  drem$llk(beta,times,sen,rec,ix,jx,px,N,M)
}
drem.mle <- function(A,N,beta,ix,jx,px) {
  fn <- function(par) {
    p <- array(0,c(1,1,length(par)))
    p[1,1,] <- par
    -drem.llk(A,par,ix,jx,px)
  }
  optim(as.vector(beta),fn)$par
}
brem.lrm <- function(A,N,K,z,beta,px) {
  M <- nrow(A)
  lrm <- array(0,c(M,N,N))
  z1 <- z[A[,2]]
  z2 <- z[A[,3]]
  zs <- lapply(1:K,function(k) which(z==k))
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      #ix <- which(z1 == k1 & z2 == k2)
      #if (length(ix) > 0 & length(zs[[k1]]) > 0 & length(zs[[k2]]) > 0) {
        #B <- A[ix,]
        #lrm[ix,,] <- lrm[ix,,] + drem.lrm(B,N,beta[k1,k2,],zs[[k1]],zs[[k2]],px)
      #}
      lrm <- lrm + drem.lrm(A,N,beta[k1,k2,],zs[[k1]],zs[[k2]],px)
    }
  }
  lrm
}
brem.llk <- function(A,N,K,z,beta,px) {
  llks <- matrix(0,K,K)
  z1 <- z[A[,2]]
  z2 <- z[A[,3]]
  zs <- lapply(1:K,function(k) which(z==k))
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      ix <- which(z1 == k1 & z2 == k2)
      if (length(ix) > 0 & length(zs[[k1]]) > 0 & length(zs[[k2]]) > 0) {
        B <- A[ix,]
        llks[k1,k2] <- drem.llk(B,N,beta[k1,k2,],zs[[k1]],zs[[k2]],px)
      }
    }
  }
  llks
}
brem.lpost <- function(A,N,K,z,beta,px) {
  llks <- brem.llk(A,N,K,z,beta,px)
  lprior <- sum(dnorm(beta,0,1,log=TRUE)) + N * log(1/K)
  sum(llks)+lprior
}

brem.mcmc <- function(A,N,K,px,niter=5,model.type="full",mcmc.sd=.1,beta=NULL,z=NULL,gibbs=TRUE,mh=TRUE) {
  llks <- rep(0,niter)
  M <- nrow(A)
  P <- length(px)
  param <- array(0,c(niter,K,K,P))
  mu <- log(M/A[M,1]/N/N)  # MLE for beta for uniform rate across all events
  sigma <- .5
  current <- array(rnorm(K^2*P,mu,sigma),c(K,K,P))
  if (is.null(z))    z <- sample(1:K,N,replace=TRUE)
  if (!is.null(beta)) current <- beta
  
  for (iter in 1:niter) {
    
    # For each effect sample via MH
    first.iter <- (iter==1)
    current <- brem.mh(A,N,K,P,z,current,px,model.type,mcmc.sd,first.iter)
    
    # Sample empty cluster parameters from prior
    for (k in 1:K) {
      if (length(which(z==k))==0) {
        current[k,,] <- rnorm(K*P,0,5)
        current[,k,] <- rnorm(K*P,0,5)
      }
    }
    
    if (gibbs) {
      # Gibbs sample assignments
      for (i in 1:N) {
        ps <- rep(0,K)
        for (k in 1:K) {
          z[i] <- k
          ps[k] <- brem.lpost(A,N,K,z,current,px)
        }
        ps <- exp(ps - max(ps))
        z[i] <- sample(1:K,size=1,prob=ps)
      }
    }
    
    param[iter,,,] <- current
    llks[iter] <- brem.lpost(A,N,K,z,current,px)
    
    cat("iter",iter,":",llks[iter],"z:",table(z),"\n")
    
  }
  return(list(z=z,llks=llks,param=param,beta=current))
}
brem.mh <- function(A,N,K,P,z,current,px,model.type="baserates",mcmc.sd=.1,first.iter=FALSE) {
  if (first.iter) olp <- -Inf
  else olp <- brem.lpost(A,N,K,z,current,px)
  cand <- current
  if (model.type=="baserates") {
    cand[,,1] <- cand[,,1] + rnorm(K^2,0,mcmc.sd)
    cand[,,-1] <- 0
    px <- c(1,rep(0,7))
    clp <- brem.lpost(A,N,K,z,cand,px)
    if (clp - olp > log(runif(1))) {
      current <- cand
      olp <- clp
    }
  }
  if (model.type=="diag.rem") {
    a <- cand[1,1,] + rnorm(P,0,mcmc.sd)
    b <- cand[1,2,] + rnorm(P,0,mcmc.sd)
    for (k1 in 1:K) {
      for (k2 in 1:K) {
        if (k1==k2) {
          cand[k1,k2,] <- a
        } else {
          cand[k1,k2,] <- b
        }
      }
    }
    clp <- brem.lpost(A,N,K,z,cand,px)
    if (clp - olp > log(runif(1))) {
      current <- cand
      olp <- clp
    }
  }
  if (model.type=="full") {
    for (p in which(px==1)) {
      cand <- current
      cand[,,p]  <- cand[,,p] + rnorm(K^2,0,mcmc.sd)
      clp <- brem.lpost(A,N,K,z,cand,px)
      if (clp - olp > log(runif(1))) {
        current <- cand
        olp <- clp
      }
    }
  }
  return(current)
}