llk_indiv <- function(a,lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 0:(M-1)) {
    i <- sen[m+1]
    j <- rec[m+1]
    if (i==a | j==a | m==(M-1)) {
      llks[m+1] <- lrm[m+1,i+1,j+1]
      for (r in 0:(N-1)) {
        if (r != i & r!=j) {
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * sum(exp(lrm[m+1,i+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * sum(exp(lrm[m+1,r+1,i+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * sum(exp(lrm[m+1,j+1,r+1]))
          llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * sum(exp(lrm[m+1,r+1,j+1]))
        }
      }
      llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * sum(exp(lrm[m+1,i+1,r+1]))
      llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,r+1,i+1]))
    }
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  return(llks)
}

llk_slow <- function(lrm,times,sen,rec) {
  sen <- sen+1
  rec <- rec+1
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 1:M) diag(lrm[m,,]) <- -Inf
  llks[1] = lrm[1,sen[1],rec[1]]
  for (m in 2:(M-1)) {
    llks[m] <- lrm[m,sen[m],rec[m]] - (times[m]-times[m-1]) * sum(exp(lrm[m,,]))
  }
  return(llks)
}

test_taus_from_s <- function(times,sen,rec,N,M,P) {
  s <- new(brem$Stat,times,sen,rec,N,M,P)
  s$precompute()
  taus <- array(0,c(M,N,N))
  for (m in 0:(M-1)) {
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i != j) {
          taus[m+1,i+1,j+1] <- s$get_tau(m,i,j)
        }
      }
    }
  }
  return(taus)
}

test_taus <- function(lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- matrix(0,N,N)
  M <- length(times)
  taus <- array(0,c(M,N,N))
  for (m in 0:(M-2)) {
    for (i in 0:(N-1)) {
      for (j in 0:(N-1)) {
        if (i!=j) {
          taus[m+1,i+1,j+1] <- times[mp[i+1,j+1]+1]
        }
      }
    }
    i <- sen[m+1]
    j <- rec[m+1]
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
    
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        taus[m+1,i+1,j+1] <- times[mp[i+1,j+1]+1]
      }
    }
  }
  return(taus)
}
llk_fast_last <- function(lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- matrix(0,N,N)
  taus <- matrix(0,N,N)
  for (m in 1:(M-2)) {
    i <- sen[m]
    j <- rec[m]
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        taus[i+1,j+1] <- times[mp[i+1,j+1]+1]
        llks[i+1,j+1] <- (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  return(list(taus=taus,llks=llks,mp=mp))
}

# 0 based indexing on sen and rec
llk_slow <- function(lrm,times,sen,rec) {
  sen <- sen+1
  rec <- rec+1
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 1:M) diag(lrm[m,,]) <- -Inf
  llks[1] = lrm[1,sen[1],rec[1]]
  for (m in 2:(M-1)) {
    llks[m] <- lrm[m,sen[m],rec[m]] - (times[m]-times[m-1]) * sum(exp(lrm[m,,]))
  }
  return(llks)
}
llk_fast <- function(lrm,times,sen,rec) {
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  llks[1] <- lrm[1,sen[1]+1,rec[1]+1]
  for (m in 1:(M-2)) {
    i <- sen[m+1]
    j <- rec[m+1]
    llks[m+1] <- lrm[m+1,i+1,j+1]
    for (r in 0:(N-1)) {
      if (r != i & r!=j) {
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,r+1]+1]) * exp(lrm[m+1,i+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,i+1]+1]) * exp(lrm[m+1,r+1,i+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,r+1]+1]) * exp(lrm[m+1,j+1,r+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[r+1,j+1]+1]) * exp(lrm[m+1,r+1,j+1])
      }
    }
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[j+1,i+1]+1]) * exp(lrm[m+1,i+1,j+1])
    llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * exp(lrm[m+1,j+1,i+1])
    
    mp[i+1,] <- m
    mp[,i+1] <- m
    mp[j+1,] <- m
    mp[,j+1] <- m
  }
  m <- M-1
  for (i in 0:(N-1)) {
    for (j in 0:(N-1)) {
      if (i!=j) {
        #cat(i," ",j," ",times[mp[i+1,j+1]+1])
        llks[m+1] <- llks[m+1] - (times[m+1]-times[mp[i+1,j+1]+1]) * sum(exp(lrm[m+1,i+1,j+1]))
      }
    }
  }
  return(llks)
}

precomputeTau <- function(A,N) {
  M <- nrow(A)
  sen <- A[,2]
  rec <- A[,3]
  tau <- array(0,c(M,N,N))
  for (m in 1:(M-1)) {
    i <- sen[m]
    j <- rec[m]
    tau[(m+1):M,,i] <- m-1
    tau[(m+1):M,,j] <- m-1
    tau[(m+1):M,i,] <- m-1
    tau[(m+1):M,j,] <- m-1
  }
  return(tau)
}
get.indices <- function(A,N) {
  # Fast enough for our scale of problems
  # Returns 0-based index
  # Includes M-1 in all lists
  M <- nrow(A)
  xs <- vector("list",N)
  for (m in 1:M) {
    i <- A[m,2]
    j <- A[m,3]
    xs[[i]] <- c(xs[[i]],m - 1)
    xs[[j]] <- c(xs[[j]],m - 1)
  }
  for (i in 1:N) 
    xs[[i]] <- c(xs[[i]],M - 1)
  return(xs)
}

ranks <- function(edgelist,ratemats,...) {
  M <- nrow(edgelist)
  n <- dim(ratemats)[2]
  r <- rep(0,M)
  for (i in 1:M) {
    if (length(dim(ratemats))==2) rmat <- ratemats
    else rmat <- ratemats[i,,]
    rk <- matrix(rank(rmat,...),n,n)
    r[i] <- rk[edgelist[i,2],edgelist[i,3]]
  }
  return(r)
}
recall <- function(rs,top=c(1:20)) {
  data.frame(k=top,recall=sapply(top,function(k) mean(rs <= k)))
}

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
