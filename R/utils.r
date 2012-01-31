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
