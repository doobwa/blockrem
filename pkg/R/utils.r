plotmat <- function(mat,labels=c("Sender","Receiver"),limits=c(0,max(mat[,3])),color="black") {
  require(ggplot2)
  colnames(mat)[3] <- "value"
  qq <- ggplot(mat,aes(x=X2,y=X1)) +
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low="white",high=color,limits=limits)+
    theme_bw() + labs(x=labels[2],y=labels[1],fill="Probability")+
    coord_equal(ratio=1) +
    opts(legend.position = "none",
         panel.grid.minor=theme_blank(),
         panel.grid.major=theme_blank())
  qq
}

#' Full lograte array using computeLambdaFast.  lrm[m,i,j] is lambda_ij prior m'th event.  
#' All lambda=0 for m=0.
lrm_slow <- function(beta,z,s,M,N,K,P) {
  lrms <- array(0,c(M,N,N))
  for (m in 2:M) {
    for (i in 1:N) {
      for (j in 1:N) {
        if (i!=j) {
          x <- s$get_s(m-2,i-1,j-1)  # why m-2?  because this should be for the rateprior to current event m-1
          lrms[m,i,j] <- compute_lambda_fast(i,j,z[i],z[j],x,beta,N,K,P)
        } 
      }
    }
  }
  return(lrms)
}
  
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



test_taus_from_s <- function(times,sen,rec,N,M,P) {
  s <- new(RemStat,times,sen,rec,N,M,P)
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

test_taus <- function(lrm,times,sen,rec,M,N) {
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
llk_slow <- function(lrm,times,sen,rec,M,N) {
  sen <- sen+1
  rec <- rec+1
  mp <- matrix(0,N,N)
  llks <- rep(0,M)
  for (m in 1:M) diag(lrm[m,,]) <- -Inf
  llks[1] = lrm[1,sen[1],rec[1]]
  for (m in 2:M) {
    llks[m] <- lrm[m,sen[m],rec[m]] - (times[m]-times[m-1]) * sum(exp(lrm[m,,]))
  }
  return(llks)
}
llk_fast <- function(lrm,times,sen,rec,M,N) {
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

#' Compute the rank of the each event in the observed edgelist given an array of intensities
#' @edgelist Mx3 matrix of event times, senders and recipients
#' @ratemats MxNxN array of intensities.  
ranks.fast <- function(edgelist,ratemats,...) {
  M <- nrow(edgelist)
  N <- dim(ratemats)[2]
  r <- rep(0,M)
  for (m in 1:M) {
    i <- edgelist[m,2]
    j <- edgelist[m,3]
    k <- N*(j-1) + i
    r[m] <- rank(ratemats[i,,],...)[k]
  }
  return(r)
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

#' Compute the recall metric using a vector of ranks
#' @rs vector of observed ranks
#' @top vector of cutpoints to compute
recall <- function(rs,top=c(1:20)) {
  data.frame(k=top,recall=sapply(top,function(k) mean(rs <= k)))
}

dyad.ps <- function(A,N) {
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
rcatlog <- function(ps) {
  ps <- ps - max(ps)
  ps <- exp(ps)
  ps <- ps/sum(ps)
  ps <- cumsum(ps)
  u <- runif(1)
  for (i in 1:length(ps)) {
    if (u < ps[i]) return(i)
  }
}
plot.blockmodel <- function(B,z) {
  nmap <- order(z)
  s <- match(B[,2],nmap)
  r <- match(B[,3],nmap)
  plot(s,r,pch=".",xlab="sender",ylab="recipient")
  cpoints <- sapply(1:K,function(k) min(which(sort(z)==k)))
  abline(v=cpoints[-1])
  abline(h=cpoints[-1])
}

#' Compute an (M x N x N) array of previous events for each dyad.  
#' @edgelist event history matrix
#' @n total number of nodes
ratemat.online <- function(edgelist,n) {
  M <- nrow(edgelist)
  rms <- array(0,c(M,n,n))
  for (i in 2:M) {
    rms[i:M,edgelist[i-1,2],edgelist[i-1,3]] <- rms[i:M,edgelist[i-1,2],edgelist[i-1,3]] + 1
    diag(rms[i,,]) <- -Inf
  }
  return(rms)
}

#' Compute an (M.test x N x N) array of intensities for each dyad using the marginals from a training set.
#' @train event history
#' @test  event history to compute the array for
#' @N     number of nodes
ratemat.from.marginals <- function(train,test,N) {
  x <- table(factor(train[,2],1:N),factor(train[,3],1:N))
  rowrates <- rowSums(x)
  colrates <- colSums(x)
  r <- rowrates %*% t(colrates)
  r <- r/sum(r)
  diag(r) <- -Inf
  M <- nrow(test)
  lrm <- array(0,c(M,N,N))
  for (i in 1:M) lrm[i,,] <- r
  return(lrm)
}

#' Compute likelihood using the log intensity array.
#' @lrm log intensity array
#' @times 
#' @sen
#' @rec
#' @z
#' @N
#' @M
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
