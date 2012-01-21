# Deprecated
build.brem <- function(A,N) {
  implementedEffects <- c("RSndSnd","RRecSnd","PSAB-BA","PSAB-BY","PSAB-XA","PSAB-XB","PSAB-XY","PSAB-AY")
  M <- nrow(A)
  P <- length(implementedEffects)
  s <- array(0,c(M,N,N,P))
  dimnames(s) <- list(NULL,NULL,NULL,implementedEffects)
  rrs <- rss <- matrix(0,N,N)
  for (m in 1:(M-1)) {
   
    a <- A[m,2]
    b <- A[m,3]
    
    # P-shift effects
    s[m,b,a,"PSAB-BA"] <- 1
    s[m,b,-a,"PSAB-BY"] <- 1
    s[m,-b,a,"PSAB-XA"] <- 1
    s[m,-a,b,"PSAB-XB"] <- 1
    s[m,-c(a,b),-c(a,b),"PSAB-XY"] <- 1
    s[m,a,-b,"PSAB-AY"] <- 1
    
    # Recency effects
    ix <- which(rss[a,] > 0)
    if (length(ix) > 0) {
      s[m,a,ix,"RSndSnd"] <- 1/(rss[a,ix]) 
    }
    ix <- which(rrs[a,] > 0)
    if (length(ix) > 0) {
      s[m,a,ix,"RRecSnd"] <- 1/(rrs[a,ix])
    }

    
    # Update recency statistics: a just sent to b
    # rrs[a,b]: rank of b in a's list of people recently received from
    jx <- which(rrs[b,] > 0)  # b just received from a
    rrs[b,jx] <- rrs[b,jx] + 1  # src to sent to rec
    rrs[b,a] <- 1
    # rss[a,b]: rank of b in a's list of people recently sent to
    jx <- which(rss[a,] > 0)  # a just sent to b
    rss[a,jx] <- rss[a,jx] + 1
    rss[a,b] <- 1
  }
  return(s)
}
block.ps <- function(A,z) {
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
simulate.brem <- function(M,N,z,beta) {
  implementedEffects <- c("Intercept","RSndSnd","RRecSnd","PSAB-BA","PSAB-BY","PSAB-XA","PSAB-XB","PSAB-XY","PSAB-AY")
  if (dim(beta)[3] != length(implementedEffects)) stop("wrong dimensions for parameter array")
  P <- length(implementedEffects)
  s <- array(0,c(M,N,N,P))
  dimnames(s) <- list(NULL,NULL,NULL,implementedEffects)
  rrs <- rss <- matrix(0,N,N)
  
  zmemb <- lapply(1:max(z),function(x) which(z==x))
  
  # Use baserate to initialize lambda
  lambda <- beta[z,z,1]
  x <- array(0,c(N,N,P))
    
  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  # s is sufficient statistics immediately following last event
  for (m in 1:(M-1)) {
    
    a <- A[m,2]
    b <- A[m,3]
    
    # Compute most recent statistics
    s[m,,,"Intercept"] <- 1
    
    # P-shift effects
    s[m,b,a,"PSAB-BA"] <- 1
    s[m,b,-a,"PSAB-BY"] <- 1
    s[m,-b,a,"PSAB-XA"] <- 1
    s[m,-a,b,"PSAB-XB"] <- 1
    s[m,-c(a,b),-c(a,b),"PSAB-XY"] <- 1
    s[m,a,-b,"PSAB-AY"] <- 1
    
    # Recency effects
    ix <- which(rss[a,] > 0)
    if (length(ix) > 0) {
      s[m,a,ix,"RSndSnd"] <- 1/(rss[a,ix]) 
    }
    ix <- which(rrs[a,] > 0)
    if (length(ix) > 0) {
      s[m,a,ix,"RRecSnd"] <- 1/(rrs[a,ix])
    }
    
    ix <- unique(zmemb[[z[a]]],zmemb[[z[b]]])
    for (i in ix) {
      for (j in ix) {
        x[i,j,] <- s[m,i,j,]  # s_ij^* in notes
        lambda[i,j] <- beta[z[i],z[j],] %*% x[i,j,]
      }
    }
       
    diag(lambda) <- -Inf
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    src <- cells[drawcell,1]
    rec <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    A <- rbind(A,c(time,src,rec))
    
    # Update recency statistics: a just sent to b
    # rrs[a,b]: rank of b in a's list of people recently received from
    jx <- which(rrs[b,] > 0)  # b just received from a
    rrs[b,jx] <- rrs[b,jx] + 1  # src to sent to rec
    rrs[b,a] <- 1
    # rss[a,b]: rank of b in a's list of people recently sent to
    jx <- which(rss[a,] > 0)  # a just sent to b
    rss[a,jx] <- rss[a,jx] + 1
    rss[a,b] <- 1
  }
  return(A)
}