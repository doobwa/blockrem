# Deprecated
build.brem <- function(A,N) {
  implementedEffects <- c("RSndSnd","RRecSnd","PSAB-BA","PSAB-BY","PSAB-XA","PSAB-XB","PSAB-XY","PSAB-AY")
  M <- nrow(A)
  P <- length(implementedEffects)
  s <- array(0,c(M,N,N,P))
  dimnames(s) <- list(NULL,NULL,NULL,implementedEffects)
  rrs <- rss <- mp <- matrix(0,N,N)
  for (m in 1:(M-1)) {
   
    a <- A[m,2]
    b <- A[m,3]
    
    # P-shift effects
    s[m,b,a,"PSAB-BA"] <- 1
    s[m,b,-a,"PSAB-BY"] <- 1
    s[m,-b,a,"PSAB-XA"] <- 1
    s[m,-a,b,"PSAB-XB"] <- 1
    s[m,a,-b,"PSAB-AY"] <- 1
    
  }
  return(s)
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
simulate.brem <- function(M,N,z,beta) {
  implementedEffects <- c("Intercept","RSndSnd","RRecSnd","PSAB-BA","PSAB-BY","PSAB-XA","PSAB-XB","PSAB-AY")
  if (dim(beta)[3] != length(implementedEffects)) stop("wrong dimensions for parameter array")
  P <- length(implementedEffects)
  
  # Use baserate to initialize lambda
  
  computeLambda <- function(i,j,a,b,beta) {
    lam = beta[1]
    if (i==b & j==a) { # ab-ba
       lam  <- lam + beta[2];  
    }
    if (i==b & j!=a) { # ab-by
       lam <- lam + beta[3]; 
    }
    if (i!=b & j==a) { # ab-xaF
       lam <- lam + beta[4];
    }
    if (i!=a & j==b) { # ab-xb
       lam <- lam + beta[5];
    }
    if (i==b & j!=b) { # ab-ay
       lam <- lam + beta[6];
    }
    return(lam)
  }
    
  # start with a event from 1 to 2
  time <- 0
  A <- matrix(c(time,1,2),1,3)
  lambda <- matrix(beta[1],N,N)
  # s is sufficient statistics immediately following last event
  for (m in 1:(M-1)) {
    
    # Draw the next event
    diag(lambda) <- -Inf
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    i <- cells[drawcell,1]
    j <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    A <- rbind(A,c(time,i,j))
    
    # Compute changes to lambda
    for (r in 1:N) {
      lambda[r,j] <- computeLambda(r,j,i,j,beta[z[r],z[j],])
    }
    for (r in 1:N) {
      lambda[i,r] <- computeLambda(i,r,i,j,beta[z[i],z[r],])
    }
    
  }
  return(A)
}