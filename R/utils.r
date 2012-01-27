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