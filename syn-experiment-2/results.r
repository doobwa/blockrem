load("pkg/experiment/res.rdata")

niter <- 100
ix <- (1:nrow(s))#[-ix]
lps <- lapply(ix,function(i) {
  x <- data.frame(iter=1:niter,s[i,],lps=res[[i]]$lps)
  return(x)
})
df <- do.call(rbind,lps)
df$case <- factor(df$case)
qplot(iter,lps,data=df,geom="line",colour=case)+facet_grid(sig~do.sm + do.extra)

df <- subset(df,sig==1)
qplot(iter,lps,data=df,geom="line",colour=case)+facet_grid(do.sm ~ do.extra)

tmp <- lapply(1:nrow(s),function(i) {
  x <- s[i,]
  x$mean.k <- mean(sapply(res[[i]]$samples,function(x) dim(x$phi)[2]))
  x$final.lp <- res[[i]]$lps[niter]
  x$acc.rate <- mean(res[[i]]$acc)
  return(x)
})
tmp <- do.call(rbind,tmp)

x <- subset(tmp,!do.sm)

y <- subset(tmp,do.sm==TRUE)
