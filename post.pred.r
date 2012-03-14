library(brem)
source("pkg/R/brem.r")
require(ggplot2)
post.pred.times <- function(sims) {
  d <- melt(sapply(sims,function(s) s$edgelist[,1]))
  colnames(d) <- c("event","sim","value")
  return(d)
}
post.pred.pshift <- function(sims) {
  require(relevent)
  d <- mclapply(1:length(sims), function(i) {
    M <- nrow(sims[[i]]$edgelist)
    pshifts <- accum.ps(sims[[i]]$edgelist)[M,c(1,3,8:10,13)]
    data.frame(sim=i,stat=names(pshifts),value=pshifts)
  })
  d <- do.call(rbind,d)
  rownames(d) <- c()
  return(d)
}
post.pred.degree <- function(sims) {
  require(network)
  require(sna)
  degrees <- mclapply(1:length(sims), function(i) {
    net <- network(sims[[i]]$edgelist[,2:3])
    N.i <- sims[[i]]$N
    rbind(
      data.frame(sim=i,stat="indegree",actor=1:N.i,obs=degree(net,cmode="indegree")),
      data.frame(sim=i,stat="outdegree",actor=1:N.i,obs=degree(net,cmode="outdegree")),
      data.frame(sim=i,stat="degree",actor=1:N.i,obs=degree(net)))
  })
  do.call(rbind,degrees)
}
post.pred.stat <- function(sims,stat) {
  switch(stat,
         "global"=post.pred.times(sims),
         "pshift"=post.pred.pshift(sims),
         "degree"=post.pred.degree(sims))
}

# Eckmann
load("data/eckmann-small.rdata")
load("results/eckmann-small-bk2/full.2.rdata")
#load("results/eckmann-small/full.2.rdata")
A <- as.matrix(A)
dimnames(A) <- NULL
library(multicore)
sims <- simulate(fit,nsim=50)
obs <- list(list(edgelist=A[1:fit$M,],N=N))

stats <- list()
stats$sim <- stats$obs <- list()

# Plot overall event vs. time
stats$sim$time <- post.pred.stat(sims,"global")
stats$obs$time <- post.pred.stat(obs,"global")
stats$sim$pshift <- post.pred.stat(sims,"pshift")
stats$obs$pshift <- post.pred.stat(obs,"pshift")
stats$sim$degree <- post.pred.stat(sims,"degree")
stats$obs$degree <- post.pred.stat(obs,"degree")

save(stats,file="figs/eckmann-small/postpred.rdata")

load("figs/eckmann-small/postpred.rdata")

ggplot(stats$sim$time) + geom_line(alpha=.5,aes(x=event,y=value,group=sim)) + geom_line(data=stats$obs$time,aes(x=event,y=2*value),colour="red",size=2) + theme_bw() + labs(y="total elapsed time")

ggplot(stats$sim$pshift) + geom_boxplot(aes(x=stat,y=value)) + geom_point(data=stats$obs$pshift,aes(x=stat,y=value),colour="red") + theme_bw()

