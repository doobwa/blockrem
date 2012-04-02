library(brem)
source("pkg/R/brem.r")
require(ggplot2)
ppc.times <- function(sims) {
  d <- melt(sapply(sims,function(s) s$edgelist[,1]))
  colnames(d) <- c("event","sim","value")
  return(d)
}
ppc.pshift <- function(sims) {
  require(relevent)
  d <- mclapply(1:length(sims), function(i) {
    M <- nrow(sims[[i]]$edgelist)
    pshifts <- accum.ps(sims[[i]]$edgelist)[M,c(1,3,8:9,13)] # 10 is ab-xy
    data.frame(sim=i,stat=names(pshifts),value=pshifts)
  })
  d <- do.call(rbind,d)
  rownames(d) <- c()
  return(d)
}
ppc.degree <- function(sims) {
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
ppc.rem <- function(sims,stat) {
  switch(stat,
         "global"=ppc.times(sims),
         "pshift"=ppc.pshift(sims),
         "degree"=ppc.degree(sims))
}


# # Eckmann
# load("data/eckmann-small.rdata")
# load("results/eckmann-small-bk2/full.2.rdata")
#load("results/eckmann-small/full.2.rdata")
load(paste("data/",opts$dataset,".rdata",sep=""))
load(paste("results/",opts$dataset,"/full.2.rdata",sep=""))
fit$niter <- max(which(fit$llks!=0))
A <- as.matrix(A)
dimnames(A) <- NULL
library(multicore)

# What if we had been sampling more?
# fit$niter <- 50
# fit$param <- fit$param[150:200,,,]
# for (i in 1:50) fit$param[i,,,] + rnorm(P*K*K)

sims <- simulate(fit,nsim=50)

sims <- simulate(fit,train,M=nrow(test),nsim=100)
dimnames(test) <- NULL
obs <- list(list(edgelist=test,N=N))

# sims <- simulate(fit,nsim=1,conditioned=train,M=nrow(test))
# sims2 <- simulate(fit,nsim=1,return.lrm=TRUE)
obs <- list(list(edgelist=A[1:fit$M,],N=N))

tail(train)
tail(test)
tail(sims[[1]]$edgelist)

stats <- list()
stats$sim <- stats$obs <- list()

# Plot overall event vs. time
stats$sim$time <- ppc.rem(sims,"global")
stats$obs$time <- ppc.rem(obs,"global")
stats$sim$pshift <- ppc.rem(sims,"pshift")
stats$obs$pshift <- ppc.rem(obs,"pshift")
stats$sim$degree <- ppc.rem(sims,"degree")
stats$obs$degree <- ppc.rem(obs,"degree")

ggplot(stats$sim$time) + geom_line(alpha=.3,aes(x=event,y=value,group=sim)) + geom_line(data=stats$obs$time,aes(x=event,y=value),colour="red",size=2) + theme_bw() + labs(y="total elapsed time")
ggsave(paste("figs/",opts$dataset,"/postpred-time.pdf",sep=""),width=5,height=5)

ggplot(stats$sim$pshift) + geom_boxplot(aes(x=stat,y=value)) + geom_point(data=stats$obs$pshift,aes(x=stat,y=value),colour="red") + theme_bw()
ggsave(paste("figs/",opts$dataset,"/postpred-pshift.pdf",sep=""),width=5,height=5)


# OLD CODE FOR TINKERING
# 
# sapply(sims,function(s) s$edgelist[nrow(s$edgelist),1])
# bad <- which(sapply(sims,function(s) length(s)) != 4)
# if (length(bad) > 0) sims <- sims[-bad]
# 
# pred <- get.pred(train,A,test.ix,fit)
# lamhats <- sapply(1:2000,function(i) sum(exp(pred$lrm$train[i,,])))
# lamhats.test <- sapply(1:nrow(test),function(i) sum(exp(pred$lrm$test[i,,])))
# lamhats.train2 <- sapply(1:nrow(test),function(i) sum(exp(sims2[[1]]$lrm[i,,])))
# summary(c(pred$lrm$train[10,,]))
# summary(c(sims2[[1]]$lrm[10,,]))
