
ppc.rem <- function(sims,obs) {
  obs <- list(list(edgelist=obs,N=N))
  stats <- list()
  stats$sim <- stats$obs <- list()
  stats$sim$time <- ppc.stat(sims,"global")
  stats$obs$time <- ppc.stat(obs,"global")
  stats$sim$pshift <- ppc.stat(sims,"pshift")
  stats$obs$pshift <- ppc.stat(obs,"pshift")
#  stats$sim$degree <- ppc.stat(sims,"degree")
#  stats$obs$degree <- ppc.stat(obs,"degree")
  return(stats)
}

ppc.stat <- function(sims,stat) {
    switch(stat,
         "global"=ppc.times(sims),
         "pshift"=ppc.pshift(sims),
         "degree"=ppc.degree(sims))
}
  
ppc.times <- function(sims) {
  d <- melt(sapply(sims,function(s) s$edgelist[,1]))
  colnames(d) <- c("event","sim","value")
  return(d)
}
ppc.pshift <- function(sims,chosen=c(1,3,8:9,10,13)) {
  require(relevent)
  d <- mclapply(1:length(sims), function(i) {
    M <- nrow(sims[[i]]$edgelist)
    pshifts <- accum.ps(sims[[i]]$edgelist)[M,chosen] # 10 is ab-xy
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
