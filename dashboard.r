#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
  make_option(c("-d", "--dataset"), 
              help="Name of dataset with data at /data/[dataset].rdata 
                    and results at /results/[dataset]/."),
  make_option(c("-p", "--predictions"), default=FALSE,
              help="compute prediction experiments?"),
  make_option(c("-b", "--save.dashboard"), default=FALSE,
              help=""),
  make_option(c("-s", "--save.figs"),default=FALSE,
              help="Save each individual figure to figs/[dataset]/ [default %default]."))
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opts   <- parse_args(OptionParser(option_list=option_list))

#' Perform analysis for a given dataset given the the current results of model fits.  
#' @dataset name of dataset.  Searches /results/[dataset]/ for results, and saves figures to /figs/[dataset]/
#' 
#' Each "res" object should have:
#' @z latent class assignments
#' @beta values of parameters in the last iteration
#' @llks log posterior at each MCMC iteration
#' @param array of parameter values at each MCMC iteration
#' @zs list of latent class assignment vectors at each MCMC iteration
#' The synthetic dataset should also have true.lpost variable available as well as the true beta parameters.
#' 
options(verbose=FALSE)
library(brem)
library(ggplot2)
chosen.model <- "3.FALSE.3"
opts <- list(dataset="synthetic-1",predictions=TRUE)
opts <- list(dataset="eckmann-small",predictions=TRUE)

# Pull data from each saved file and grab the name of the fit
results.dir <- paste("results/",opts$dataset,sep="")
load(paste("data/",opts$dataset,".rdata",sep=""))
dir.create(paste("results/",opts$dataset,"/ranks",sep=""),showWarn=FALSE)
dir.create(paste("results/",opts$dataset,"/llks",sep=""),showWarn=FALSE)
dir.create(paste("results/",opts$dataset,"/final",sep=""),showWarn=FALSE)

fs <- list.files(results.dir,full.names=TRUE)
fs <- fs[-grep("ranks",fs)]  # TODO: Bug if one of these not present
fs <- fs[-grep("final",fs)]
fs <- fs[-grep("llks",fs)]
ix <- grep("counts",fs)
if (length(ix) > 0) fs <- fs[-grep("counts",fs)]
fnames <- unlist(strsplit(fs,".rdata"))
fitnames <- sapply(fnames,function(f) strsplit(f,paste(results.dir,"/",sep=""))[[1]][2])
fits <- lapply(fs,function(f) {
  load(f)
  return(fit)
})
names(fits) <- fitnames

cat("Plotting log posterior during MCMC.\n")
llks <- lapply(1:length(fits),function(i) {
  data.frame(model=names(fits)[i],llk=fits[[i]]$lps,iter=1:fits[[i]]$niter)
})
llks <- do.call(rbind,llks)

llks <- subset(llks,iter>50)
if (opts$dataset == "synthetic") {
  load("data/synthetic.rdata")
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model))+ geom_abline(intercept=true.lpost,slope=0) + labs(x="iteration",y="log posterior",colour="model") + theme_bw() #+ limits(c(min(llks$llk),0),"y")
} else {
  q1 <- qplot(iter,llk,data=llks,geom="line",colour=factor(model)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()
}

cat("Trace plot of beta.\n")
library(coda)
load(paste("results/",opts$dataset,"/",chosen.model,".rdata",sep=""))
table(fit$params$z)
betas <- lapply(fit$samples,function(s) s$beta)
r <- melt(betas)

# Trace plot of intercept and pshifts
r <- subset(r, Var2 < 6 & Var3 < 6)
q1a <- qplot(L1,value,data=subset(r, L1 > 10 & Var1 %in% c(1,2,3,6)), colour=factor(Var1),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw() + facet_grid(Var2~Var3)#,scales="free")
q1a

# Trace plot of degree effects
q1b <- qplot(L1,value,data=subset(r, L1 > 10 & Var1 %in% c(8,10,12)), colour=factor(Var1),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw() + facet_grid(Var2~Var3,scales="free")
q1b

# Trace plot of all effects
q1b <- qplot(L1,value,data=subset(r, L1 > 40 & Var1 %in% fit$priors$px), colour=factor(Var1),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw() + facet_grid(Var2~Var3,scales="free")
q1b

# Trace plot of upper level means
mus <- lapply(fit$samples,function(s) {
  as.matrix(s$mu)
})
r <- melt(mus)
q1c <- qplot(L1,value,data=subset(r, Var1 %in% fit$priors$px), colour=factor(Var1),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw()# + facet_grid(Var2~Var3,scales="free")
q1c

# Trace plot of upper level sigmas
sigmas <- lapply(fit$samples,function(s) {
  as.matrix(s$sigma)
})
r <- melt(sigmas)
q1d <- qplot(L1,value,data=subset(r, Var1 %in% fit$priors$px), colour=factor(Var1),geom="line") + labs(colour="parameters for\n each block",x="iteration") + theme_bw()# + facet_grid(Var2~Var3,scales="free")
q1d

cat("Plot of posterior probabilities at each level of the hierarchy")
plot.posterior(train,N,fit)

cat("Plot overall rate\n")
if (FALSE) {
  s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,nrow(A),P)
  s$precompute()
  lrm <- LogIntensityArrayPc(fit$beta,fit$z-1,s$ptr(),K)
  for (m in 1:nrow(A)) diag(lrm[m,,]) <- -Inf
  total.lambda  <- sapply(1:nrow(A),function(i) sum(exp(lrm[i,,])))
  pdf("figs/online-debug-time1.pdf",width=7,height=7)
  par(mfrow=c(1,1))
  time.hat  <- cumsum(1/total.lambda)
  time.fixed <- cumsum(rep(A[nrow(A),1]/nrow(A),nrow(A)))
  time.obs <- A[,1]
  plot(time.obs,type="l")
  lines(time.hat,col="red")
  lines(time.fixed,col="green")
  dev.off()
}

#qplot(iter,log(-llk + 1),data=subset(llks,iter < 20 & iter > 15),geom="line",colour=factor(model)) + labs(x="iteration",y="log posterior",colour="model") + theme_bw()

cat("Plotting progress of z at each MCMC iteration.\n")
zs <- lapply(fits,function(f) {
  z <- lapply(f$samples,function(s) s$z)
  do.call(rbind,z)
})
names(zs) <- names(fits)
zs <- melt(zs)
q2 <- qplot(Var1,Var2,data=zs,geom="tile",fill=factor(value)) + facet_wrap(~L1) + labs(x="iteration",y="node")

cat("Compute distribution of activity using results from the full model.\n")
tb <- table(factor(c(train[,2],train[,3]),1:N))
fx <- grep(chosen.model,names(fits))
z <- fits[[fx]]$params$z
df <- data.frame(group=z,count=tb)
q3a <- qplot(log(count.Freq),data=df,geom="histogram")+facet_grid(group~.,scales="free")+labs(y="number of users",x="log(number of events)") + theme_bw()

cat("Compute dyad counts for each pshift using results from the true model.\n")
fx <- grep(chosen.model,fnames)
if (length(fx) > 0) {
  if (!opts$dataset %in% c("synthetic","synthetic-1")) {
    z <- fits[[fx]]$params$z
  } else  {
    z <- c(1,1,1,1,1,2,2,2,2,2)
  }
  
  df <- dyad.ps(train,N)
  df <- melt(df)
  df$i <- z[df$Var1]
  df$j <- z[df$Var2]
  df <- subset(df,Var3 != "AB-AB")
  df$Var3 <- factor(as.character(df$Var3))
  levels(df$Var3) <- rev(list("AB-BA"="ab-ba","AB-BY"="ab-by","AB-XA"="ab-xa","AB-XB"="ab-xb","AB-AY"="ab-ay"))#,"AB-AB"="ab-ab"))
  q3 <- qplot(Var3,value,data=df,geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="",y="count for a given dyad") + coord_flip()# +opts(axis.text.x=theme_text(angle=90)) 
  q3
  ggsave("figs/synthetic/counts.pdf",width=4,height=4)
} else {
  q3 <- qplot(0,0,label="not available",geom="text")
}


# Plot stats for particular example

  df <- dyad.ps(train,N)
dimnames(df)[[3]] <- c("ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab","total")
  df <- melt(df)
  df$i <- z[df$Var1]
  df$j <- z[df$Var2]
  df <- subset(df,Var3 %in% c("total","ab-ba","ab-ay","ab-by"))
  df$Var3 <- factor(as.character(df$Var3),c("total","ab-ba","ab-ay","ab-by"))
  df <- subset(df,(i==1 & j==3) | (i==3 & j==1))
#  levels(df$Var3) <- rev(list("AB-BA"="ab-ba","AB-BY"="ab-by","AB-XA"="ab-xa","AB-XB"="ab-xb","AB-AY"="ab-ay"))#,"AB-AB"="ab-ab"))
  q3 <- qplot(Var3,value,data=df,geom="boxplot",outlier.size=0.5) + facet_wrap( ~ i+j) + theme_bw() + labs(x="",y="count for a given dyad") + scale_y_continuous(limits=c(0,35))
  q3

effs <- c("intercept","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab")
b <- lapply(fit$samples,function(s) s$beta)
b <- melt(b)
b <- subset(b,Var1 %in% fit$priors$px & L1 > 20 &
            ((Var2 == 1 & Var3 == 3) | (Var2 ==3 & Var3 == 1)))
b$stat <- factor(effs[b$Var1],effs)

qplot(stat,value,data=b,geom="boxplot",outlier.size=0.5)+facet_wrap(~Var2 + Var3) + theme_bw() + xlab("parameter")



cat("Recall experiment on training data.\n")

load(paste("data/",opts$dataset,".rdata",sep=""))
folder <- paste("results/",opts$dataset,"/ranks/",sep="")
rks <- lapply(dir(folder,full.names=TRUE),function(f) {
  load(f)
  return(list(train=rk.train,test=rk.test))
})
fnames <- unlist(strsplit(dir(folder),".rdata"))
names(rks) <- fnames

ds <- lapply(rks,function(r) {
  d <- list(train=list(recall.200 = recall(r$train,top=1:100),
                       recall.all = recall(r$train,top=seq(1,N^2,length.out=100))),
            test =list(recall.200 = recall(r$test,top=1:100),
                       recall.all = recall(r$test,top=seq(1,N^2,length.out=100))))
  return(melt(d,id.vars="k",measure.vars="recall"))
})
for (i in 1:length(ds)) ds[[i]]$model <- fnames[i]
ds <- do.call(rbind,ds)
rownames(ds) <- c()

# Recall at k
q4 <- qplot(k,value,data=subset(ds,L2=="recall.200"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q4)
#ggsave(paste("figs/",opts$dataset,"/recall.200.pdf",sep=""),width=10,height=8)

q5 <- qplot(k,value,data=subset(ds,L2=="recall.all"),geom="line",colour=factor(model),group=factor(model)) + facet_grid(L1 ~ L2,scales="free") + theme_bw() + labs(x="cutpoint k",y="recall",colour="model")
print(q5)
#ggsave(paste("figs/",opts$dataset,"/recall.all.pdf",sep=""),width=10,height=8)

# # Experimental code for looking at the relation between the rank and the observed activity
# tbtrain <- table(factor(c(train[,2],train[,3]),1:N))
# tbtest <- table(factor(c(test[,2],test[,3]),1:N))
# rks.all <- lapply(1:length(rks),function(i) 
#   rbind(
#   data.frame(model=fnames[i],
#              iter=1:length(rks[[i]]$train),
#              rk=rks[[i]]$train,type="train",
#              tc=tbtrain[train[,2]] + tbtrain[train[,3]]),
#   data.frame(model=fnames[i],
#              iter=1:length(rks[[i]]$test),
#              rk=rks[[i]]$test,type="test",
#              tc=tbtest[test[,2]] + tbtest[test[,3]])
#   )
# )
# rks.all <- do.call(rbind,rks.all)
# qplot(iter,rk,data=subset(rks.all,type=="train"),colour=factor(model),alpha=.2)+facet_grid(model~.)+theme_bw()
# ggsave(paste("figs/",opts$dataset,"/rank.train.pdf",sep=""),width=20,height=12)
# qplot(iter,rk,data=subset(rks.all,type=="test"),colour=factor(model),alpha=.2)+facet_grid(model~.)+theme_bw()
# ggsave(paste("figs/",opts$dataset,"/rank.test.pdf",sep=""),width=20,height=12)
# qplot(tc,rk,data=rks.all,colour=factor(model),alpha=.1)+facet_grid(type ~ model)+theme_bw() + labs(x="sum of sender and recipient degree",y="rank of observed event")
# ggsave(paste("figs/",opts$dataset,"/rank.v.count.pdf",sep=""),width=20,height=12)

cat("Plotting parameter estimates.\n")
betas <- lapply(fits,function(f) f$params$beta)
names(betas) <- names(fits)
b <- melt(betas)
q6 <- qplot(Var1,value,data=b,geom="point",colour=factor(L1)) + facet_grid(Var2~Var3) + theme_bw() #+ labs(colour="model")
#qplot(X1,value,data=subset(b,L1=="full.3"),geom="point") + facet_grid(X2~X3) + theme_bw() #+ labs(colour="model")

# Get names
if (opts$dataset=="twitter-small") {
  lapply(1:3,function(k) nmap[which(z==k)])
}

# b <- melt(fits[['results/twitter-small/full.3']]$beta)
# b$block <- paste(b$X2,b$X3)
# qplot(factor(X1),value,data=b,geom="point",colour=factor(block))+theme_bw()+labs(x="parameter")

if (opts$predictions) {
  # # Compute out of sample log posterior
  load(paste("data/",opts$dataset,".rdata",sep=""))
  folder <- paste("results/",opts$dataset,"/llks/",sep="")
  fs <- dir(folder,full.names=TRUE)
  llks.test <- lapply(fs,function(f) {
    load(f)
    return(llk.test)
  })
  llks.train <- lapply(fs,function(f) {
    load(f)
    return(llk.train)
  })
  mllks.test <- lapply(fs,function(f) {
    load(f)
    return(llkm.test)
  })
  mllks.train <- lapply(fs,function(f) {
    load(f)
    return(llkm.train)
  })
  names(llks.test) <- names(llks.train) <- names(mllks.test) <- names(mllks.train) <- strsplit(dir(folder),".rdata")
  
  llks.train <- melt(llks.train)
  llks.train$event <- 1:nrow(train)
  llks.test <- melt(llks.test)
  llks.test$event <- 1:nrow(test)
  mllks.train <- melt(mllks.train)
  mllks.train$event <- 1:nrow(train)
  mllks.test <- melt(mllks.test)
  mllks.test$event <- 1:nrow(test)
  
  if (length(grep("full.2",fs) & length(grep("online",fs)) > 0)) {
    tmp <- subset(llks.test,L1 %in% c("full.2","online"))
    tmp <- cast(tmp, event ~ L1)
    q7 <- ggplot(tmp) + geom_point(aes(x=online,y=full.2),alpha=.1) + geom_abline(intercept=0,slope=1,colour="red") + theme_bw()
    pdf(paste("figs/",opts$dataset,"/lpost-compare-2-online.pdf",sep=""),height=4,width=4)
    print(q7)
    dev.off()
  }
  
  # Make results table
  df <- rbind(
        cbind(likelihood="rem",type="train",
              ddply(llks.train,.(L1),summarise,value=mean(value))),
        cbind(likelihood="rem",type="test",
              ddply(llks.test,.(L1),summarise,value=mean(value))),
        cbind(likelihood="mult",type="train",
              ddply(mllks.train,.(L1),summarise,value=mean(value))),
        cbind(likelihood="mult",type="test",
              ddply(mllks.test,.(L1),summarise,value=mean(value)))
              )
 # TODO: Fix this up 
#  df$L1 <- factor(df$L1,c("uniform","marg","online","full.1","full.2","full.3","dp","truth"))
  df$dataset <- opts$dataset
  save(df,file=paste("results/",opts$dataset,"/final/results.rdata",sep=""))
}

cat("Creating dashboard.\n")
if (save.dashboard) {
  library(gridExtra)
  pdf(paste("figs/",opts$dataset,"/dashboard.pdf",sep=""),width=25,height=10)
  blankPanel <- grid.rect(gp=gpar(col="white"))
  grid.arrange(q1, q1a, q2, q3, q3a, q4, q5, q6,ncol=4)
  dev.off()
}

if (opts$dataset=="eckmann-small") {
#n  load("data/eckmann-small.rdata")
  load(paste("results/eckmann-small/",chosen.model,".rdata",sep=""))
  b <- lapply(fit$samples,function(x) x$beta)
  b <- melt(b)
  colnames(b) <- c ("p","k1","k2","value","iter")
  b <- subset(b, value != 0)# & iter > 50)
  pnames <- c("intercept","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab","sen. outdegree","rec. outdegree","sen. indegree","rec. indegree","dyad count","total")
  b$p <- pnames[b$p]
  d <- ddply(b,.(p,k1,k2),function(x) c(mean=mean(x$value),quantile(x$value,c(.025,.2,.8,.975))))
  
  d$p <- factor(as.character(d$p),rev(pnames))

  pdf(paste("figs/",opts$dataset,"/params-estimates.pdf",sep=""),width=30,height=30)
  ggplot(d) + geom_point(aes(x=p,y=mean),colour="black") +  geom_linerange(aes(x=p,ymin=`20%`,ymax=`80%`),colour="black") + geom_linerange(aes(x=p,ymin=`2.5%`,ymax=`97.5%`),colour="black")  +facet_grid(k1~k2) + theme_bw() + labs(x="",y="value") + coord_flip()
  dev.off()
}

if (opts$dataset=="synthetic") {

  load("data/synthetic-1.rdata")
  load(paste("results/synthetic-1/",chosen.model,".rdata",sep=""))
  b <- lapply(fit$samples,function(x) x$beta)
  b <- melt(b)
  colnames(b) <- c ("p","k1","k2","value","iter")
  b <- subset(b, iter > 50 & value != 0)
  pnames <- c("intercept","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab","sen. outdegree","rec. outdegree","sen. indegree","rec. indegree","dyad count","total")
  dimnames(beta)[[1]] <- pnames
  tb <- melt(beta)
  colnames(b) <- c ("p","k1","k2","value","iter")
  colnames(tb) <- c("p","k1","k2","value")
 # tb <- subset(tb,!(p %in% c("ab-ab","total")))
  d <- ddply(b,.(p,k1,k2),function(x) c(mean=mean(x$value),quantile(x$value,c(.025,.2,.8,.975))))
  
  # Fix group assignments
  if (fit$params$z[1]==2) {
    d$k1 <- c(2,1)[d$k1]
    d$k2 <- c(2,1)[d$k2]
  }
  
  # Fix parameter names
  d$p <- pnames[d$p]
  ignore <- c("ab-ab","ab-xa","ab-xb","total","rec. indegree","sen. indegree","rec. outdegree","sen. outdegree","dyad count")
  d <- subset(d,! p %in% ignore)
  tb <- subset(tb,! p %in% ignore)
  d$p <- factor(as.character(d$p),rev(pnames[-c(7:13)]))
  tb$p <- factor(as.character(tb$p),rev(pnames[-c(7:13)]))
  
  ggplot(d) + geom_point(aes(x=p,y=mean),colour="white") +  geom_linerange(aes(x=p,ymin=`20%`,ymax=`80%`),colour="white") + geom_linerange(aes(x=p,ymin=`2.5%`,ymax=`97.5%`),colour="white") + geom_point(colour="red",aes(x=p,y=value),data=tb) +facet_grid(k1~k2) + theme_bw() + labs(x="",y="value") + coord_flip()
  ggsave("figs/synthetic/params-true.pdf",width=5,height=4)
  ggplot(d) + geom_point(aes(x=p,y=mean),colour="black") +  geom_linerange(aes(x=p,ymin=`20%`,ymax=`80%`),colour="black") + geom_linerange(aes(x=p,ymin=`2.5%`,ymax=`97.5%`),colour="black") + geom_point(colour="red",alpha=.5,aes(x=p,y=value),data=tb) +facet_grid(k1~k2) + theme_bw() + labs(x="",y="value") + coord_flip()
  ggsave("figs/synthetic/params-estimates.pdf",width=5,height=4)
  
}

cat("N x N plot of mean beta_p,z_i,z_j for a given p.\n")
load(paste("data/",opts$dataset,".rdata",sep=""))
load(paste("results/",opts$dataset,"/3.FALSE.3.rdata",sep=""))
plot(fit$lps[-1])

iters <- 20:45#5:20#40:80#f1fit$niter
mats <- vector("list",length=P)
for (p in 1:P) {
  mats[[p]] <- 0
}
for (i in iters) {
  s <- fit$samples[[i]]
  for (p in 1:P) {
    mats[[p]] <- mats[[p]] + s$beta[p,s$z,s$z]
  }
}
for (p in 1:P) {
  mats[[p]] <- mats[[p]]/length(iters)
}
beta.pm <- mats

## Get upper level mean and variance for each effect p
mus <- do.call(rbind,lapply(fit$samples[iters],function(s) s$mu))
mu.hat <- colMeans(mus)
sigmas <- do.call(rbind,lapply(fit$samples[iters],function(s) s$sigma))
sigma.hat <- colMeans(sigmas)

## Order rows according to z from last sample
z <- fit$params$z
o <- order(z)

dir.create(paste("figs/",opts$dataset,sep=""),showWarn=FALSE)
dir.create(paste("figs/",opts$dataset,"/parmat/",sep=""),showWarn=FALSE)
px <- which(sigma.hat!=0)
pdf(paste("figs/",opts$dataset,"/parmat/all.pdf",sep=""),height=20,width=20)

par(mfrow=c(3,3),mar=c(3,1,1,1))

# Plot observed data with cols nad rows sorted
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
image(log(tb[o,o]+1),xaxt="n",yaxt="n",col=grey.colors(100))
#image(mats[[1]][o,o],xaxt="n",yaxt="n",col=grey.colors(100))
#par(mar=c(0,0,0,0))
px <- fit$priors$px

# Plot blocks
plot(1,xlim=c(1,N),ylim=c(1,N))
zs <- sort(z)
for (i in 2:N) {
  if (zs[i] != zs[i-1]) {
    print(i)
    text(i,i,label=paste(zs[i],zs[i]))
    abline(v=i,h=i)
  }
}

# Plot matrix of each effect
for (p in px) {
  mat <- (mats[[p]] - mu.hat[p]) / sigma.hat[p]
  mat <- mat[o,o]
  #pdf(paste("figs/",opts$dataset,"/parmat/",p,".pdf",sep=""),height=4,width=4)
  image(mat,xaxt="n",yaxt="n",col=grey.colors(100),main=pnames[p])
  #dev.off()
}

dev.off()

# Plot pshifts
tb <- table(z[train[,2]],z[train[,3]])
df <- dyad.ps(train,N,ego=1)
df <- melt(df)
stats1 <- cast(subset(df,Var3=="AB-BA"), Var1 ~ Var2)

df$zi <- z[df$Var1]
df$zj <- z[df$Var2]
table(df$zi,df$zj)
df <- subset(df,Var3 != "AB-AB")
df$Var3 <- factor(as.character(df$Var3))
#levels(df$Var3) <- rev(list("AB-BA"="ab-ba","AB-BY"="ab-by","AB-XA"="ab-xa","AB-XB"="ab-xb","AB-AY"="ab-ay"))#,"AB-AB"="ab-ab"))
q3 <- qplot(Var3,value,data=subset(df,i < 4 & j < 4),geom="boxplot",outlier.size=0.1) + facet_grid(i ~ j) + theme_bw() + labs(x="",y="count for a given dyad") + coord_flip()# +opts(axis.text.x=theme_text(angle=90)) 
  q3

# Plot of each (i,j) of beta_p,k1,k2 from a given (k1,k2)
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
tb <- melt(tb)
colnames(tb) <- c("i","j","count")
tb$i <- as.numeric(as.character(tb$i))
tb$j <- as.numeric(as.character(tb$j))
tb$zi <- z[tb$i]
tb$zj <- z[tb$j]
tb$pair <- paste(tb$zi,tb$zj)
tb <- subset(tb,pair %in% c("1 3","3 1"))

tb$beta1 <- beta.pm[[1]][as.matrix(tb[,1:2])]
tb$beta2 <- beta.pm[[2]][as.matrix(tb[,1:2])]

qplot(beta1,count,data=tb) +facet_wrap(~pair)

# LATEST: Plot of beta_p,k1,k2 across samples for two blocks
# Counts
qplot(pair,count,data=tb,geom="boxplot")  # observed
b <- lapply(fit$samples,function(s) s$beta)
b <- melt(b)
b <- subset(b,Var1 %in% fit$priors$px & L1 > 20 &
            ((Var2 == 1 & Var3 == 3) | (Var2 ==3 & Var3 == 1)))
qplot(factor(Var1),value,data=b,geom="boxplot")+facet_wrap(~Var2 + Var3)

# Pshifts
tmp <- df
tmp$zi <- z[as.numeric(as.character(tmp$Var1))]
tmp$zj <- z[as.numeric(as.character(tmp$Var2))]
tmp$pair <- paste(tmp$zi,tmp$zj)
tmp <- subset(tmp,pair %in% c("1 3","3 1"))
qplot(pair,value,data=tmp,geom="boxplot")+facet_grid(~Var3)





#for (i in 1:length(k1)) {
  k <- k1[i]
  l <- k2[i]
  x$pair <- paste(x$zi,x$zj)
  x <- subset(x,pair %in% c("1 3","3 1") & L2 %in% fit$priors$px)
  qplot(pair,value,data=x,geom="boxplot",outlier.size=0)+scale_y_continuous(limits=c(-1,1))

  qplot(factor(j),value,data=x,geom="boxplot",outlier.size=0)


betas <- lapply(fit$samples[20:70],function(s) s$beta)
betas <- melt(betas)



k1 <- 1


k1 <- 
k2 <- 4
k3 <- 5
tb <- table(z[train[,2]],z[train[,3]])
num <- outer(table(z),table(z),"*")
lam <- tb/num
ik3 <- which(z[o] == k3)
counts <- colSums(tb[o,o][ik1,c(ik2,ik3)])

mats[[1]][ik1,ik2]
mats[[1]][ik1,k3]
mats[[2]][k1,k2]
mats[[2]][k1,k3]

k1 <- 2; k2 <- 5
beta[1,2,5]
beta[1,5,2]
beta[3,2,5]
beta[3,5,2]


## Find cases where beta_p1,k1,k2 > beta_p1,k1,k3 but beta_p2,k1,k2 < beta_p2,k1,k3.  Just use last draw
K <- dim(fit$params$beta)[2]
beta <- fit$params$beta
current <- list()
i <- 0
for (k1 in 1:K) {
  print(k1)
  for (k2 in (1:K)[-k1]) {
   for (k3 in (1:K)[-c(k1,k2)]) {
    for (p1 in 1) {
      for (p2 in px) {
        if (#abs(beta[p1,k1,k2] - beta[p1,k1,k3]) < sigma.hat[p1] &
            (beta[p1,k1,k2] - mu.hat[p1]) / sigma.hat[p1] > 2 &
            (beta[p1,k1,k3] - mu.hat[p1]) / sigma.hat[p1] > 2 &            
            ((beta[p2,k1,k2] - mu.hat[p2]) / sigma.hat[p2])*2 <
            ((beta[p2,k1,k3] - mu.hat[p2]) / sigma.hat[p2])) {
          current[[i <- i+1]] <- c(p1,p2,k1,k2,k3)
        }
      }
    }
   }
  }
}
spots <- do.call(rbind,current)
colnames(spots) <- c("p1","p2","k1","k2","k3")
spots <- data.frame(spots)
tmp <- subset(spots,p1 == 1 & p2 %in% c(3,6) &
              k1 %in% c(1:4,7,9,11,12,15,16) &
              k2 %in% c(1:4,7,9,11,12,15,16) &
              k3 %in% c(1:4,7,9,11,12,15,16))

res <- lapply(1:length(tmp),function(i) {
  profile(train,fit,tmp$p1[i],tmp$p2[i],tmp$k1[i],tmp$k2[i],tmp$k3[i])
})

profile <- function(train,fit,p1,p2,k1,k2,k3) {
  c(beta[p1,k1,k2], beta[p1,k1,k3], beta[p2,k1,k2], beta[p2,k1,k3])
  z <- fit$params$z
  ik1 <- which(z == k1)            # members of k2
  ik2 <- which(z == k2)            # members of k2
  ik3 <- which(z == k3)
  ix <- which(train[,2] %in% ik1 & train[,3] %in% ik2)
  iy <- which(train[,2] %in% ik1 & train[,3] %in% ik3)
  as.vector(c(table(train[ix,3]),table(train[iy,3])))
}

## Profile plot working off of posterior mean
tb <- table(z[train[,2]],z[train[,3]])
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
ik1 <- which(z == k1)            # members of k2
ik2 <- which(z == k2)            # members of k2
ik3 <- which(z == k3)
counts <- colSums(tb[ik1,c(ik2,ik3)])

## # Old version
## beta.ij <- lapply(fit$samples[iters],function(s) {
##   lapply(1:P,function(p) {
##     s$beta[p,s$z,s$z]
##   })
## })
## beta.ij <- melt(beta.ij)
## colnames(beta.ij) <- c("Var1","Var2","value","L2","L1")
## tmp <- ddply(beta.ij,.(Var1,Var2,L2),summarise,mean=mean(value))

## dir.create(paste("figs/",opts$dataset,sep=""),showWarn=FALSE)
## dir.create(paste("figs/",opts$dataset,"/parmat/",sep=""),showWarn=FALSE)
## px <- which(sigma.hat!=0)
## for (p in px) {
##   a <- subset(tmp,L2==p)
##   a <- data.frame(X1=a$Var1,X2=a$Var2,value=a$mean)
## #  plotmat(a,limits=c(-2,0))
##   mat <- as.matrix(cast(a,X1~X2)[,-1])
##   mat <- (mat - mu.hat[p]) / sigma.hat[p]
##   pdf(paste("figs/",opts$dataset,"/parmat/",p,".pdf",sep=""),height=4,width=4)
##   image(mat,xaxt="n",yaxt="n",col=grey.colors(100))#rainbow(100))
##   dev.off()
## }

cat("Saving figures\n")
if (opts$save.figs) {
  print(q1)
  ggsave(paste("figs/",opts$dataset,"/lpost.pdf",sep=""),height=4,width=5)
  print(q1a)
  ggsave(paste("figs/",opts$dataset,"/lpost-iter.pdf",sep=""),height=30,width=30)
  print(q2)
  ggsave(paste("figs/",opts$dataset,"/zs.pdf",sep=""),height=4,width=5)
  print(q3)
  ggsave(paste("figs/",opts$dataset,"/counts.pdf",sep=""),width=6,height=4)
  print(q3a)
  ggsave(paste("figs/",opts$dataset,"/countsa.pdf",sep=""),width=6,height=4)
  print(q4)
  ggsave(paste("figs/",opts$dataset,"/recall.pdf",sep=""),width=5,height=4)
  print(q5)
  ggsave(paste("figs/",opts$dataset,"/recall-all.pdf",sep=""),width=5,height=4)
  print(q6)
  ggsave(paste("figs/",opts$dataset,"/parameters.pdf",sep=""),width=5,height=4)
  print(q7)
  ggsave(paste("figs/",opts$dataset,"/traceplots.pdf",sep=""),width=5,height=4)
}

