
library(brem)
library(reshape2)
library(ggplot2)

load(paste("data/eckmann-small.rdata",sep=""))
load(paste("figs/eckmann-sm-fit2.rdata",sep=""))

iters <- 50:340

# Plot observed stats for total and a few pshifts for two blocks
df <- dyad.ps(train,N)
dimnames(df)[[3]] <- c("ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab","total")
df <- melt(df)
z <- fit$params$z  # use latest clustering
df$i <- z[df$Var1]
df$j <- z[df$Var2]
df <- subset(df,Var3 %in% c("total","ab-ba","ab-by"))
df$Var3 <- factor(as.character(df$Var3),c("total","ab-ba","ab-by"))
df <- subset(df,(i==1 & j==3) | (i==3 & j==1))
q3 <- qplot(Var3,value,data=df,geom="boxplot",outlier.size=0.5) + facet_wrap( ~ i+j) + theme_bw() + labs(x="",y="count for a given dyad") + scale_y_continuous(limits=c(0,35))

pdf("figs/eckmann-small/example-obs-stats.pdf",width=4,height=3)
q3
dev.off()

# Boxplot of posterior samples for each of the above blocks
effs <- c("intercept","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab")
b <- lapply(fit$samples,function(s) s$beta)
b <- melt(b)
b <- subset(b,Var1 %in% c(1,2,3) & L1 > 20 &
            ((Var2 == 1 & Var3 == 3) | (Var2 ==3 & Var3 == 1)))
b$stat <- factor(effs[b$Var1],c("intercept","ab-ba","ab-by"))

pdf("figs/eckmann-small/example-estimates.pdf",width=4,height=3)
qplot(stat,value,data=b,geom="boxplot",outlier.size=0.5)+facet_wrap(~Var2 + Var3) + theme_bw() + xlab("parameter")
dev.off()

# Get posterior mean of beta
P <- dim(fit$params$beta)[1]
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

par(mfrow=c(3,3),mar=c(3,1,1,1))

# Plot observed data with cols nad rows sorted
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
image(log(tb[o,o]+1),xaxt="n",yaxt="n",col=rev(grey.colors(100)))
#image(mats[[1]][o,o],xaxt="n",yaxt="n",col=grey.colors(100))
#par(mar=c(0,0,0,0))
px <- fit$priors$px

# Plot blocks
plot(1,xlim=c(1,N),ylim=c(1,N),xaxt="n",yaxt="n")
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
  image(mat,xaxt="n",yaxt="n",col=rev(grey.colors(100)),main=effs[p])
}

dev.off()


# Plot observed data with cols nad rows sorted
cols <- rev(grey.colors(100,start=0,end=.95,gamma=1))
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
pdf("figs/eckmann-small/parmat/observed.pdf",height=4,width=4)
par(mar=c(0,0,0,0))
image(log(tb[o,o]+1),xaxt="n",yaxt="n",col=cols)
dev.off()

# Plot matrix of each effect
for (p in px) {
  mat <- (mats[[p]] - mu.hat[p]) / sigma.hat[p]
  mat <- mat[o,o]
  pdf(paste("figs/eckmann-small/parmat/",p,".pdf",sep=""),height=4,width=4)
  par(mar=c(0,0,0,0))
  image(mat,xaxt="n",yaxt="n",col=cols)
  dev.off()
}

# Zoom in plots
ix <- o[which(z[o] %in% c(1,3))]  # members of 1 or 3
tb <- table(factor(train[,2],1:N),factor(train[,3],1:N))
pdf("figs/eckmann-small/parmat/observed-zoom.pdf",height=4,width=4)
par(mar=c(0,0,0,0))
image(log(tb[ix,ix]+1),xaxt="n",yaxt="n",col=cols)
dev.off()
for (p in px) {
  mat <- (mats[[p]] - mu.hat[p]) / sigma.hat[p]
  pdf(paste("figs/eckmann-small/parmat/",p,"-zoom.pdf",sep=""),height=4,width=4)
  par(mar=c(0,0,0,0))
  image(mat[ix,ix],xaxt="n",yaxt="n",col=cols)
  dev.off()
}
