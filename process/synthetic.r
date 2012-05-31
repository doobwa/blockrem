
library(ggplot2)
library(brem)
source("pkg/R/brem.r")

N <- 10
K <- 2
beta <- list("intercept"=matrix(c(0,-1,-1,0),K,K),
             "abba" = matrix(c(3,0,0,0),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(c(0,0,0,2),K,K),
             "abab"=matrix(c(0,0,0,0),K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(c(0,0,0,0),K,K),
             "rid"=matrix(0,K,K),
#             "dc" =matrix(c(1,0,0,2),K,K),
             "dc" =matrix(0,K,K),
             "cc" =matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)

M <- 7000
set.seed(1)
z <- c(rep(1,5),rep(2,5))

sim <- generate.brem(M,N,beta,z,ego=TRUE)

A <- sim$edgelist
plot(A[,1],type="l")

train <- sim$edgelist[1:2000,]
test.ix <- 2001:7000
test <- sim$edgelist[test.ix,]
test[,1] <- test[,1] - test[1,1]

# Compute log likelihood and log posterior of data under true model
ego <- 1
times <- train[,1]
sen <- train[,2]
rec <- train[,3]
s <- new(RemStat,times,sen-1,rec-1,N,2000,P,ego)
s$precompute()
true.llk <- sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
priors <- list(beta=list(mu=0,sigma=1))
true.lpost <- brem.lpost.fast(train,N,K,z,s,beta,priors=priors)

#expect_that(true.lpost,equals(brem.lpost(train,N,K,z,beta,priors=priors)))
# 
# tmp1 <- RemLogLikelihoodPc(beta, z - 1, s$ptr(), K)
# tmp2 <- RemLogLikelihood(beta, times,sen-1,rec-1,z-1,N,nrow(train),K,P)

#TODO: have test.ix and A's time from 0, t_M
save(A,sim,N,K,P,M,z,beta,train,test,test.ix,true.lpost,file=paste("data/synthetic-",ego,".rdata",sep=""))


mat <- table(sim$edgelist[,2],sim$edgelist[,3])
mat <- melt(as.matrix(mat))
colnames(mat) <- c("X1","X2","value")
mat$X1 <- as.character(mat$X1)
mat$X1 <- factor(mat$X1,10:1)
mat$X2 <- factor(mat$X2,1:10)
plotmat(mat)
ggsave("figs/synthetic/mat.pdf",width=3,height=3)


set.seed(2)
ix <- sample(1:10)
mat$X1 <- factor(as.character(mat$X1),rev(ix))
mat$X2 <- factor(as.character(mat$X2),ix)
plotmat(mat)
ggsave("figs/synthetic/unsorted.pdf",width=3,height=3)

bm <- matrix(0,N,N)
bm[1:5,6:10] <- bm[6:10,1:5] <- -8
diag(bm) <- -10
mat <- melt(as.matrix(bm))
colnames(mat) <- c("X1","X2","value")
mat$X1 <- factor(as.character(mat$X1),10:1)
mat$X2 <- factor(mat$X2,1:10)
plotmat(mat,limits=c(-10,3))
ggsave("figs/synthetic/bm.pdf",width=3,height=3)

library(sna)
library(RColorBrewer)
edgelist <- melt(table(sim$edgelist[,2],sim$edgelist[,3]))
colnames(edgelist) <- c("sen","rec","val")
edgelist$sen <- as.numeric(as.character(edgelist$sen))
edgelist$rec <- as.numeric(as.character(edgelist$rec))
mypal <- brewer.pal(9,"Greys")
edge.colors <- col2rgb(mypal[as.numeric(cut(edgelist$val,9))])
edge.colors <- rbind(edge.colors,100)
edgelist <- as.matrix(edgelist)
attr(edgelist,"n") <- N
net <- as.edgelist.sna(edgelist)

pdf("figs/synthetic/network.pdf",width=20,height=20)
par(mar=rep(0,4))
set.seed(5)
gplot(net,pad=0,thresh=30,vertex.cex=2,label.cex=4,vertex.sides=30,label=1:10)#,edge.col=1:5)#t(edge.colors))
dev.off()

## OLD CODE


## test_that("simulated lrm agrees with brem.lrm",{
##   tmp <- brem.lrm(sim$A,N,z,beta)
##   for (i in 1:M) diag(tmp[i,,]) <- -Inf
##   diag(sim$lrm[1,,]) <- -Inf
##   expect_that(all.equal(tmp[-1,,],sim$lrm[-1,,]),is_true())
## })

# Look at distance between true and estimated parameter vectors
fit <- fit2
dist <- function(x,y) sqrt(sum((x-y)^2))
ds <- sapply(1:niter,function(i) {
  c(dist(beta[2:6,1,1],fit$param[i,1,1,2:6]),# - beta[1,1,1]),
    dist(beta[2:6,2,2],fit$param[i,2,2,2:6]),#, - beta[1,2,2]),
    dist(beta[2:6,2,1],fit$param[i,2,1,2:6]),#, - beta[1,2,1]),
    dist(beta[2:6,1,2],fit$param[i,1,2,2:6]))# - beta[1,1,2]))
})
ds <- melt(t(ds))
qplot(X1,value,data=ds,geom="line",colour=factor(X2)) + theme_bw() + labs(x="iteration",colour="block",y="Euclidean distance to truth")
ggsave("figs/syn/bias.pdf",width=5,height=4)


## Alternative example
N <- 10
K <- 2
beta <- list("intercept"=matrix(c(1.5,0,0,2),K,K),
             "abba" = matrix(c(2,0,.5,2),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(0,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K),
             "dc" =matrix(0,K,K),
             "cc" =matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)
M <- 1000
set.seed(1)
z <- c(rep(1,5),rep(2,5))
sim <- generate.brem(M,N,beta,z,ego=TRUE)
y1 <- sim$edgelist

beta <- list("intercept"=matrix(0,K,K),
             "abba" = matrix(c(2,0,.5,2),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(0,K,K),
             "abab"=matrix(0,K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(0,K,K),
             "rid"=matrix(0,K,K),
             "dc" =matrix(0,K,K),
             "cc" =matrix(0,K,K))
beta <- abind(beta,rev.along=3)
set.seed(1)
sim <- generate.brem(M,N,beta,z,ego=TRUE)
y2 <- sim$edgelist

cols <- rev(grey.colors(100,start=0,end=.95,gamma=1))
tb <- table(factor(y1[,2],1:N),factor(y1[,3],1:N))
pdf("figs/introexample/mat1.pdf",width=3,height=3)
par(mar=c(2,2,2,2),mgp=c(1.15,0,0))# oma=c(0,0,0,0),
tb <- as.matrix(tb)
image(1:N,1:N,tb,xaxt="n",yaxt="n",col=cols,xlab="",ylab="Observed counts",axes=FALSE,cex=1.25)
axis(1,at=c(3,8),labels=c("A","B"),tick=FALSE)
axis(2,at=c(3,8),labels=c("A","B"),tick=FALSE)
dev.off()

tb <- table(factor(y2[,2],1:N),factor(y2[,3],1:N))
pdf("figs/introexample/mat2.pdf",width=3,height=3)
par(mar=c(2,2,2,2),mgp=c(1.15,0,0))
image(1:N,1:N,tb,xaxt="n",yaxt="n",xlab="",ylab="",col=cols,axes=FALSE,cex=1.25)
axis(1,at=c(3,8),labels=c("A","B"),tick=FALSE)
axis(2,at=c(3,8),labels=c("A","B"),tick=FALSE)
dev.off()

pdf("figs/introexample/dynamics1.pdf",width=3,height=3)
par(mar=c(3,3,0,0),mgp=c(2,1.15,0))
x1 <- c(2,0,.5,2)
plot(1:4,x1,xlim=c(0.5,4.5),ylim=c(-1,3),xlab="",xaxt="n",ylab="Propensity for reciprocal events",pch=19,col="red",cex=1.25)
axis(1,at=1:4,labels=c("(A,A)","(A,B)","(B,A)","(B,B)"),tick=TRUE)
dev.off()

pdf("figs/introexample/dynamics2.pdf",width=3,height=3)
par(mar=c(3,3,0,0),mgp=c(2,1.15,0))
plot(1:4,x1,xlim=c(0.5,4.5),ylim=c(-1,3),xlab="Block",xaxt="n",ylab="",pch=19,col="red",cex=1.25)
axis(1,at=1:4,labels=c("(A,A)","(A,B)","(B,A)","(B,B)"),tick=TRUE)
dev.off()

pdf("figs/introexample/all.pdf",width=5,height=5)
par(mfrow=c(2,2))
cols <- rev(grey.colors(100,start=0,end=.95,gamma=1))
tb <- table(factor(y1[,2],1:N),factor(y1[,3],1:N))
par(mar=c(1,3,2,.5),mgp=c(1.15,0,0))# oma=c(0,0,0,0),
tb <- as.matrix(tb)
image(1:N,1:N,tb,xaxt="n",yaxt="n",col=cols,xlab="",ylab="Observed counts",axes=FALSE,cex=1.25)
axis(1,at=c(3,8),labels=c("A","B"),tick=FALSE)
axis(2,at=c(3,8),labels=c("A","B"),tick=FALSE)

tb <- table(factor(y2[,2],1:N),factor(y2[,3],1:N))
par(mar=c(1,3,2,.5),mgp=c(1.15,0,0))
image(1:N,1:N,tb,xaxt="n",yaxt="n",xlab="",ylab="",col=cols,axes=FALSE,cex=1.25)
axis(1,at=c(3,8),labels=c("A","B"),tick=FALSE)
axis(2,at=c(3,8),labels=c("A","B"),tick=FALSE)

par(mar=c(3,3,2,.5),mgp=c(2,1.15,0))
x1 <- c(2,0,.5,2)
plot(1:4,x1,xlim=c(0.5,4.5),ylim=c(-1,3),xlab="",xaxt="n",ylab="Propensity for reciprocal events",pch=19,col="red",cex=1.25)
axis(1,at=1:4,labels=c("(A,A)","(A,B)","(B,A)","(B,B)"),tick=TRUE)

plot(1:4,x1,xlim=c(0.5,4.5),ylim=c(-1,3),xlab="",xaxt="n",ylab="",pch=19,col="red",cex=1.25)
axis(1,at=1:4,labels=c("(A,A)","(A,B)","(B,A)","(B,B)"),tick=TRUE)

dev.off()
