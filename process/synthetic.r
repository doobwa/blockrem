
library(ggplot2)
library(brem)

N <- 10
K <- 2
beta <- list("intercept"=matrix(c(0,-1,-1,0),K,K),
             "abba" = matrix(c(3,0,0,0),K,K),
             "abby"=matrix(0,K,K),
             "abxa"=matrix(0,K,K),
             "abxb"=matrix(0,K,K),
             "abay"=matrix(c(0,0,0,2),K,K),
             "abab"=matrix(c(0,-1,-1,0),K,K),
             "sod"=matrix(0,K,K),
             "rod"=matrix(0,K,K),
             "sid"=matrix(c(0,0,0,0),K,K),
             "rid"=matrix(0,K,K),
#             "dc" =matrix(c(1,0,0,2),K,K),
             "dc" =matrix(0,K,K),
             "cc" =matrix(0,K,K))
P <- length(beta)
beta <- abind(beta,rev.along=3)

M <- 2000
set.seed(1)
z <- c(rep(1,5),rep(2,5))

sim <- simulate.brem(M,N,z,beta)

A <- train <- sim$A
plot(A[,1],type="l")

test <- simulate.brem(5000,N,z,beta)$A

# Compute log likelihood and log posterior of data under true model
times <- sim$A[,1]
sen <- sim$A[,2]
rec <- sim$A[,3]
s <- new(brem:::RemStat,times,sen-1,rec-1,N,M,P)
s$precompute()
llk.true <- loglikelihood_fast(beta,z-1,s$ptr(),K)  
true.lpost <- brem.lpost.fast(A,N,K,z,s,beta)

A <- rbind(train,test)
save(A,sim,N,K,P,M,z,beta,train,test,true.lpost,file="data/synthetic.rdata")


mat <- table(sim$A[,2],sim$A[,3])
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
edgelist <- melt(table(sim$A[,2],sim$A[,3]))
colnames(edgelist) <- c("sen","rec","val")
mypal <- brewer.pal(9,"Greys")
edge.colors <- col2rgb(mypal[as.numeric(cut(edgelist$val,9))])
edge.colors <- rbind(edge.colors,100)
edgelist <- as.matrix(edgelist)
attr(edgelist,"n") <- N
net <- as.edgelist.sna(edgelist)
pdf("figs/synthetic/network.pdf",width=20,height=20)
par(mar=rep(0,4))
set.seed(4)
gplot(net,pad=0,thresh=10,vertex.cex=2,label.cex=4,vertex.sides=30,label=1:10)#,edge.col=1:5)#t(edge.colors))
dev.off()

test_that("simulated lrm agrees with brem.lrm",{
  tmp <- brem.lrm(sim$A,N,z,beta)
  for (i in 1:M) diag(tmp[i,,]) <- -Inf
  diag(sim$lrm[1,,]) <- -Inf
  expect_that(all.equal(tmp[-1,,],sim$lrm[-1,,]),is_true())
})





## OLD CODE
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
