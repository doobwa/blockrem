library(ggplot2)
source("R/brem.cpp.r")
source("R/brem.r")
source("R/sbm.r")

system("setenv OMP_NUM_THREADS 8")

load("data/eckmann/dyadic-small-sorted.rdata")

K <- 2
P <- 11
A <- as.matrix(A)
N <- max(unique(c(A[,2],A[,3])))
M <- nrow(A)
z <- c(rep(1,42),rep(2,N-42))  #sample(1:K,N,replace=TRUE)
beta <- array(rnorm(K^2*P),c(P,K,K))
beta[7:11,,] <- 0

niter <- 100
times <- A[,1]
sen <- A[,2]
rec <- A[,3]
s <- new(brem$Stat,times,sen-1,rec-1,N,M,P)
s$precompute()

niter <- 300
fit1 <- brem.mcmc(A,N,K,s,model.type="shared",niter=niter,gibbs="fast")
fit1.knownz <- brem.mcmc(A,N,K,s,model.type="shared",niter=niter,z=z,gibbs="none")
save(fit1,fit1.knownz,file="data/eckmann/fits.shared.rdata")

load("data/eckmann/fits.shared.rdata")
pdf("figs/eckmann/lpost.shared.pdf",width=5,height=5)
plot(fit1$llks,type="l",ylab="log posterior",xlab="iteration")
lines(fit1.knownz$llks,type="l",col="red")
dev.off()

ests <- list("gibbs"=fit1$beta,"known"=fit1.knownz$beta)
ds <- melt(ests)
qplot(factor(X1),value,data=ds,geom="point",colour=factor(L1)) + theme_bw() + labs(x="parameter",colour="block",y="estimate") + facet_grid(X2~X3)
ggsave("figs/eckmann/shared.pdf")

fit2 <- brem.mcmc(A,N,K,s,model.type="full",niter=niter,gibbs="fast")
fit2.knownz <- brem.mcmc(A,N,K,s,model.type="full",niter=niter,z=z,gibbs="none")
save(fit2,fit2.knownz,file="data/eckmann/fits.full.rdata")

fit3 <- brem.mcmc(A,N,1,s,model.type="full",niter=niter,gibbs="none")

# Visualize
plotspmat <- function(A) {
  mat <- table(A[,2],A[,3])
  mat <- melt(as.matrix(mat))
  mat <- subset(mat,value!=0)
  plot(mat[,1],mat[,2],pch=".")
}
plotspmat(A)


set.seed(4)
niter <- 100
px <- c(1,rep(0,6))
# fit0 <- brem.mcmc(A,N,K,P,px,model.type="baserates",niter=3)
z=c(rep(1,50),rep(2,38))
sbm.lpost(A,N,K,z,beta)

fit0 <- sbm.mcmc(A,N,2,niter=100)
fit1 <- sbm.mcmc(A,N,3,niter=100)
fit1 <- sbm.mcmc(A,N,10,niter=100)

# Fit with BREM
fit1 <- brem.mcmc(A,N,K,model.type="shared",niter=niter,gibbs="fast")

fit=fit1
K <- 2
nmap <- order(fit$z)
s <- match(A[,2],nmap)
r <- match(A[,3],nmap)
plot(s,r,pch=".",xlim=c(0,88),ylim=c(0,88))
cpoints <- sapply(1:K,function(k) min(which(sort(fit$z)==k)))
abline(v=cpoints[-1])
abline(h=cpoints[-1])
fit$beta

# Atlernate visualization
baseline <- mult.dir(cbind(s,r),c(N,N))
qplot(d1,d2,data=baseline,colour=count) + theme_bw()

# Fit rem to each separately
zs <- fit$z[A[,2]]
zr <- fit$z[A[,3]]
A_11 <- A[which(zs==1 & zr==1),]
A_22 <- A[which(zs==2 & zr==2),]
A_33 <- A[which(zs==3 & zr==3),]
N_2 <- sum(fit$z==2)
N_3 <- sum(fit$z==3)
px <- c(0,rep(1,6))
K <- 1
fit_22 <- brem.mcmc(A_22,N_2,K,P,px,model.type="full",niter=1000)
fit_33 <- brem.mcmc(A_33,N_3,K,P,px,model.type="full",niter=1000)
plot(fit_33$llks,type="l")

a <- melt(fit_33$param)
a$g <- "33"
b <- melt(fit_22$param)
b$g <- "22"
df <- rbind(a,b)
df$pshift <- c("baserate","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab")[df$X4]
df <- subset(df,X1>100)
qplot(X1, value, data=df,geom="line",colour=factor(pshift),group=factor(pshift)) + facet_grid(~g) + theme_bw() + labs(x="iteration",colour="pshift")

