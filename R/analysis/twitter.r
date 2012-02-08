source("R/brem.r")
source("R/brem.cpp.r")
load("data/twitter.example.rdata")
s <- new(brem$Stat,A[,1],A[,2]-1,A[,3]-1,N,M,P)
s$precompute()

K <- 2
niter <- 500
#fit0 <- sbm.mcmc(A,N,K,niter=niter,gibbs=TRUE)
fit1 <- brem.mcmc(A,N,K,s,model.type="shared",niter=niter,gibbs="fast")
fit2 <- brem.mcmc(A,N,K,s,model.type="full",niter=niter,gibbs="fast")
fit3 <- brem.mcmc(A,N,1,s,model.type="full",niter=niter,gibbs="fast",mcmc.sd=.05)
save(fit0,fit1,fit2,fit3,file="data/twitter/fits.rdata")

z <- sample(1:K,N,replace=TRUE)
beta <- matrix(rnorm(K^2),K,K)

set.seed(4)
niter <- 100
px <- c(1,rep(0,6))
sbm.lpost(B,N,K,z,beta)

fit0 <- sbm.mcmc(B,N,2,niter=100)
fit1 <- sbm.mcmc(B,N,3,niter=100)
fit1 <- sbm.mcmc(B,N,10,niter=100)

fit <- fit0
nmap <- order(fit$z)
s <- match(B[,2],nmap)
r <- match(B[,3],nmap)
pdf("figs/twitter/sorted.pdf",width=4,height=4)
plot(s,r,pch=".",xlab="sender",ylab="recipient")
cpoints <- sapply(1:K,function(k) min(which(sort(fit$z)==k)))
abline(v=cpoints[-1])
abline(h=cpoints[-1])
dev.off()
fit$beta
