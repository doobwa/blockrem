
# Create unique ids
load("data/rstats/rstats.interaction.rdata")
users <- sort(unique(c(as.character(A$s),as.character(A$r))))
A$sid <- match(as.character(A$s),users)
A$rid <- match(as.character(A$r),users)

# Visualize adjacency matrix
pdf("figs/twitter/adjmat.pdf",width=4,height=4)
plot(A[,4:5],pch=".",xlab="sender",ylab="recipient")
dev.off()

# Fit baseline model
N <- length(users)
K <- 2
P <- 7
M <- nrow(A)

# Rescale time to be in (0,1)
B <- A[,c(1,4,5)]
B[,1] <- as.numeric(B[,1])
B[,1] <- B[,1] - B[1,1]
B[,1] <- B[,1]/B[nrow(B),1]

z <- sample(1:K,N,replace=TRUE)
save(B,N,M,file="data/twitter.example.rdata")
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
