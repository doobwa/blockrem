d <- read.table("eckmann/eckmann.txt",header=FALSE)
colnames(d) <- c("from","to","size","timestamp","chron",
                 "weekday","multiplicity")
d <- subset(d,multiplicity==1)
times <- d$chron - d$chron[1]
d$ntime <- times/max(times)
all <- sort(unique(c(d$from,d$to)))
min.count <- 50
max.count <- 1000
counts <- table(c(d$from,d$to))
ix <- which(counts > min.count & counts < max.count)
length(ix)
keep <- as.numeric(names(counts[ix]))
ix <- which(d$from %in% keep & d$to %in% keep)
d <- d[ix,c("ntime","from","to")]
              plot(as.numeric(table(table(c(d$from,d$to)))),type="l")

nodes <- sort(unique(c(d$from,d$to)))
d$s <- as.numeric(factor(d$from,nodes))
d$r <- as.numeric(factor(d$to,nodes))
         
A <- d[,c("ntime","s","r")]

# Visualize
plotspmat <- function(A) {
  mat <- table(A[,2],A[,3])
  mat <- melt(as.matrix(mat))
  mat <- subset(mat,value!=0)
  plot(mat[,1],mat[,2],pch=".")
}
plotspmat(A)

# Create even smaller subset
ix <- which(A[,2] %in% 40:130 & A[,3] %in% 40:130)
A <- A[ix,]
A[,2] <- as.numeric(factor(A[,2],40:130))
A[,3] <- as.numeric(factor(A[,3],40:130))

source("rem.cpp.r")
source("fns.r")
N <- length(unique(c(A[,2],A[,3])))
K <- 2
P <- 7
M <- nrow(A)

z <- sample(1:K,N,replace=TRUE)
beta <- array(rnorm(K^2*P),c(K,K,P))
system.time(brem.llk(A,N,z,beta,px))

set.seed(4)
niter <- 100
px <- c(1,rep(0,6))
fit0 <- brem.mcmc(A,N,K,P,px,model.type="baserates",niter=3)
