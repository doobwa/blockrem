load("data/eckmann-small.rdata")
P <- 13
library(brem)
s <- new(RemStat,train[,1],train[,2]-1,train[,2]-1,N,nrow(train),P)
s$precompute()
x <- s$get_all_s()
x[[1]][[2]][1:5]

s$get_s(99,0,1)
s$get_s(100,0,1)
s$get_s(101,0,1)

s$transform()

x <- s$get_all_s()
x[[1]][[2]][1:5]

library(reshape)
y <- melt(head(x))

sij <- s$get_s(100,0,1)

xs <- list()
k <- 0
for (i in 1:N) {
  for (j in (1:N)[-i]) {
    xs[[k <- k+1]] <- do.call(rbind,x[[i]][[j]])
    xs[[k]] <- cbind(i,j,1:nrow(xs[[k]]),xs[[k]])
  }
}

ys <- do.call(rbind,xs)
colnames(ys) <- c("i","j","v",paste("s",0:12))
