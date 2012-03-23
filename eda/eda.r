load("data/eckmann-small.rdata")
P <- 13
library(brem)
s <- new(RemStat,train[,1],train[,2]-1,train[,3]-1,N,nrow(train),P)
s$precompute()
s$transform()

# Arrange statistics vectors in large matrix with P columns
x <- s$get_all_s()
x[[1]][[2]][1:5]
xs <- list()
k <- 0
for (i in 1:N) {
  for (j in (1:N)[-i]) {
    xsk <- do.call(rbind,x[[i]][[j]][-1])
    ms <- xsk[,ncol(xsk)] + 1
    timegap <- c(0,diff(train[ms,1]))
    xs[[k <- k+1]] <- cbind(i,j,1:nrow(xsk),timegap,xsk)
  }
}
ys <- do.call(rbind,xs)
colnames(ys) <- c("i","j","v","timegap",paste("s",0:12,sep=""))

# Whether or not that dyad was observed
ms <- ys[,ncol(ys)]   # index of changepoints
obs <- (train[ms+1,2] == ys[,1]) & (train[ms+1,3] == ys[,2])
ys <- cbind(obs,ys)

# Sanity check
ix <- head(which(obs))
ys[ix,1:5]
train[ms[ix]+1,]

# Remove event ids
tmp <- ys[,-c(2,3,4,ncol(ys))]

# Kmeans
K <- 10
km <- kmeans(tmp,K)

clusters <- factor(km$cluster,1:K)
sen <- factor(ys[,2],1:N)
rec <- factor(ys[,3],1:N)
plot.phi(table(clusters,sen) + table(clusters,rec))
plot.phi(table(clusters,sen))
plot.phi(table(clusters,rec))

# Starting with most popular cluster, get list of actors with higher than average assignments.
plot.phi <- function(phi) {
    o.k <- order(rowSums(phi),decreasing=TRUE)
    ixs <- lapply(1:K,function(i) {
      o.a <- order(phi[o.k[i],],decreasing=TRUE)
      ix <- which(phi[o.k[i],o.a] > mean(phi[o.k[i],]))
      o.a[ix]
    })
    o <- unique(unlist(ixs))

    phi <- apply(phi,1,function(p) p/colSums(phi))

    library(ggplot2)
    phi <- melt(phi)
    colnames(phi)[1] <- "actor"
    phi$actor <- factor(phi$actor,o)
    phi$clusters <- factor(phi$clusters,o.k)
    q <- ggplot(phi, aes(x = clusters, y = actor)) + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "black")
    return(q)
}

