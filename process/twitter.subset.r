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


# Only chose events involving people with 2 or more evnets
users <- c(as.character(A$s),as.character(A$r))
users <- factor(users,unique(users))
tb <- table(users)
table(tb)
ix <- which(tb > 1)
chosen <- names(tb)[ix]
ix <- which(A$s %in% chosen & A$r %in% chosen)
A <- A[ix,]
N <- length(chosen)
bad <- which(as.character(A$s) == as.character(A$r))
A <- A[-bad,]
M <- nrow(A)

# Rescale time to be in (0,1)
A[,1] <- as.numeric(A[,1])
A[,1] <- A[,1] - A[1,1]
times <- A[,1]/A[nrow(A),1]

nmap <- sort(chosen)
sen <- match(A[,2],nmap)
rec <- match(A[,3],nmap)
A <- cbind(times,sen,rec)
N <- length(nmap)
train.ix <- 1:2000
test.ix <- 2001:nrow(A)
train <- A[train.ix,]
test  <- A[test.ix,]

save(A,train.ix,test.ix,train,test,N,M,nmap,file="data/twitter-small.rdata")
