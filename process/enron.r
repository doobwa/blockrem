x <- read.csv("data/enron/enron-edgelist.csv")

# Get dyadic events among enron employees
tb <- table(x$mid)
ix <- which(tb > 1)
bad <- as.numeric(names(tb))[ix]
ix <- which(x$mid %in% bad)
y <- x[-ix,]

a <- as.matrix(table(y$s,y$r))
image(log(a))

A <- y[,c(3,4,5)]
N <- max(A[,2:3])
A[,1] <- A[,1] - A[1,1]
A[,1] <- A[,1]/A[nrow(A),1]

# Smaller
A <- y[1001:5000,c(3,4,5)]
A[,1] <- A[,1] - A[1,1]
A[,1] <- A[,1]/A[nrow(A),1] * 1000
rownames(A) <- c()
A <- as.matrix(A)
train.ix <- 1:3000
test.ix  <- 3001:4000
train <- A[train.ix,]
test  <- A[test.ix,]
save(A,N,train,test,train.ix,test.ix,file="data/enron-small.rdata")
