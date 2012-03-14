d <- read.table("data/eckmann/eckmann.txt",header=FALSE)
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
N <- length(unique(c(d$s,d$r)))
A <- d[,c("ntime","s","r")]
train <- A[1:20000,]
test  <- A[20001:nrow(A),]
save(A,train,test,N,file="data/eckmann.rdata")

# Create even smaller subset
ix <- which(A[,2] %in% 40:130 & A[,3] %in% 40:130)
A <- A[ix,]
A[,2] <- as.numeric(factor(A[,2],40:130))
A[,3] <- as.numeric(factor(A[,3],40:130))

save(A,file="data/eckmann/dyadic-small-sorted.rdata")

# Reorder names
B <- A
all <- sort(unique(c(A[,2],A[,3])))
all <- sample(all)
A[,2] <- match(A[,2],all)
A[,3] <- match(A[,3],all)
A[,1] <- 100 * A[,1] # rescale time
train <- as.matrix(A[1:2000,])
test.ix <- 2001:nrow(A)
test  <- as.matrix(A[test.ix,])
A <- as.matrix(A)
N <- length(all)
save(A,N,train,test,test.ix,file="data/eckmann-small.rdata")