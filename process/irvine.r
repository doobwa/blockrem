x <- read.table("data/irvine/contactData.txt")
colnames(x) <- c("s","r","day")

plot(x[,3],type="l")

tb <- table(c(x[,1],x[,2]))
plot(table(table(c(x[,1],x[,2]))))

a <- unique(c(x[,1],x[,2]))
chosen <- as.numeric(names(tb)[which(tb > 30)])

A <- x[,c(3,1,2)]
ix <- which(A[,2] %in% chosen & A[,3] %in% chosen)
A <- A[ix,]
N <- length(chosen)
s <- as.numeric(factor(A[,2],chosen))
r <- as.numeric(factor(A[,3],chosen))
plot(s,r,pch=".")
dim(A)

train <- A[1:5000,]
test  <- A[5001:nrow(A),]
save(A,train,test,N,file="data/irvine.rdata")