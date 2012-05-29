## load("/extra/duboisc0/highschool2/working.data/setup.sm.rdata")
## save(sids,edgelists,N,Z,file="/extra/duboisc0/highschool2/working.data/classrooms.rdata")

load("data/classrooms.rdata")

Ms <- sapply(edgelists,nrow)
Ns <- N

tbs <- lapply(1:length(edgelists),function(i) {
  el <- edgelists[[i]]
  table(c(el[,2],el[,3]))
})

# Sort by entropy.  Don't want to model those dominated by a few people
h <- function(x) {
  p <- x/sum(x)
  mean(p * log(p))
}
hs <- sapply(tbs,h)

ix <- which(Ms > 400)
o <- order(hs[ix],decreasing=TRUE)
ix <- ix[o][1:5]

# Double check
hs[ix]  # smaller to larger
Ms[ix]  # over 400

for (i in ix) {
  A <- edgelists[[i]]
  N <- Ns[i]
  train.ix <- 1:300
  test.ix  <- 301:nrow(A)
  train <- A[train.ix,]
  test  <- A[test.ix,]
  save(A,N,train,test,train.ix,test.ix,file=paste("data/classroom-",i,".rdata",sep=""))
}
