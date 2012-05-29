x <- read.csv("data/realitymining/voice-all.edges",sep=";",header=FALSE)
#unix_timestamp;caller_id;receiver_id;duration_in_second;communication_type
colnames(x) <- c("time","sen","rec","dur","type")
a <- read.csv('data/realitymining/all.attributes.txt',sep=";",header=FALSE)
colnames(a) <- c("id","core")

nodes <- a$id[which(a$core == 1)]
A <- x[which(x$sen %in% nodes & x$rec %in% nodes),]

bad <- which(A$sen == A$rec)
A <- A[-bad,]

A$sctime <- A$time - A$time[1]
A$sctime <- A$sctime/max(A$sctime) * 1000
A$sen <- A$sen + 1
A$rec <- A$rec + 1
A <- as.matrix(A[,c("sctime","sen","rec")])
N <- length(nodes)
train.ix <- 1:1500
test.ix <- 1501:nrow(A)
train <- A[train.ix,]
test  <- A[test.ix,]
save(A,N,train,test,test.ix,file="data/realitymining-sm.rdata")


# Exploratory
A$datetime <- as.POSIXlt(A$time, origin="1970-01-01", tz="America/New_York")
dt <- strptime(A$datetime, "%Y-%m-%d %H:%M:%S", tz="UTC")
datetimes <- as.data.frame(lapply(dt,unclass))
datetimes$year <- 1900 + datetimes$year

plot(datetimes$yday[1:2000],type="l")

hist(datetimes$hour + datetimes$min/60,breaks=100)
hist(datetimes$wday)

md <- paste(datetimes$mon,datetimes$mday)
tmp <- table(md)
