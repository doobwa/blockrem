library(hash)

ts <- sort(runif(10000,0,100))
x <- rnorm(10000)
a <- list()
for (i in 1:length(ts)) {
  a[as.character(ts[i])] <- x[i]
}
b <- hash(key=ts,values=x)

head(a)
system.time(a[8000:9000])
system.time(b[ts[8000:9000]])