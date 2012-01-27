# Grab all #rstats tweets between 5/11/2009 and 1/26/2012
# Select tweets beginning with @
# Use first mention as recipient
# 28337 tweets in total
# 3926 directed among 1079 users

require(multicore)
require(rjson)
require(twitteR)
require(ROAuth)
require(RCurl)
source("pull.fns.r")
options(cores=4)
load("heroic.oauth.rdata")

# xml files obtained from google reader by hand by motifying the continuation code repeatedly
# http://www.google.com/reader/atom/feed%2Fhttp%3A%2F%2Fsearch.twitter.com%2Fsearch.atom%3Fq%3D%2523rstats?n=1000&c=CPHpgbSjipwC
# reference: http://code.google.com/p/pyrfeed/wiki/GoogleReaderAPI

# Xpath code found here: http://heuristically.wordpress.com/2011/04/08/text-data-mining-twitter-r/

grab.tweets.from.xml <- function(f) {
  require(XML)
  doc <- xmlParseDoc(f)
  urls <- xpathSApply(doc, '//s:link[@href]', xmlGetAttr,"href",
                      namespaces =c('s'='http://www.w3.org/2005/Atom'))
  urls <- urls[-1]  # remove google reader url
  titles <- xpathSApply(doc, '//s:title', xmlValue,
                        namespaces =c('s'='http://www.w3.org/2005/Atom'))
  ix <- grep("search.twitter.com",urls)
  urls <- urls[-ix]
  titles <- titles[-ix]
  dates <- xpathSApply(doc, '//s:published', xmlValue,
                       namespaces =c('s'='http://www.w3.org/2005/Atom'))
  tweets <- gsub("http://twitter.com/","",urls)
  users <- sapply(tweets,function(tweet) {
    strsplit(tweet,"/")[[1]][1]
  })
  names(users) <- c()
  df <- data.frame(date=dates,user=users,text=titles)
  return(df)
}

fs <- paste("data/rstats.",0:33,".xml",sep="")
df <- lapply(fs,function(f) grab.tweets.from.xml(f))
df <- do.call(rbind,df)
save(df,file="data/rstats.tweets.rdata")

# Select users and download their timelines
load("data/rstats.tweets.rdata")
K <- 1000
tb <- sort(table(df$user),decreasing=TRUE)[1:K]
plot(log(tb),type="l")
chosen <- names(tb)

setup <- expand.grid(page=1:16,user=chosen)
urls <- sapply(1:nrow(setup), function(i) {
    paste("https://api.twitter.com/1/statuses/user_timeline.json?include_entities=true&include_rts=true&count=200&screen_name=",setup$user[i],"&page=",setup$page[i],sep="")
})
save(urls,file="data/rstats.urls.rdata")


# Try grabbing directed messages among #rstats tweets
load("data/rstats.tweets.rdata")
get.mentions <- function(tweet) {
  w <- strsplit(as.character(tweet)," ")[[1]]
  fl <- sapply(w,substr,0,1)
  ix <- which(fl=="@")
  mentions <- w[ix]
  mentions <- gsub("@","",mentions)
  mentions <- gsub(":","",mentions)
  return(mentions) 
}
get.recipient <- function(tweet) {
  w <- strsplit(as.character(tweet)," ")[[1]][1]
  if (substr(w,0,1)=="@") {
    w <- gsub("@","",w)
    w <- gsub(":","",w)
  } else {
    w <- ""
  }
  return(w) 
}

mentions <- lapply(df$text,get.mentions)
save(mentions,file="data/rstats.mentions.rdata")

rec <- lapply(df$text,get.recipient)
length(rec)
nrow(df)
sum(rec!="")  # number messages
ix <- which(rec!="")
s <- as.character(df$user[ix])
r <- unlist(rec[ix])
datetime <- strptime(df$date[ix],format="%Y-%m-%dT%H:%M:%S")
A <- data.frame(datetime=datetime,s=s,r=r)
save(A,file="data/rstats.interaction.rdata")

# Create unique ids
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

# Rescale time to be in (0,1)
B <- A[,c(1,4,5)]
B[,1] <- as.numeric(B[,1])
B[,1] <- B[,1] - B[1,1]
B[,1] <- B[,1]/B[nrow(B),1]

z <- sample(1:K,N,replace=TRUE)
beta <- matrix(rnorm(K^2),K,K)

set.seed(4)
niter <- 100
px <- c(1,rep(0,6))
sbm.lpost(B,N,K,z,beta)

fit0 <- sbm.mcmc(B,N,2,niter=100)
fit1 <- sbm.mcmc(B,N,3,niter=100)
fit1 <- sbm.mcmc(B,N,10,niter=100)

fit <- fit0
nmap <- order(fit$z)
s <- match(B[,2],nmap)
r <- match(B[,3],nmap)
pdf("figs/twitter/sorted.pdf",width=4,height=4)
plot(s,r,pch=".",xlab="sender",ylab="recipient")
cpoints <- sapply(1:K,function(k) min(which(sort(fit$z)==k)))
abline(v=cpoints[-1])
abline(h=cpoints[-1])
dev.off()
fit$beta
)