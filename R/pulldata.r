# keyword <- "collar"
  
require(multicore)
require(rjson)
require(twitteR)
require(ROAuth)
require(RCurl)
source("pull.fns.r")
options(cores=1)
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

timelines <- download.timeline(chosen[1:2])
save(timelines,file="data/rstats.timelines.rdata")

# Grab messages from a list of tweets in json
get_messages <- function(tweets) {
  tmp <- sapply(tweets,function(tweet) {
    if (!is.null(tweet$to_user_id)) {
      tweet[c("id_str","created_at","to_user","to_user_id","from_user","from_user_id","text")]
    }
  })
  do.call(rbind,tmp)  
}


pages <- which(sapply(rstats,function(r) length(r$results)) == 0)
# check query limit
remaining()
