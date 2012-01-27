download <- function(urls) {
  ids <- as.vector(sapply(urls,function(u) strsplit(u,"/statuses/")[[1]][2]))
  mclapply(1:length(ids),function(i) {
    if (i %% 100 == 0) cat(i,"\n") 
    url <- paste("https://api.twitter.com/1/statuses/retweets/",ids[i],".json",sep="")
    req <- try(cred$OAuthRequest(url,"GET"))
    if (class(req) != "try-error") {
      return(fromJSON(req))
    } else {
      return(list("Error occurred"))
    }
  })
}
download.search <- function(hashtag,pages=1:10) {
  urls <- paste("http://search.twitter.com/search.json?q=%23",hashtag,"&include_entities=true&include_rts=true&rpp=100&page=",pages,sep="")
  res <- list()
  for (j in pages){
    toGet<-urls[j]
    req <- try(cred$OAuthRequest(toGet,"GET"))
    if (class(req) != "try-error") {
      res[[j]] <- fromJSON(req)
    } else {
      res[[j]] <- list("Error occurred")
    }
    if (j %% 10 == 0) cat("done with page ",j,"\n") 
  }
  return(res)
}
download.timeline <- function(usernames) {
  mclapply(1:length(usernames),function(i) {
    if (i %% 100 == 0) cat("done with user ",i,"\n") 
    urls <- paste("https://api.twitter.com/1/statuses/user_timeline.json?include_entities=true&include_rts=true&count=200&screen_name=",usernames[i],"&page=",1:16,sep="")
    res <- vector("list", 16)
    for (j in 1:16){
      toGet<-urls[j]
      req <- try(cred$OAuthRequest(toGet,"GET"))
      if (class(req) != "try-error") {
        res[[j]] <- fromJSON(req)
      } else {
        res[[j]] <- list("Error occurred")
      }
      if (j %% 10 == 0) cat("done with page ",j,"\n") 
    }
    return(res)
  })
}


backoff <- function(req) {
  # According to this: https://dev.twitter.com/docs/streaming-api (not correct I know...)
  Sys.sleep(.01)
}
chunk <- function(urls,keyword,size=10000) {
  if (length(urls) > size) {
    ms <- as.numeric(cut(1:length(urls),breaks=length(urls)/size))
  } else {
    ms <- rep(1,length(urls))
  }
  for (m in unique(ms)) {
    while (remaining() < size) {
      cat("queries left:",remaining())
      cat("waiting 10 minutes")
      Sys.sleep(60 * 10)
    }
    ix <- which(ms==m)
    d <- download(urls[ix])
    cat("Saving chunk of retweets.  Queries left:",remaining(),"\n")
    save(d,urls,keyword,file=paste("data/retweets/",keyword,"-",m,".rdata",sep=""))
  }
}
retweet.stats <- function(d) {
  ix <- sapply(d,function(x) length(x)>0 & is.list(x))
  cat("Number retweets:",sum(ix),"\n")
  return(ix)
}
get.retweet.df <- function(d) {
  # Get indices that have retweets
  ix <- which(sapply(d,function(x) length(x)>0 & is.list(x)))
  d <- d[ix]
  
}
getstatus <- function(id) {
  url <- paste("https://api.twitter.com/1/statuses/show/",id,".json",sep="")
  req <- cred$OAuthRequest(url,"GET")
  fromJSON(req)
}


# url <- "https://api.twitter.com/1/statuses/show/91226161736187907.json"
# reqs <- ds <- list()
# for (i in 1:100) {
#   url <- paste("https://api.twitter.com/1/statuses/retweets/",ids[i],".json",sep="")
#   req <- try(cred$OAuthRequest(url,"GET"))
#   reqs[[i]] <- req
#   if (class(req) != "try-error") {
#     ds[[i]] <- fromJSON(req)
#   } else {
#     ds[[i]] <- "Error occurred"
#   }
# }


remaining <- function() {
  url <- "http://api.twitter.com/1/account/rate_limit_status.json"
  fromJSON(cred$OAuthRequest(url,"GET"))$remaining_hits
}
