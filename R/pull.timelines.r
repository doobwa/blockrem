require(multicore)
require(rjson)
require(twitteR)
require(ROAuth)
require(RCurl)
source("R/pull.fns.r")
options(cores=1)  # seems necessary. otherwise twitter blocks.
load("data/heroic.oauth.rdata")
load("data/rstats.urls.rdata")
req <- try(cred$OAuthRequest(urls[1],"GET"))
tmp <- fromJSON(req)
if (length(tmp) != 200) print("Timeline retrieval broken")

get.timelines <- function(urls) {
  lapply(1:length(urls),function(i) {
    req <- try(cred$OAuthRequest(urls[i],"GET"))
    if (class(req) != "try-error") {
      req <- fromJSON(req)
    } else {
      req <- list("Error occurred")
    }
    if (i %% 10 == 0) cat("done with url ",i,"\n")
    if (length(req) == 200)
      save(req,file=paste("data/timelines/",i,".rdata",sep=""))
    return(req)
  })
}
reqtmp <- get.timelines(urls[1:50])
if (length(reqtmp[[1]]) != 200) print("Timeline retrieval broken")
if (length(reqtmp[[33]]) != 200) print("Timeline retrieval broken")
sapply(reqtmp,function(r) length(r))

# Get list of urls still needed
fs <- list.files("data/timelines")
fs <- gsub(".rdata","",fs)
ix <- (1:length(urls))[-as.numeric(fs)]
ix[1:50]
reqs <- get.timelines(urls[58:100])
save(reqs,file="data/rstats.timelines.rdata")




# Grab messages from a list of tweets in json
get_messages <- function(tweets) {
  tmp <- sapply(tweets,function(tweet) {
    if (!is.null(tweet$to_user_id)) {
      tweet[c("id_str","created_at","to_user","to_user_id","from_user","from_user_id","text")]
    }
  })
  do.call(rbind,tmp)  
}

