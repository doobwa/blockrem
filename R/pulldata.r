# keyword <- "collar"
  
require(multicore)
require(rjson)
require(twitteR)
require(ROAuth)
source("pull.fns.r")
options(cores=8)
load("heroic.oauth.rdata")

load(paste("data/rawdata/",keyword,"_orig.Rdata",sep=""))

# Samples 20000
set.seed(1)
dim(kwd.orig)
M <- min(20000,nrow(kwd.orig))
ix <- sample(1:nrow(kwd.orig),M)
kwd.samp <- kwd.orig[ix,]
save(kwd.samp,file=paste("data/sampled/",keyword,".rdata",sep=""))

load(file=paste("data/sampled/",keyword,".rdata",sep=""))
urls <- kwd.samp$tweetURL

# check query limit
remaining()

for (i in 0:3) {
  start <- 1 + 5000*i
  end <- 5000 * (i+1)
  # check rate limit
  rl<-remaining()
  if (rl>=5000){
    d <- download(urls[start:end])
    save(d,file=paste("data/retweets/",keyword,"-",i+1,".rdata",sep=""))
    inds<-unlist(lapply(d,is.character)) # inds where connection error occured
    while(sum(inds)>0){
      if (length(inds)>0) {
		rl<-remaining()
		if (rl > sum(inds)){
        	d2 <- download(urls[start:end][inds])
        	d[inds]<-d2
        	inds<-unlist(lapply(d,is.character))
		} else {
 			cat("sleeping \n")
   			Sys.sleep(15*60)
		}
      } 
      cat(".")
    }
    save(d,file=paste("data/retweets/",keyword,"-",i+1,".rdata",sep=""))
  } else {
    cat("sleeping \n")
    Sys.sleep(60*60)
  }
  }
}


# # Combine chunks
# chunk(urls,keyword,size=5000)
a <- list()
for (i in 1:length(grep(keyword,list.files("data/retweets")))) {
  load(paste("data/retweets/",keyword,"-",i,".rdata",sep=""))
  a <- c(a,d)
}

save(a,file=paste("data/retweets/",keyword,".rdata",sep=""))
