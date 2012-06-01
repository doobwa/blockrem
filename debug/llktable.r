library(brem)
#load("results/llk.rdata")  # using last draw
load("results/allpred.rdata") # averaged over draws

dataset <- "classroom-16"
models <- c("kinit2.kmax1.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha1000.xsigbeta1000",
            "kinit2.kmax2.sm0.nb0.pshift1.deg1.trans1.collapse1.xsigalpha1000.xsigbeta1000")

  load(paste("data/",dataset,".rdata",sep=""))
  M <- nrow(A)
  P <- 13
  ego <- 1
  train.ix <- 1:nrow(train)
  test.ix  <- (1:nrow(A))[-train.ix]
  s <- new(RemStat,A[,1],as.integer(A[,2])-1,as.integer(A[,3])-1,N,M,P,ego)
  s$precompute()
  s$transform()

folder <- paste("results/",dataset,"/fits/",sep="")
llks <- lapply(models,function(model) {
  load(paste(folder,model,".rdata",sep=""))
  K <- dim(fit$params$beta)[2]
  list(train = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, train.ix - 1)[train.ix],
       test  = RemLogLikelihoodPcSubset(fit$params$beta, fit$params$z - 1, s$ptr(),  K, test.ix  - 1)[test.ix])
})
names(llks) <- models

# Averaged over samples
a <- allpred[[dataset]][[models[1]]]
b <- allpred[[dataset]][[models[2]]]
mean(a$llk$test)
mean(b$llk$test)

# Using last sample
mean(llks[[models[1]]]$test)
mean(llks[[models[2]]]$test)

