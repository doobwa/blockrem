
ix <- subset(llks,kinit==10 & iter > 10 & iter < 200 & coll==1 & llk < -13500 & deg==1)

load("results/eckmann-small/fits/kinit10.sm0.nb0.pshift1.deg1.trans1.collapse1.rdata")
load("data/eckmann-small.rdata")

# TEMP
ego <- TRUE
transform <- TRUE
  M <- nrow(train)
  P <- 13
  ego <- ego*1  # RemStat doesn't want boolean
  s <<- new(RemStat,
           train[,1],
           as.integer(train[,2])-1,
           as.integer(train[,3])-1,
           N,M,P,ego)
  s$precompute()
  if (transform) s$transform()

fit$llks[30:32]
fit$lps[30:32]
priors <- fit$priors
num.extra <- 50

# Before gibbs
params <- fit$samples[[30]]
lposterior(params,priors)$all
lposterior(params,priors)$y

# After just gibbs
params$z <- fit$samples[[31]]$z
lposterior(params,priors)$all
lposterior(params,priors)$y

# Double check
z <- fit$samples[[30]]$z
beta <- fit$samples[[30]]$beta
mu <- fit$samples[[30]]$mu
sigma <- fit$samples[[30]]$sigma

      prs <- priors
      prs$phi <- list(mu=mu,sigma=sigma)
      for (j in 1:num.extra) {
        beta <- add_cluster(beta)
        beta <- sample_cluster_from_prior(beta,dim(beta)[2],prs)
      }

h <- gibbs(beta,z,1:N,llk_node,lpz,N,priors)
table(fit$samples[[31]]$z,h$z)

a <- which(h$z==57)
K <- dim(beta)[2]
z[a] <- 2
sum(RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K))
sum(RemLogLikelihoodPc(beta, z-1, s$ptr(), K))
z[a] <- 57
sum(RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K))
sum(RemLogLikelihoodPc(beta, z-1, s$ptr(), K))

ix <- which(A[,2]==a | A[,3]==a)


params <- fit$samples[[30]]
params$z <- h$z
lposterior(params,priors)#[c("all","y")]


sum(RemLogLikelihoodPcSubset(params$beta, params$z - 1, s$ptr(), K,ix))

t(sapply(28:32,function(i) fit$samples[[i]]$mu))
fit$samples[[30]]$beta[1,,]
fit$samples[[31]]$beta[1,,]
table(fit$samples[[30]]$z,fit$samples[[31]]$z)
table(fit$samples[[30]]$z,fit$samples[[31]]$z)
z <- fit$samples[[30]]$z
table(z[train[,2]],z[train[,3]])
