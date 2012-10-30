#' At each time t:
#' for each (i,j) 
#'   compute lambda_ij(t|.) = sum_kl exp{ phi_ik + phi_jl + beta_kl s_ij(t|.)}
#' draw (i,j) ~ lambda_ij/sum_ij lambda_ij
#' draw next t ~ Exp(sum_ij lambda_ij)

CLogLambda <- function(i,j,beta,phi,s) {
  lam <- 0
  for (k in 1:K) {
    for (l in 1:K) {
      lam <- lam + 
    }
  }
}

simulate.cbrem.cond <- function(beta,phi,train,M,N,ego=FALSE) {
  K <- dim(beta)[2]
  P <- dim(beta)[1]
  edgelist <- train[nrow(train),,drop=FALSE]
  time <- edgelist[1]
  s <- InitializeStatisticsArray(N,P);
  for (m in 1:nrow(train)) {
    i <- train[m,2]
    j <- train[m,3]
    s <- UpdateStatisticsArray(s,m-1,i-1,j-1,N,P)
  }
  lambda <- matrix(0,N,N)
  for (i in 1:N) {
    for (j in 1:N) {
      lambda[i,j] <- LogLambda(i-1,j-1,z[i]-1,z[j]-1,s,beta,N,K,P)
    }
  }
  for (m in 1:M) {
    for (r in 1:N) {
      lambda[j,r] <- LogLambda(j-1,r-1,z[j]-1,z[r]-1,s,beta,N,K,P)
      lambda[i,r] <- LogLambda(i-1,r-1,z[i]-1,z[r]-1,s,beta,N,K,P)
      if (!ego) {
        lambda[r,j] <- LogLambda(r-1,j-1,z[r]-1,z[j]-1,s,beta,N,K,P)
        lambda[r,i] <- LogLambda(r-1,i-1,z[r]-1,z[i]-1,s,beta,N,K,P)
      }
    }
    diag(lambda) <- -Inf
    cells <- cbind(as.vector(row(lambda)), as.vector(col(lambda)), exp(as.vector(lambda)))
    if (any(is.nan(cells[,3]))) {
      stop("Intensity explosion: infinite lambda_ij present")
    }
    drawcell <- sample(1:NROW(cells),1,prob=cells[,3])
    i <- cells[drawcell,1]
    j <- cells[drawcell,2]
    time <- time + rexp(1,sum(cells[,3]))
    edgelist <- rbind(edgelist,c(time,i,j))

    s <- UpdateStatisticsArray(s,m-1,i-1,j-1,N,P)
  }
  dimnames(edgelist) <- NULL
  edgelist <- edgelist[-1,]
  return(list(edgelist=edgelist,s=s,N=N))
}


lpz <- function(count,priors) {
  lp <- switch(priors$type,
               "crp" = log(count + priors$alpha),
               "nb"  = dnbinom(count,priors$nb$shape,mu=priors$nb$mu,log=TRUE),
               "bcrp"= 3 * log(count + priors$alpha))
}

##' RemStat object s required in environment
##' @title 
##' @param params 
##' @param priors 
##' @param collapse.sigma 
##' @return 
##' @author chris
lposterior <- function(params,priors,collapse.sigma=TRUE) {
  K <- dim(params$beta)[2]
  P <- dim(params$beta)[1]
  if (s$get_P() != P) stop("Mismatch in parameter vector lengths")

  # Likelihood
  y <- sum(RemLogLikelihoodPc(params$beta,params$z-1,s$ptr(),K))

  # Prior on beta
  if (collapse.sigma) {
    beta <- lprior.nosigma.hier.gaussian(params$beta,params$mu,priors,grad=FALSE)
  } else {
    beta <- lprior.hier.gaussian(params$beta,params$mu,params$sigma,priors,grad=FALSE)
  }

  # Prior on z
  tb <- table(factor(params$z,1:K))
  z <- sapply(tb,function(count) count * lpz(count,priors))

  # Prior on sigma
  sigma <- dgamma(params$sigma,priors$sigma$alpha,priors$sigma$beta,log=TRUE)

  r <- list(y=y,beta=beta,z=z,sigma=sigma)
  r$all <- sum(r$sigma[priors$px]) + y + beta + sum(z)
  return(r)
}

llk_node <- function(a,beta,z,priors) {
  K <- dim(beta)[2]
  sum(RemLogLikelihoodActorPc(a-1,beta,z-1,s$ptr(),K))
#  sum(RemLogLikelihoodPc(beta,z-1,s$ptr(),K))
}

##' Requires s to be in environment
##' @title Sample block level parameters beta
##' @param beta 
##' @param z 
##' @param mu 
##' @param sigma 
##' @param priors 
##' @param kx 
##' @return 
##' @author chris
sample_beta <- function(beta,z,mu,sigma,priors,kx=NULL,collapse.sigma=TRUE) {
  K <- dim(beta)[2]
  if (is.null(kx)) kx <- 1:K
  for (k1 in kx) {
    for (k2 in kx) {
      cat(".")
      olp <- NULL
      for (p in priors$px) {
        brem.lpost.pk1k2 <- function(x) {
          k1nodes <- which(z==k1)
          k2nodes <- which(z==k2)
          beta[p,k1,k2] <- x
          pr.y <- sum(RemLogLikelihoodBlockPc(k1-1,k2-1,k1nodes-1,k2nodes-1,beta,z-1,s$ptr(),K))
#          pr.beta <- sum(dnorm(beta[,k1,k2],mu,sigma,log=TRUE))
          if (collapse.sigma) {
            pr.beta <- lprior.nosigma.hier.gaussian(beta[,k1,k2],mu,priors,grad=FALSE)
          } else {
            pr.beta <- lprior.hier.gaussian(beta[,k1,k2],mu,sigma,priors,grad=FALSE)
          }
          return(pr.y + pr.beta)
        }
        res <- slice(beta[p,k1,k2],brem.lpost.pk1k2,m=20,olp=olp)
        olp <- attr(res,"log.density")
        beta[p,k1,k2] <- res
      }
    }
  }
  return(list(beta=beta,z=z))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param train edgelist
##' @param N maximum id of actors
##' @param priors alpha, sigma.proposal, phi=list(mu=0,sigma=1), mu=list(mu=0,sigma=1), sigma=list(alpha=2,beta=2))
##' @param kinit 
##' @param kmax 
##' @param effects character vector that is a subset of "intercept","abba","abby","abxa","abxb","abay","abab","sen_outdeg","rec_outdeg","sen_indeg","rec_indeg","dyad_count","changepoint_count"
##' @param ego 
##' @param transform 
##' @param do.sm perform split merge moves
##' @param num.extra number of extra clusters to propose from the prior each iteration
##' @param niter number of iterations
##' @param collapse.sigma 
##' @param verbose 
##' @param outfile save progress to this file
##' @param K initial number of clusters
##' @return 
##' @author chris
brem_continuous <- function(train,N,priors,k,effects=c("intercept","abba","abby","abay"),ego=TRUE,transform=FALSE,do.sm=FALSE,num.extra=2,niter=20,collapse.sigma=TRUE,verbose=TRUE,outfile=NULL) {
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

  enam <- c("intercept","abba","abby","abxa","abxb","abay","abab",
            "sen_outdeg","rec_outdeg","sen_indeg","rec_indeg",
            "dyad_count","changepoint_count")
  priors$px <- match(effects,enam)

  # Initialize parameters
  K <- kinit
  beta <- array(0,c(P,K,K))
  beta[1,,] <- rnorm(K^2)
  
  z     <- sample(1:K,N,rep=T)
  mu    <- rep(0,P)
  sigma <- rinvgamma(P,priors$sigma$alpha,priors$sigma$beta)
  beta  <- sample_beta(beta,phi,mu,sigma,priors,collapse.sigma=collapse.sigma)$beta

  samples <- vector('list',length=niter)
  lps <- llks <- acc <- rep(0,niter)
  for (iter in 1:niter) {

    ## Gibbs sample assignments
    h <- gibbs(beta,phi,1:N,llk_node,lpz,N,priors)
    z <- h$z

    ## Sample phi
    beta <- sample_beta(beta,phi,mu,sigma,priors)$beta
    
    ## Fit hierarchical portion
    if (K != 1) {
      theta <- t(array(beta,dim=c(P,K^2)))
      pr <- list(theta=theta,mu=mu,sigma=sigma)
      mu <- gibbs.mu.hier.gaussian(pr,priors)$mu
      sigma <- gibbs.sigma.hier.gaussian(pr,priors)$sigma
    }
    
    ## Save progress
    # not included in likelihood, so set to 0 for display purposes
    mu[-priors$px] <- sigma[-priors$px] <- 0
    params <- list(beta=beta,z=z,mu=mu,sigma=sigma,ego=ego)
    lp <- lposterior(params,priors)
    lps[iter] <- lp$all
    llks[iter] <- lp$y
    cat("\n",iter,": llk",llks[iter]," lp",lps[iter],"\n")
    if (verbose) cat(z,"\n")
    samples[[iter]] <- params
    fit <- list(params=params,samples=samples,ego=ego,transform=transform,priors=priors,lps=lps,llks=llks,effects=effects,lp=lp,iter=iter,niter=niter)
    class(fit) <- "brem"
    if (!is.null(outfile)) save(fit,file=outfile)
  }

  return(fit)
}
