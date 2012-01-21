
beta <- c(1,0)
effects <- names(beta) <- c("PSAB-BA","PSAB-BY")
A <- simulate.rem(1000,10,effects,beta)
s <- build.brem(A,N)

N <- 10
K <- 2
z <- c(rep(1,N/2),rep(2,N/2))
P <- 9
beta <- array(0,c(K,K,P))
beta[1,1,] <- c(1,0,0,4,0,0,0,0,0)
beta[2,1,] <- c(1,0,0,0,4,0,0,0,0)
beta[1,2,] <- c(1,0,0,0,4,0,0,0,0)  # AB-
beta[2,2,] <- c(1,0,0,0,0,0,0,0,2)  # AB-AY


M <- 1000
A <- simulate.brem(M,N,z,beta)

df <- block.ps(A,z)
qplot(pshift,value,data=df,geom="bar") + facet_grid(i~j) + theme_bw()

precomp <- rem.dyad.setup(A,10,effects=effects)
x <- rem.dyad(A,10,effects,ordinal=FALSE)

value <- list(ms=precomp$modstats,pv=c(1,0))
rem.llk(value)