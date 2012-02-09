library(inline)
library(Rcpp)

settings <- getPlugin("Rcpp")
settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)

src <- paste(readLines("brem.cpp"),collapse="\n")
fx <- cxxfunction(,"",includes=src, plugin="Rcpp",settings=settings)
brem <- Module("brem",getDynLib(fx))
