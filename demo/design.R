## illustrates sequential design and optimization
## by active learning heuristics; the default is
## ALM, but you can also try ALC and EI

## load the libraries
library(dynaTree)
## dyn.load("../dynaTree/src/dynaTree.so")
library(tgp)

## describing the (x,y) data
f1d <- function(x, sd=0.1){
  return( sin(x) - dcauchy(x,1.6,0.15) + rnorm(1,0,sd))
} 
rect <- c(0,7)

## initial (x,y) data
start <- 10
X <- dopt.gp(start, Xcand=lhs(10*start, rect))$XX
y <- f1d(X)

## size of predictive grid and type of AS
ngrid <- 100
method <- "alm" ## also try "alc", "ei", or "ieci"
prec <- 0.1 ## for ei

## PL fit to initial data
obj <- dynaTree(X=X, y=y, model="linear")

## determining the number of adaptive sampling rounds
end <- 100

##
## Do the sequential design
##

track <- NULL

for(t in start:end){

  ## random predictive grid
  XX <- lhs(ngrid, rect)

  ## predict arguments
  ieci <- alc <- ei <- FALSE
  if(method == "alc") { alc <- TRUE
  } else if(method == "ei") { ei <- TRUE 
  } else if(method == "ieci") ieci <- TRUE
  
  ## predict at the XX locations
  obj <- predict(obj, XX, quants=FALSE, ei=ei)
  if(alc) { obj <- alc(obj, XX, rect=rect)
  } else if(ieci) obj <- ieci(obj, XX)

  ## extract via ALM, ALC, EI-prec
  al <- alcalc(obj, method, prec)
  m <- which.max(al)
  track <- c(track, al[m])
  xstar <- drop(obj$XX[m,])
  ystar <- f1d(xstar)

  ## plot the progress
  plotprogress(obj, xstar, ystar, method, track, f1d)
  
  ## update the fit for the next round
  obj <- update(obj, xstar, ystar)
}

## free the particle cloud on the C-side
deletecloud(obj); obj$num <- NULL

##dev.off()
