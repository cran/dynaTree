#*******************************************************************************
#
# Dynamic trees for learning, design, variable selection, and sensitivity
# Copyright (C) 2011, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## dynaTree:
##
## Initialization and PL for dynamic tree models

dynaTree <-
  function(X, y, N=1000, model=c("constant", "linear", "class", "prior"),
           nu0s20=c(0,0), ab=c(0.95, 2), minp=NULL, sb=NULL,
           nstart=minp, icept=c("implicit", "augmented", "none"),
           rprop=c("luvar", "luall", "reject"), verb=round(length(y)/10))
  {
    ## extract vitals of X, and check dims
    X <- as.matrix(X)
    m <- ncol(X)
    T <- nrow(X)
    if(T != length(y)) stop("dim of X and Y mismatch")

    ## check model and encode as integer
    model <- match.arg(model)
    if(model == "constant") imodel <- 1
    else if(model == "linear") imodel <- 2
    else if(model == "class") { 
      imodel <- 3
      y <- round(y)-1 ## check for sain class labels
      if(any(y < 0)) stop("class labels must start at 1")
      if(length(setdiff(y, 0:max(y))) != 0)
        warning("y without one label in each class")
    } else imodel <- 4 ## or sample from prior

    ## check splitmin and basemax
    if(!is.null(sb)) {
      if(model != "linear") stop("sb only makes sense for linear model\n")
      if(length(sb) != 2) stop("must have length(sb) = 2")
      if(sb[1] <= 0 || sb[1] > m) stop("must have 1 <= sb[1] <= ncol(X)")
      if(sb[2] <= 0 || sb[2] > m) stop("must have 1 <= sb[2] <= ncol(X)")
    } else sb <- c(1, m)

    ## default minimum number of data points in each parition
    if(is.null(minp)) {
      if(model == "constant") minp <- 4
      else if(model == "linear") minp <- 2*sb[2] + 4
      else minp <- 1 ## for classify or prior
    } else if(length(minp) != 1 || minp <= 0) stop("minp must be a positive integer")

    ## check nstart
    if(is.null(nstart)) nstart <- 2*minp
    else if(length(nstart) != 1 || nstart < minp) stop("nstart must be >= minp")
    if(length(y) <= nstart) stop("must have more than nstart X-y pairs")
    
    ## check intercept
    icept <- match.arg(icept)
    if(icept != "implicit" && model != "linear")
      stop("icept != \"implicit\" only valid for linear model")
    if(icept == "augmented") {
      X <- cbind(rep(1, T), X)
      m <- m + 1; sb <- sb + 1
    }
    if(icept == "implicit") icepti <- 1
    else icepti <- 0

    ## check for missing data
    tXNA <- Xna <- NULL; NAX <- is.na(X)
    if(any(NAX)) {
      o  <- order(apply(NAX, 1, sum))
      X <- X[o,]; y <- y[o]; NAX <- NAX[o,]   ## NOT SURE THIS IS NECESSARY
      Xna <- apply(NAX, 1, any) 
      tXNA <- t(NAX[Xna,])
      X[is.na(X)] <- -12345 ## sanity checking code
    }

    ## double check that minp is largest that longest initial run
    if(model != "class" && model != "prior") {
      if(length(unique(y[1:minp])) == 1)
        stop("initial minp run in y must have at 2+ unique values")
    }

    ## check variance prior parameters 
    if(length(nu0s20) != 2 || nu0s20[1] < 0 || nu0s20[2] < 0)
      stop("must have nu0s20[1] >= 0 and nu0s20[2] >= 0")
    if(nu0s20[1] == 0 && nu0s20[2] != 0)
      stop("must have nu0s20[2] == 0 when nu0s20[1] == 0")
    
    ## check tree prior parameters 
    if(length(ab) != 2 || ab[1] < 0 || ab[1] >= 1 || ab[2] <= 0)
      stop("must have 0 <= ab[1] < 1 and ab[2] > 0 ")

    ## check rprop
    rprop <- match.arg(rprop)
    if(rprop == "luall") irprop <- 1
    else if(rprop == "luvar") irprop <- 2
    else irprop <- 3

    ## collect all parameters
    params <- c(nu0s20, ab, minp, sb, icepti, irprop)

    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## call the C routine to build up the PL object
    obj <- .C("dynaTree_R",
              m = as.integer(m),
              T = as.integer(T),
              N = as.integer(N),
              X = as.double(t(X)),
              bna = as.integer(any(NAX)),
              Xna = as.integer(Xna),
              XNA = as.integer(tXNA),
              y = as.double(y),
              model = as.integer(imodel),
              params = as.double(params),
              nstart = as.integer(nstart),
              verb = as.integer(verb),
              lpred = double(T),
              num = integer(1),
              PACKAGE = "dynaTree")

    ## end timing
    obj$time <- proc.time()[3] - p1
    
    ## put non-transposed X back and model
    if(!is.null(tXNA)) {
      X[NAX] <- NA
      obj$XNA <- t(tXNA)
    } else obj$XNA <- NULL
    obj$X <- X
    obj$Xna <- Xna
    obj$T <- NULL
    obj$bna <- NULL
    obj$model <- model
    if(obj$model == "class") {
      obj$y <- y + 1
      obj$classes <- sort(unique(obj$y))
    }
    obj$icept <- icept
    
    ## assign class and return
    class(obj) <- "dynaTree"
    invisible(obj)
  }


## update:
##
## update the particle cloud to include new (x,y) pairs

update.dynaTree <- function(object, X, y, verb=round(length(y)/10), ...)
  {
    ## for timing purposes
    p1 <- proc.time()[3]

    ## sanity check
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## extract vitals of X
    m <- object$m
    X <- as.matrix(X)
    T <- nrow(X)
    if(object$icept == "augmented") X <- cbind(rep(1,T), X)
    if(ncol(X) != m) stop("column mismatch for X")

    ## check for missing data
    tXNA <- Xna <- NULL; NAX <- is.na(X)
    if(any(NAX)) {
      o  <- order(apply(NAX, 1, sum))
      X <- X[o,]; y <- y[o]; NAX <- NAX[o,]  ## NOT SURE THIS IS NECESSARY
      Xna <- apply(NAX, 1, any) 
      tXNA <- t(NAX[Xna,])
      X[is.na(X)] <- -12345 ## sanity checking code
    }

    ## ensure new class labels are valid
    if(object$model == "class") {
      y <- round(y)-1 ## check for sain class labels
      if(any(y < 0)) stop("class labels must start at 1")
      if(any(y >= max(object$classes)))
        stop("class labels must be <= max(object$classes)")
    }

    ## echo something to the screen
    if(verb > 0 && T > verb) cat("updating with", T, "new (x,y) pairs\n");

    ## call the C routine to build up the PL object
    object2 <- .C("update_R",
                  cloud = as.integer(object$num),
                  m = as.integer(m),
                  T = as.integer(T),
                  X = as.double(t(X)),
                  bna = as.integer(any(NAX)),
                  Xna = as.integer(Xna),
                  tXNA = as.integer(tXNA),
                  y = as.double(y),
                  verb = as.integer(verb),
                  lpred = double(T),
                  PACKAGE = "dynaTree")
    
    ## remove cloud
    object2$cloud <- NULL
    
    ## put non-transposed X back, and combine
    if(object$model == "class") y <- y + 1
    if(!is.null(tXNA)) {
      X[NAX] <- NA
      obj$XNA <- rbind(object$XNA, t(tXNA))
    }
    object$X <- rbind(object$X, X)
    object$Xna <- c(object$Xna, Xna)
    object$y <- c(object$y, y)
    object$lpred <- c(object$lpred, object2$lpred)

    ## update time
    object$time <- object$time + proc.time()[3] - p1

    ## assign class and return
    class(object) <- "dynaTree"
    invisible(object)
  }


## retire.dynaTree:
##
## move the specified indices into the prior, retireing them,
## i.e., removing them from the marginal likelihood calculation

setGeneric("retire",
            function(object, indices, lambda=1, verb=0)
            standardGeneric("retire")
            )

retire.dynaTree <- function(object, indices, lambda=1, verb=0)
  {
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## must have explicit intercept or none
    if(object$model == "linear" && object$icept == "implicit")
      stop("must use explicit intercept (i.e., augmented or none)")
       
    ## check to make sure we're not removing non-existant indices
    n <- nrow(object$X)
    indices <- unique(indices)
    if(length(union(1:n, indices)) > n)
      stop("indices must lie in 1:nrow(object$X)")

    ## check lambda
    if(length(lambda) != 1 || lambda <= 0 || lambda > 1)
      stop("lambda must be a postive scalar proportion")

    ## new data size
    removed <- length(indices)
    nnew <- n - removed
    m <- ncol(object$X)
    
    ## check that we're not removing too many
    if(nnew == 0) stop("cannot remove all rows of object$X")

    out <- .C("retire_R",
              cloud = as.integer(object$num),
              indicies = as.integer(indices-1),
              ilen = as.integer(removed),
              lambda = as.double(lambda),
              verb = as.integer(verb),
              X = double(nnew * m),
              y = double(nnew),
              PACKAGE = "dynaTree")

    ## copy new X and y into object
    object$X <- matrix(out$X, ncol=m, byrow=TRUE)
    if(object$model == "class") out$y <- out$y + 1
    object$y <- out$y
    if(is.null(object$removed)) object$removed <- removed
    else object$removed <- object$removed + removed

    ## update time
    object$time <- object$time + proc.time()[3] - p1
    
    return(object)
  }

setMethod("retire", "dynaTree", retire.dynaTree)


## delete.cloud:
##
## deletes the C-side in a particular 

deletecloud <- function(obj)
  {
    if(is.null(obj$num)) stop("no cloud number in object")
    .C("delete_cloud_R",
       num = as.integer(obj$num),
       PACKAGE = "dynaTree")
    invisible(NULL)
  }


## deleteclouds:
##
## deletes all dynatree clouds on the C side

deleteclouds <- function()
  {
    .C("delete_clouds_R", PACKAGE="dynaTree")
    invisible(NULL)
  }


## copy.dynaTree:
##
## copyies the entire object, also duplicating the clouds
## on the C side

setGeneric("copy",
            function(obj)
            standardGeneric("copy")
            )

copy.dynaTree <- function(obj)
  {
    if(is.null(obj$num)) stop("no cloud number in object")
    r <- .C("copy_cloud_R",
            num = as.integer(obj$num),
            newnum = integer(1),
            PACKAGE="dynaTree")
    obj$num <- r$newnum
    return(obj)
  }

setMethod("copy", "dynaTree", copy.dynaTree)


## rejuvenate.dynaTree:
##
## re-initializes a particle set and combines it with the old
## one

setGeneric("rejuvenate",
            function(object, odr=order(runif(length(object$y))),
                                verb=round(length(object$y)/10))
            standardGeneric("rejuvenate")
            )

rejuvenate.dynaTree <- function(object, odr=order(runif(length(object$y))),
                                verb=round(length(object$y)/10))
  {
    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## check the cloud
    if(is.null(object$num)) stop("no cloud number in object")

    ## sanity check o is a reordering of 1:length(object$y)
    n <- length(odr)
    if(!is.null(odr)) {
      odr <- round(odr)
      if(n != length(object$y)) 
        stop("odr should be a length(object$y) vector")
      if(min(odr) <= 0 || max(odr) > n || length(unique(odr)) != n)
        stop("odr should be a reordering of 1:length(object$y)")
    }

    ## perhaps print something
    ## if(verb > 0) cat("rejuvenating particles\n")
    
    ## call C-side rejuvenate
    r <- .C("rejuvenate_R",
            num = as.integer(object$num),
            odr = as.integer(odr-1),
            n = as.integer(n),
            verb = as.integer(verb),
            lpred = double(length(object$y)),
            PACKAGE="dynaTree")

    ## update time
    object$time <- object$time + proc.time()[3] - p1
    
    return(object)
  }

setMethod("rejuvenate", "dynaTree", rejuvenate.dynaTree)


## dynaTrees:
##
## calls dynaTree and then predict R times in order to asses
## the Monte Carlo error of the PL procedure and aggregate the
## predictive distributions of many re-roderings of the data

"dynaTrees"<-
  function(X, y, N=1000, R=10, sub=length(y),
           model=c("constant", "linear", "class", "prior"),
           nu0s20=c(0,0), ab=c(0.95, 2), minp=NULL, sb=NULL,
           nstart=minp,icept=c("implicit", "augmented", "none"),
           rprop=c("luvar", "luall", "reject"), XX=NULL, yy=NULL,
           varstats=FALSE, lhs=NULL, plotit=FALSE, proj=1, rorder=TRUE,
           verb=round(sub/10), pverb=round(N/10),  ...)
  {
    ## use dynaTree and predict by themselves when R = 1
    if(R <= 1) stop("R should be >= 2")

    ## coerse X
    X <- as.matrix(X)
    n <- nrow(X)

    ## check model
    model <- match.arg(model)
    if(model == "prior" && !is.null(XX))
      stop("cannot predict at XX for prior model")
    
    ## check rorder
    if(length(rorder) > 1) {
      if(nrow(rorder) != sub && ncol(rorder) != R)
        stop("bad rorder argument")
      else o <- rorder
    } else o <- apply(matrix(runif(sub*(R-1)), ncol=R-1), 2, order)
    o <- cbind((1:sub), o)

    ## check varstats
    if(length(varstats) != 1 || !is.logical(varstats))
      stop("varstats should be a scalar logical")

    ## build the first model
    if(verb > 0) cat("\nround 1:\n")
    obj <- dynaTree(X[1:sub,], y[1:sub], N, model, nu0s20, ab, minp, sb, nstart, icept, rprop, verb)

    ## predict or perform sensitivity analysis
    if(!is.null(XX)) {
      if(is.character(XX) && XX == "sens") {
        if(!is.null(yy)) warning("yy ignored in sensitivity analysis")
        obj <- sens.dynaTree(obj, lhs=lhs, verb=pverb, ...)
      } else obj <- predict(obj, XX, yy, verb=pverb, ...)
    }

    ## maybe accumulate variable use proportions
    if(varstats) {
      obj$vpu <- varpropuse(obj)
      obj$vpt <- varproptotal(obj)
      obj <- relevance(obj, verb=pverb, ...)
    }
    
    ## delete cloud
    deletecloud(obj); obj$num <- NULL

    ## possibly plot in 1d case
    if(plotit) {
      if(is.null(XX))
        warning("cannot plot without XX predictive grid", immediate.=TRUE)
      if(is.character(XX) && XX == "sens")
        warning("sens plots not implemented yet", immediate.=TRUE)
      else plot(obj, proj=proj)
    }

    ## now do the same ting R-1 more times and combine outputs
    for(r in 2:R) {

      ## build the Rth model on a the random re-ordering
      if(verb > 0) cat("\nround ",  r, ":\n", sep="")
      obj2 <- dynaTree(X[o[,r],], y[o[,r]], N, model, nu0s20, ab, minp, sb,
                       nstart, icept, rprop, verb)

      ## predict or perform sensitivity analysis
      if(!is.null(XX)) {
        if(is.character(XX) && XX == "sens")
          obj2 <- sens.dynaTree(obj2, lhs=lhs, verb=pverb, ...)
        else obj2 <- predict(obj2, XX, yy, verb=pverb, ...)
      }

      ## maybe accumulate variable use proportions
      if(varstats) {
        obj2$vpu <- varpropuse(obj2)
        obj2$vpt <- varproptotal(obj2)
        obj2 <- relevance(obj2, verb=pverb, ...)
      }

      ## delete cloud
      deletecloud(obj2); obj2$num <- NULL

      ## possibly add to the plot in the 1d/non-sens case
      if(plotit) {
        if(!is.null(XX) && !(is.character(XX) && XX == "sens"))
          plot(obj2, add=TRUE, proj=proj)
      }

      ## combine the PL bits of the object
      obj$N <- obj$N + obj2$N
      obj$lpred <- cbind(obj$lpred, obj2$lpred)

      ## combine times
      obj$time <- c(obj$time, obj2$time)

      ## combine the predictive bits
      if(!is.null(XX)) {
        if(is.character(XX) && XX == "sens") {
          ## sensitivity collecting
          if(model != "class") { ## regression
            obj$MEmean <- ((r-1)*obj$MEmean + obj2$MEmean)/r
            obj$MEq1 <- ((r-1)*obj$MEq1 + obj2$MEq1)/r
            obj$MEq2 <- ((r-1)*obj$MEq2 + obj2$MEq2)/r
            obj$S <- rbind(obj$S, obj2$S)
            obj$T <- rbind(obj$T, obj2$T)
          } else { ## classification
            for(i in obj2$sens.class) {
              obj$MEmean[[i]] <- ((r-1)*obj$MEmean[[i]] + obj2$MEmean[[i]])/r
              obj$MEq1[[i]] <- ((r-1)*obj$MEq1[[i]] + obj2$MEq1[[i]])/r
              obj$MEq2[[i]] <- ((r-1)*obj$MEq2[[i]] + obj2$MEq2[[i]])/r
              obj$S[[i]] <- rbind(obj$S[[i]], obj2$S[[i]])
              obj$T[[i]] <- rbind(obj$T[[i]], obj2$T[[i]])
            }
          }
        } else {
          if(obj$model != "class") {
            ## regression collecting
            obj$mean <- cbind(obj$mean, obj2$mean)
            obj$vmean <- cbind(obj$vmean, obj2$vmean)
            obj$var <- cbind(obj$var, obj2$var)
            obj$df <- cbind(obj$df, obj2$df)
            obj$q1 <- cbind(obj$q1, obj2$q1)
            obj$q2 <- cbind(obj$q2, obj2$q2)
            obj$alc <- cbind(obj$alc, obj2$alc)
          } else { ## classification averaging
            obj$p <- ((r-1)*obj$p + obj2$p)/r
            obj$entropy <- cbind(obj$entropy, obj2$entropy)
          }
          if(!is.null(yy)) obj$yypred <- cbind(obj$yypred, obj2$yypred)
        }
      }

      ## combine the variable use bits
      if(varstats) {
        obj$vpu <- rbind(obj$vpu, obj2$vpu)
        obj$vpt <- rbind(obj$vpt, obj2$vpt)
        obj$relevance <- rbind(obj$relevance, obj2$relevance)
      }
    }
    
    ## print for next round
    if(verb > 0) cat("\n")

    ## assign R and class and return
    obj$R <- R
    class(obj) <- "dynaTree"
    invisible(obj)
  }


## intervals:
##
## returns the upper and lower bounds for the column
## of all tree partitions used by any X[index,] in the object

intervals.dynaTree <- function(object, index, var)
  {
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")

    ## check index
    n <- nrow(object$X)
    if(length(index) != 1 || index <= 0 || index > n)
      stop("index must be scalar in  in 1:nrow(object$X)")

    ## check var
    m <- ncol(object$X)
    if(length(var) != 1 || var < 1)
      stop("var must be scalar >= 1 and <= ncol(X)")
    if(object$icept == "augmented") var <- var + 1
    if(var > m) stop("var must be scalar >= 1 and <= ncol(X)")

    out <- .C("intervals_R",
              cloud = as.integer(object$num),
              index = as.integer(index),
              var = as.integer(var),
              a = double(object$N),
              b = double(object$N),
              PACKAGE = "dynaTree")

    return(data.frame(a=out$a, b=out$b))
  }


## treestats:
##
## returns the upper and lower bounds for the column
## of all tree partitions used by any X[index,] in the object

setGeneric("treestats",
            function(object)
            standardGeneric("treestats")
            )

treestats.dynaTree <- function(object)
  {
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")

    out <- .C("treestats_R",
              cloud = as.integer(object$num),
              avgheight = double(1),
              avgleaves = double(1),
              avgsize = double(1),
              avgretire = double(1),
              PACKAGE = "dynaTree")

    return(data.frame(avgheight=out$avgheight, avgleaves=out$avgleaves,
                      avgsize=out$avgsize, avgretire=out$avgretire))
  }

setMethod("treestats", "dynaTree", treestats.dynaTree)
