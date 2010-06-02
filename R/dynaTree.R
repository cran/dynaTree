## dynaTree:
##
## Initialization and PL for dynamic tree models

dynaTree <- function(X, y, N=1000, model=c("constant", "linear", "class"),
                     alpha=0.95, beta=2, minp=NULL, verb=10)
  {
    ## extract vitals of X
    X <- as.matrix(X)
    m <- ncol(X)
    T <- nrow(X)

    ## check model and encode as integer
    model <- match.arg(model)
    if(model == "constant") imodel <- 1
    else if(model == "linear") imodel <- 2
    else { ## or classify
      imodel <- 3
      y <- round(y)-1 ## check for sain class labels
      if(any(y < 0)) stop("class labels must start at 1")
      if(length(setdiff(y, 0:max(y))) != 0)
        warning("y without one label in each class")
    }

    ## default minimum number of data points in each parition
    if(is.null(minp)) {
      if(model == "constant") minp <- 3
      else if(model == "linear") minp <- 2*ncol(X) + 3
      else minp <- 1 ## for classify
    }

    ## check tree prior parameters and length of data
    if( alpha < 0 || alpha >= 1 || beta <=0 || T != length(y) )
      stop("bad Params")
    params <- c(alpha, beta, minp)

    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## call the C routine to build up the PL object
    obj <- .C("dynaTree_R",
              m = as.integer(m),
              T = as.integer(T),
              N = as.integer(N),
              X = as.double(t(X)),
              y = as.double(y),
              model = as.integer(imodel),
              params = as.double(params),
              verb = as.integer(verb),
              lpred = double(T),
              num = integer(1),
              PACKAGE = "dynaTree")

    ## end timing
    obj$time <- proc.time()[3] - p1
    
    ## put non-transposed X back and model
    obj$X <- X
    obj$model <- model
    if(obj$model == "class") obj$y <- y + 1
    
    ## assign class and return
    class(obj) <- "dynaTree"
    invisible(obj)
  }


## update:
##
## update the particle cloud to include new (x,y) pairs

update.dynaTree <- function(object, X, y, verb=10, ...)
  {
    ## sanity check
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## extract vitals of X
    X <- matrix(as.matrix(X), ncol=object$m)
    m <- ncol(X)
    T <- nrow(X)

    ## ensure new class labels are valid
    if(object$model == "class") {
      y <- round(y)-1 ## check for sain class labels
      if(any(y < 0)) stop("class labels must start at 1")
      if(any(y > max(object$y)))
        stop("class labels must be <= max(object$y)")
    }
    
    ## echo something to the screen
    if(verb > 0 && T > verb) cat("updating with", T, "new (x,y) pairs\n");

    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## call the C routine to build up the PL object
    object2 <- .C("update_R",
               cloud = as.integer(object$num),
               m = as.integer(m),
               T = as.integer(T),
               X = as.double(t(X)),
               y = as.double(y),
               verb = as.integer(verb),
               lpred = double(T),
               PACKAGE = "dynaTree")

    ## update time
    object$time <- object$time + proc.time()[3] - p1
    
    ## remove cloud
    object2$cloud <- NULL
    
    ## put non-transposed X back, and combine
    if(object$mode == "class") y <- y + 1
    object$X <- rbind(object$X, X)
    object$y <- c(object$y, y)
    object$pred <- c(object$lpred, object2$lpred)

    ## assign class and return
    class(object) <- "dynaTree"
    invisible(object)
  }


## delete.cloud:
##
## deletes the C-side in a particular 

deletecloud <- function(obj)
  {
    if(is.null(obj$num)) stop("NULL cloud obj$num")
    .C("delete_cloud_R",
       num = as.integer(obj$num),
       PACKAGE = "dynaTree")
    invisible(NULL)
  }


## delete.clouds:
##
## deletes all dynatree clouds on the C side

deleteclouds <- function()
  {
    .C("delete_clouds_R", PACKAGE="dynaTree")
    invisible(NULL)
  }


## predict.dynaTree:
##
## generic prediction function for new data XX --
## uses the existing obj$num C-side cloud which
## must not have been deleted

predict.dynaTree <- function(object, XX, quants=TRUE, alc=FALSE, ei=FALSE,
                             verb=0, ...)
  {
    ## extract the vitals of XX
    XX <- as.matrix(XX)
    nn <- nrow(XX)
    if(ncol(XX) != object$m) stop("XX has bad dimensions");

    ## special handling for classification
    if(object$model == "class") {

      ## helpful prints
      if(alc || ei)
        warning("options alc or ei not relevant for classification")
      
      ## number of classification labels
      nc <- max(object$y)
      
      pred <- .C("predclass_R",
                 cloud = as.integer(object$num),
                 XX = as.double(t(XX)),
                 nn = as.integer(nn),
                 verb = as.integer(verb),
                 p = double(nn*nc),
                 ##pvar = double(nn*nc),
                 entropy = double(nn),
                 PACKAGE = "dynaTree")

      ## combune pred with object
      pred$p <- matrix(pred$p, ncol=nc)
      colnames(pred$p) <- paste("c", 1:nc, sep="")
      object$p <- pred$p
      object$entropy <- pred$entropy
      
    } else { ## typical handling for regression
      
      ## call the C-side predict routine
      pred <- .C("predict_R",
                 cloud = as.integer(object$num),
                 XX = as.double(t(XX)),
                 nn = as.integer(nn),
                 verb = as.integer(verb),
                 mean = double(nn),
                 var = double(nn),
                 q1 = double(nn*quants),
                 q2 = double(nn*quants),
                 alc = double(nn*alc),
                 ei = double(nn*ei),
		 PACKAGE = "dynaTree")
      
      ## combine pred with object
      object$mean <- pred$mean
      object$var <- pred$var
      object$q1 <- pred$q1
      object$q2 <- pred$q2
      object$alc <- pred$alc
      object$ei <- pred$ei

      ## possibly remove alc and/or quants
      if(!ei) object$ei <- NULL
      if(!alc) object$alc <- NULL
      if(!quants) object$q1 <- object$q2 <- NULL
    }

    ## put the un-transposed XX back
    object$XX <- XX

    ## assign class and return
    class(pred) <- "dynaTree"
    invisible(object)
  }


## dynaTrees:
##
## calls dynaTree and then predict R times in order to asses
## the Monte Carlo error of the PL procedure and aggregate the
## predictive distributions of many re-roderings of the data

"dynaTrees" <- function(X, y, N=1000, R=10,
                        model=c("constant", "linear", "class"),
                        alpha=0.95, beta=2, minp=NULL, XX=NULL,
                        plotit=FALSE, rorder=TRUE, verb=10, ...)
  {
    ## use dynaTree and predict by themselves when R = 1
    if(R <= 1) stop("R should be >= 2")

    ## coerse X
    X <- as.matrix(X)
    n <- nrow(X)
    
    ## check rorder
    if(length(rorder) > 1) {
      if(nrow(rorder) != nrow(X) && ncol(rorder) != R)
        stop("bad rorder argument")
      else o <- rorder
    } else o <- apply(matrix(runif(nrow(X)*(R-1)), ncol=R-1), 2, order)
    o <- cbind((1:n), o)

    ## build the first model, based on the original ordering,
    ## and predict at XX
    if(verb > 0) cat("\nround 1:\n")
    obj <- dynaTree(X, y, N, model, alpha, beta, minp, verb)
    if(!is.null(XX)) obj <- predict(obj, XX, ...)
    deletecloud(obj); obj$num <- NULL
    obj$num <- NULL

    ## possibly plot in 1d case
    if(plotit) {
      if(ncol(X) != 1)
        warning("cannot plot when ncol(X) >= 1", immediate.=TRUE)
      if(is.null(XX))
        warning("cannot plot without XX predictive grid", immediate.=TRUE)
      plot(obj)
    }

    ## now do the same ting R-1 more times and combine outputs
    for(r in 2:R) {

      ## build the Rth model on a the random re-ordering
      ## and predict at XX
      if(verb > 0) cat("\nround ",  r, ":\n", sep="")
      obj2 <- dynaTree(X[o[,r],], y[o[,r]], N, model, alpha, beta, minp, verb)
      if(!is.null(XX)) obj2 <- predict(obj2, XX, ...)
      deletecloud(obj2); obj2$num <- NULL

      ## possibly add to the plot in the 1d case
      if(plotit && ncol(X) == 1) plot(obj2, add=TRUE)

      ## combine the PL bits of the object
      obj$N <- obj$N + obj2$N
      obj$lpred <- cbind(obj$lpred, obj2$lpred)

      ## combine times
      obj$time <- c(obj$time, obj2$time)

      ## combine the predictive bits
      if(!is.null(XX)) {
        if(obj$model != "class") {
          ## regression collecting
          obj$mean <- cbind(obj$mean, obj2$mean)
          obj$var <- cbind(obj$var, obj2$var)
          obj$q1 <- cbind(obj$q1, obj2$q1)
          obj$q2 <- cbind(obj$q2, obj2$q2)
          obj$alc <- cbind(obj$alc, obj2$alc)
        } else { ## classification averaging
          obj$p <- ((r-1)*obj$p + obj2$p)/r
          obj$entropy <- cbind(obj$entropy, obj2$entropy)
        }
      }
    }
    if(verb > 0) cat("\n")

    ## assign R and class and return
    obj$R <- R
    class(obj) <- "dynaTree"
    invisible(obj)
  }


## alcalc:
##
## extract the active learning heuristic over
## the XX space

alcalc <- function(obj, method=c("alm", "alc", "ei"), prec=1)
  {
    ## check the method argument
    method <- match.arg(method)

    if(method == "alc") return(obj$alc)
    else if(method == "ei")  return(obj$ei + sqrt(obj$var)/prec)
    else return(obj$var)
  }


## getBF:
##
## get the (approximate) Bayes Factor for obj1 and obj2
## setting early non-comparible time steps to zero

getBF <- function(obj1, obj2)
{
  ## sanity check
  if(obj1$R != obj2$R) stop("obj1$R != obj2$R")
  if(nrow(obj1$X) != nrow(obj2$X))
    stop("number of rows does not match")

  ## extract lpred and build BF (for repeats)
  R <- obj1$R
  bf <- matrix(0, nrow=nrow(obj1$X), ncol=R)
  lp1 <- matrix(obj1$lpred, ncol=R)
  lp2 <- matrix(obj2$lpred, ncol=R)

  ## extract the BF for each repeat
  for(r in 1:R) {
    l1 <- which(lp1[,r] != 0)[1] - 1
    l2 <-  which(lp2[,r] != 0)[1] - 1
    l <- max(l1, l2)
    rem <- -(1:l)
    bf[rem,r] <- cumsum(lp1[rem,r]) - cumsum(lp2[rem,r])
  }

  ## return the bayes factor
  return(bf)
}
