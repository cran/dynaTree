## predict.dynaTree:
##
## generic prediction function for new data XX --
## uses the existing obj$num C-side cloud which
## must not have been deleted

predict.dynaTree <- function(object, XX, yy=NULL, quants=TRUE,
                             ei=FALSE, verb=0, ...)
  {
    ## for timing purposes
    ## p1 <- proc.time()[3]

    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## extract the vitals of XX
    XX <- as.matrix(XX)
    nn <- nrow(XX)
    if(object$icept == "augmented") XX <- cbind(rep(1,nn), XX)
    if(ncol(XX) != object$m) stop("XX has bad dimensions");

    ## check yy if not NULL
    if(!is.null(yy)) {
      if(length(yy) != nn) stop("length(yy) must match nrow(XX)")
      if(object$model == "class") {
        yy <- round(yy)-1 ## check for sain class labels
        if(any(yy < 0)) stop("class labels must start at 1")
        if(any(yy >= max(object$classes)))
          stop("class labels must be <= max(object$y)")
      }
    }

    ## special handling for classification
    if(object$model == "class") {

      ## helpful prints
      if(ei) warning("ei not relevant for classification")
      
      ## number of classification labels
      nc <- max(object$classes)
      
      pred <- .C("predclass_R",
                 cloud = as.integer(object$num),
                 XX = as.double(t(XX)),
                 yy = as.integer(yy),
                 nn = as.integer(nn),
                 verb = as.integer(verb),
                 p = double(nn*nc),
                 yypred = double(length(yy)),
                 entropy = double(nn),
                 PACKAGE = "dynaTree")

      ## combine pred with object
      pred$p <- matrix(pred$p, ncol=nc)
      colnames(pred$p) <- paste("c", 1:nc, sep="")
      object$p <- pred$p
      object$entropy <- pred$entropy
      object$yypred <- pred$yypred
      
    } else { ## typical handling for regression
      
      ## call the C-side predict routine
      pred <- .C("predict_R",
                 cloud = as.integer(object$num),
                 XX = as.double(t(XX)),
                 yy = as.double(yy),
                 nn = as.integer(nn),
                 verb = as.integer(verb),
                 mean = double(nn),
                 var = double(nn),
                 q1 = double(nn*quants),
                 q2 = double(nn*quants),
                 yypred = double(length(yy)),
                 ei = double(nn*ei),
		 PACKAGE = "dynaTree")
      
      ## combine pred with object
      object$mean <- pred$mean
      object$var <- pred$var
      object$q1 <- pred$q1
      object$q2 <- pred$q2
      object$yypred <- pred$yypred
      object$ei <- pred$ei

      ## possibly remove ei and/or quants
      if(!ei) object$ei <- NULL
      if(!quants) object$q1 <- object$q2 <- NULL
    }

    ## put originals back
    object$XX <- XX

    ## update time
    ## object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }


## classprobs.dynaTree:
##
## return the probabilities of a particular class, or a Boolean
## sample of that class

setGeneric("classprobs",
            function(object, ...)
            standardGeneric("classprobs")
            )

classprobs.dynaTree <- function(object, class, XX, output=c("p", "c", "both"),
                                verb=0, ...)
  {
    ## for timing purposes
    ## p1 <- proc.time()[3]

    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## sanity check
    if(object$model != "class") stop("only valid for classification models")

    ## number of classification labels
    nc <- max(object$classes)
    if(length(class) != 1 || class <= 0 || class > nc)
      stop("class must be in 1:nc=max(object$classes)")

    ## extract the vitals of XX
    XX <- as.matrix(XX)
    nn <- nrow(XX)
    if(object$icept == "augmented") XX <- cbind(rep(1,nn), XX)
    if(ncol(XX) != object$m) stop("XX has bad dimensions");
    
    ## check output argument
    output <- match.arg(output)
    pcls <- cls <- FALSE
    if(output == "p") pcls <- TRUE
    else if (output == "c") cls <- TRUE
    else pcls <- cls <- TRUE
      
    pred <- .C("classprobs_R",
               cloud = as.integer(object$num),
               class = as.integer(class - 1),
               XX = as.double(t(XX)),
               nn = as.integer(nn),
               verb = as.integer(verb),
               pcls = double(nn*object$N*pcls),
               cls = double(nn*object$N*cls),
               PACKAGE = "dynaTree")

    ## combine pred with object
    if(pcls) object$pcls[[class]] <- matrix(pred$pcls, nrow=object$N, byrow=TRUE)
    if(cls) object$cls[[class]] <- matrix(pred$cls, nrow=object$N, byrow=TRUE)
      
    ## put originals back
    object$XX <- XX

    ## update time
    ## object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }

setMethod("classprobs", "dynaTree", classprobs.dynaTree)


## ieci.dynaTree:
##
## ieci calculation at new XX locations based on reference
## locations in Xref -- uses the existing
## obj$num C-side cloud which must not have been deleted

setGeneric("ieci",
            function(object, ...)
            standardGeneric("ieci")
            )

ieci.dynaTree <- function(object, XX, Xref=NULL, probs=NULL, verb=0)
  {
    ## for timing purposes
    ## p1 <- proc.time()[3]

    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## sanity check
    if(object$model == "class") stop("not for use with classification models")
    
    ## extract the vitals of XX
    XX <- as.matrix(XX)
    nn <- nrow(XX)
    if(object$icept == "augmented") XX <- cbind(rep(1,nn), XX)
    if(ncol(XX) != object$m) stop("XX has bad dimensions");
    
    ## checking Xref and rect
    if(!is.null(Xref)) { ## numerical ALC with Xref

      ## check formatting of Xref
      if(is.null(nrow(Xref))) Xref <- matrix(Xref, ncol=1)
      nref <- nrow(Xref)
      if(object$icept == "augmented") Xref <- cbind(rep(1,nref), Xref)
      if(nref > 0 && ncol(Xref) != object$m)
        stop("Xref has bad dimensions")
      Xref <- t(Xref)

    } else nref <- 0 ## using XX as Xref

    ## check probs
    if(!is.null(probs) && (is.null(dim(probs)) ||
                           nrow(probs) != object$N || ncol(probs) != nref))
      stop("probs must be a object$N * nrow(Xref) matrix")
    if(!is.null(probs)) probs <- t(probs)

    ## call the C-side alc routine
    pred <- .C("ieci_R",
               cloud = as.integer(object$num),
               XX = as.double(t(XX)),
               nn = as.integer(nn),
               Xref = as.double(Xref), 
               nref = as.integer(nref),
               probs = as.double(probs), ## transposed above
               verb = as.integer(verb),
               ieci = double(nn),
               PACKAGE = "dynaTree")
      
    ## combine pred with object
    object$ieci <- pred$ieci

    ## put originals back
    object$XX <- XX
    object$Xref <- Xref
    if(!is.null(probs)) object$probs <- t(probs)

    ## update time
    ## object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }

setMethod("ieci", "dynaTree", ieci.dynaTree)


## alc.dynaTree:
##
## alc calculation at new XX locations based on reference
## locations in Xref, or an analytic calculation using a
## rectangle in Xref (when rect=TRUE) -- uses the existing
## obj$num C-side cloud which must not have been deleted

setGeneric("alc",
            function(object, ...)
            standardGeneric("alc")
            )

alc.dynaTree <- function(object, XX, rect=NULL, Xref=NULL, probs=NULL, verb=0)
  {
    ## for timing purposes
    ## p1 <- proc.time()[3]

    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")

    ## sanity check
    if(object$model == "class") stop("not for use with classification models")
    
    ## extract the vitals of XX
    XX <- as.matrix(XX)
    nn <- nrow(XX)
    if(object$icept == "augmented") XX <- cbind(rep(1,nn), XX)
    if(ncol(XX) != object$m) stop("XX has bad dimensions");
    
    ## checking Xref and rect
    if(!is.null(Xref)) { ## numerical ALC with Xref

      ## sanity check
      if(!is.null(rect)) stop("Xref only valid when rect is NULL")

      ## check formatting of Xref
      if(is.null(nrow(Xref))) Xref <- matrix(Xref, ncol=1)
      nref <- nrow(Xref)
      if(object$icept == "augmented") Xref <- cbind(rep(1,nref), Xref)
      if(nref > 0 && ncol(Xref) != object$m)
        stop("Xref has bad dimensions")
      Xref <- t(Xref)

    } else if(length(rect) == 1 && rect == FALSE) {
      nref <- 0 ## using XX as Xref
    } else {
      ## sanity check
      if(!is.null(probs)) stop("must use Xref to use probs")
      
      ## Xref is null, using rect
      if(is.null(rect)) ## automatic rect
        rect <- t(apply(rbind(object$X, XX), 2, range)) 
      else { ## specified rect
        if(is.null(nrow(rect))) rect <- matrix(rect, nrow=1)
        if(object$icept == "augmented") rect <- rbind(rep(1,2), rect)
      }
      
      ## post sanity check
      if(ncol(rect) != 2 && nrow(rect) != object$m)
        stop("rect has bad dimensions")
      
      ## for sending to C
      nref <- -1; Xref <- rect 
    }

    ## check probs
    if(!is.null(probs) && (is.null(dim(probs)) ||
                           nrow(probs) != object$N || ncol(probs) != nref))
      stop("probs must be a object$N * nrow(Xref) matrix")
    if(!is.null(probs)) probs <- t(probs)

    ## call the C-side alc routine
    pred <- .C("alc_R",
               cloud = as.integer(object$num),
               XX = as.double(t(XX)),
               nn = as.integer(nn),
               Xref = as.double(Xref), 
               nref = as.integer(nref),
               probs = as.double(probs),
               verb = as.integer(verb),
               alc = double(nn),
               PACKAGE = "dynaTree")
      
    ## combine pred with object
    object$alc <- pred$alc

    ## put originals back
    object$XX <- XX
    object$Xref <- Xref
    if(!is.null(probs)) object$probs <- t(probs)

    ## update time
    ## object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }

setMethod("alc", "dynaTree", alc.dynaTree)


## alcX.dynaTree:
##
## calculate the ALC statistic at the data locations
## based on a rectangle of reference locations; 
## uses the existing obj$num C-side cloud which
## must not have been deleted

setGeneric("alcX",
            function(object, ...)
            standardGeneric("alcX")
            )

alcX.dynaTree <- function(object, rect=NULL, verb=0)
  {
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no could number in object")
    
    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## extract/check the vitals of rect
    if(object$model == "class") stop("not for use with classification models")
    if(is.null(rect)) rect <- t(apply(object$X, 2, range)) ## automatic rect
    else { ## specified rect
      if(is.null(nrow(rect))) rect <- matrix(rect, nrow=1)
      if(object$icept == "augmented") rect <- rbind(rep(1,2), rect)
    }
    if(ncol(rect) != 2 && nrow(rect) != object$m)
      stop("rect has bad dimensions");

    ## call the C-side predict routine
    pred <- .C("alcX_R",
               cloud = as.integer(object$num),
               rect = as.double(rect),
               verb = as.integer(verb),
               alcX = double(nrow(object$X)),
               PACKAGE = "dynaTree")
      
    ## combine pred with object
    object$alcX <- pred$alcX
    object$rect <- rect

    ## update time
    object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }

setMethod("alcX", "dynaTree", alcX.dynaTree)


## entropyX.dynaTree:
##
## calculate Entropy statistics at the data locations
## uses the existing obj$num C-side cloud which
## must not have been deleted

setGeneric("entropyX",
            function(object, ...)
            standardGeneric("entropyX")
            )

entropyX.dynaTree <- function(object, verb=0)
  {
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## for timing purposes
    p1 <- proc.time()[3]
    
    ## extract/check the vitals of rect
    if(object$model != "class") stop("only used by classification models")

    ## call the C-side predict routine
    pred <- .C("entropyX_R",
               cloud = as.integer(object$num),
               verb = as.integer(verb),
               entropyX = double(nrow(object$X)),
               PACKAGE = "dynaTree")
      
    ## combine pred with object
    object$entropyX <- pred$entropyX

    ## update time
    object$time <- object$time + proc.time()[3] - p1
    
    ## return
    invisible(object)
  }

setMethod("entropyX", "dynaTree", entropyX.dynaTree)


## sens.dynaTree:
##
## sensity analysis for dynaTree models
## uses the existing obj$num C-side cloud which
## must not have been deleted; a bunch of this C code
## was copied from the tgp package

setGeneric("sens",
            function(object, ...)
            standardGeneric("sens")
            )

sens.dynaTree <- function(object, class=NULL, nns=1000, nME=100, span=0.3,
                           method=c("lhs", "boot"), lhs=NULL, categ=NULL,
                           verb=0)
  {  
    ## make sure object$num is defined
    if(is.null(object$num)) stop("no cloud number in object")
    
    ## process the X data and extract dimensions & intercept
    X <- as.matrix(object$X)
    T <- nrow(X)
    if(object$icept == "augmented") X <- X[,-1]
    d <- ncol(X)
    
    ## check nns and nME
    nns <- round(nns)
    if(length(nns) != 1 || nns <= 0) stop("nns should be a positive scalar integer")
    nME <- round(nME)
    if(length(nME) != 1 || nME <= 0) stop("nME should be a positive scalar integer")
    
    ## process the method argument
    method <- match.arg(method)
    if(method == "boot" && !is.null(lhs)) {
      warning("lhs should be NULL for boot methd")
      lhs <- NULL
    }
    
    ## default categ argument
    if(is.null(categ)) {
      categ <- apply(X, 2, function(x) { setequal(unique(x), c(0,1)) })      
    } else if(verb > 0) { cat("treating as categorical:\n"); print(categ) } 
    
    ## check categ argument (not sure about the second check)
    if(length(categ) != d && ! is.logical(categ))
      stop("categ should be NULL or a length-d logical vector")
    if(method == "boot" && any(categ)) {
      if(any(apply(as.matrix(X[,categ]), 2, 
                   function(x) { !setequal(unique(x), c(0,1)) })))
        stop("X[,categ] can only have entries in {0,1}")
    }
    
    ## deal with LHS list
    if(is.null(lhs)) lhs <- list(rect=NULL, shape=NULL, mode=NULL)
    if(method == "lhs") {
      
      ## process the rect argument
      if(is.null(lhs$rect)) lhs$rect <- apply(X,2,range)
      else if(nrow(lhs$rect) != d || ncol(lhs$rect) != 2)
        stop(paste("rect should be a ", d, "x2-vector", sep=""))
      
      ## check the shape LHS parameter vector
      if(is.null(lhs$shape)) lhs$shape <- as.numeric(!categ)
      else if(length(lhs$shape) != d || !all(lhs$shape >= 0)) { 
        print(lhs$shape)
        stop(paste("lhs$shape should be a non-negative ", d, "-vector", sep=""))
      }
      
      ## check the mode LHS parameter vector
      if(is.null(lhs$mode)) lhs$mode <- apply(X,2,mean)
      else if(length(lhs$mode) != d) {
        print(lhs$mode)
        stop(paste("lhs$mode should be a ", d, "-vector", sep=""))
      }
      
      ## check the LHS rectangle in the categorical variable context
      for(i in 1:d){
        if(categ[i] == TRUE){
          if(lhs$rect[1,i] != 0 || lhs$rect[2,i] != 1){
            print(lhs$rect[i,])
            stop(paste("lhs$rect must be [0,1] for categorical variables (i=",
                       i,", lhs$shape[i]=", lhs$shape[i],").", sep=""))
          }
          if(any(lhs$shape[categ] != 0)) 
            stop("shape must be zero for categorical inputs")
        }
      }
      
      ## transpose rect
      lhs$rect <- t(lhs$rect)
    }
    
    ## check span and create Main Effect grid
    if(length(span) != 1 || ((span > 1) || (span < 0)))
      stop("Bad smoothing span; must be scalar in (0,1).")
    MEgrid <- matrix(ncol=d, nrow=nME)
    if(! is.null(lhs$rect)) MErect <- t(lhs$rect)
    else MErect <- apply(X,2,range)
    for(i in 1:d){ MEgrid[,i] <- seq(MErect[1,i], MErect[2,i], length=nME) }
    
    ## checks for classification
    if(object$model == "class") {
      
      ## default to all classes
      if(is.null(class)) class <- sort(unique(object$y))
      else { ## check the class argument
        for(cls in class) {
          uy <- unique(object$y)
          if(all(cls != uy)) {
            print(sort(uy))
            stop("classes must be one of the above choices")
          }
        }
      }
    } else class <- -1 ## special coding for regression
    
    ## looping over classes, or single run for regression
    for(cls in class) {
      
      ## dummy class for regression
      if(cls < 0) {
        cls <- NULL
        if(verb > 0) cat("regression sensitivity analysis:\n")
      } else if(verb > 0) cat("class", cls, "sensitivity analysis:\n")
      
      ## call the C-side sens routine
      sens <- .C("sens_R",
                 cloud = as.integer(object$num),
                 cls = as.integer(cls-1),
                 nns = as.integer(nns),
                 aug = as.integer(object$m - d), 
                 rect = as.double(lhs$rect),
                 shape = as.double(lhs$shape),
                 mode = as.double(lhs$mode),
                 categ = as.integer(categ),
                 nME = as.integer(nME),
                 span = as.double(span),
                 MEgrid = as.double(MEgrid), ## do not transpose
                 verb = as.integer(verb),
                 MEmean = double(nME*d),
                 MEq1 = double(nME*d),
                 MEq2 = double(nME*d),
                 S = double(d*object$N),
                 T = double(d*object$N),
                 PACKAGE = "dynaTree")
      
      ## check S and T are in the right range
      o <- sens$S > 1 | sens$S < 0
      if(any(o)) sens$S[o] <- NA
      o <- sens$T < 0 | sens$T > 1
      if(any(o)) sens$T[o] <- NA
      
      ## names for columns of the outputs
      dn <- list(NULL, colnames(object$X))
      
      ## combine pred with object
      if(object$model == "class") { ## for classification
        object$MEmean[[cls]] <- matrix(sens$MEmean, ncol=d, dimnames=dn)
        object$MEq1[[cls]] <- matrix(sens$MEq1, ncol=d, dimnames=dn)
        object$MEq2[[cls]] <- matrix(sens$MEq2, ncol=d, dimnames=dn)
        object$S[[cls]] <- matrix(sens$S, ncol=d, byrow=TRUE, dimnames=dn)
        object$T[[cls]] <- matrix(sens$T, ncol=d, byrow=TRUE, dimnames=dn)
      } else { ## for regression
        object$MEgrid <- matrix(sens$MEgrid, ncol=d, dimnames=dn)
        object$MEmean <- matrix(sens$MEmean, ncol=d, dimnames=dn)
        object$MEq1 <- matrix(sens$MEq1, ncol=d, dimnames=dn)
        object$MEq2 <- matrix(sens$MEq2, ncol=d, dimnames=dn)
        object$S <- matrix(sens$S, ncol=d, byrow=TRUE, dimnames=dn)
        object$T <- matrix(sens$T, ncol=d, byrow=TRUE, dimnames=dn)
      }
    } ## end class loop
    
    ## remember which classes sens was used on, if any
    if(all(class != -1)) object$sens.class <- class
    object$MEgrid <- MEgrid
    colnames(object$MEgrid) <- dn[[2]]
    
    ## assign class and return
    class(sens) <- "dynaTree"
    invisible(object)
  }

setMethod("sens", "dynaTree", sens.dynaTree)
    
## varpropuse:
##
## calculate the proportion of particles that use each column
## of X to make a treed partition

varpropuse <- function(object)
  {
    vc <- .C("varpropuse_R",
             cloud = as.integer(object$num),
             counts = double(ncol(object$X)),
             PACKAGE = "dynaTree")$counts

    vc[vc < 0] <- NA;
    if(object$icept == "augmented") vc <- vc[-1]

    return(vc)
  }


## varproptotal:
##
## calculate the proportion of particles that use each column
## of X to make a treed partition

varproptotal <- function(object)
  {
    vc <- .C("varproptotal_R",
             cloud = as.integer(object$num),
             counts = double(ncol(object$X)),
             PACKAGE = "dynaTree")$counts
    
    vc[vc < 0] <- NA;
    if(object$icept == "augmented") vc <- vc[-1]
    
    return(vc)
  }


## alcalc:
##
## extract the active learning heuristic over
## the XX space

alcalc <- function(obj, method=c("alm", "alc", "ei", "ieci"), prec=1)
  {
    ## check the method argument
    method <- match.arg(method)

    if(method == "alc") return(obj$alc)
    else if(method == "ei" && prec > 0) return(obj$ei + sqrt(obj$var)/prec)
    else if(method == "ei") return(obj$ei)
    else if(method == "ieci") return(obj$ieci)
    else return(obj$var)
  }


## getBF:
##
## get the (approximate) Bayes Factor for obj1 and obj2
## setting early non-comparible time steps to zero

getBF <- function(obj1, obj2)
{
  ## one repeat default
  if(is.null(obj1$R)) obj1$R <- 1
  if(is.null(obj2$R)) obj2$R <- 1
  
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