## plot.dynaTree:
##
## plot the contents of a single dynaTree object that
## contains predictive information

plot.dynaTree <-
  function(x, add=FALSE, ylim=NULL, col=2, lwd=1,
           ptype=c("each", "mean"), ...)
  {
    ## sanity check
    if(is.null(x$XX)) stop("nothing to plot; call predict first")

    ## match ptype
    ptype = match.arg(ptype)

    ## usually plot the data first
    if(add==FALSE) {
      if(is.null(ylim)) {
        if(!is.null(x$q1)) ylim <- range(c(x$q1, x$q2))
        else ylim <- range(x$y)
      }
      plot(x$X, x$y, pch=20, ylim=ylim, ...)
    }

    ## for plotting in orer
    o <- order(x$XX)

    ## plot a single set of (3) lines
    if(is.null(x$R)) {
      lines(x$XX[o], x$mean[o], col=col)
      if(!is.null(x$q1)) {
        lines(x$XX[o], x$q1[o], col=col, lty=2)
        lines(x$XX[o], x$q2[o], col=col, lty=2)
      }
    } else { ## or plot R sets

      if(ptype == "each") {
        ## plotting each pair of lines
        for(r in 1:x$R) {
          lines(x$XX[o], x$mean[o,r], col=col)
          if(!is.null(x$q1)) {
            lines(x$XX[o], x$q1[o,r], col=col, lty=2)
            lines(x$XX[o], x$q2[o,r], col=col, lty=2)
          }
        }
      } else if(ptype == "mean") {
        ## plotting the average of the lines
        lines(x$XX[o], apply(x$mean, 1, mean)[o],
              col=col, lwd=lwd)
        if(!is.null(x$q1)) {
          lines(x$XX[o], apply(x$q1, 1, mean)[o],
                col=col, lty=2, lwd=lwd)
          lines(x$XX[o], apply(x$q2, 1, mean)[o],
                col=col, lty=2, lwd=lwd)
        }
      } 
    }
  }



## plotprogress:
##
## function to plot the progress of the sequential
## design algorithm in one dimension

plotprogress <- function(x, xstar, ystar, method="alm",
                         track, f1d, ...)
{  
  ## set up to plot
  par(mfrow=c(2,2), mai=c(.8,.7,.2,.1))
  
  ## plot current data and newly chosen point
  plot(x, ylab="y", xlab="x", ...)
  points(xstar, ystar, col=3, pch=20)
  
  ## plot the truth
  o <- order(x$XX)
  lines(x$XX[o], f1d(x$XX[o], sd=0), col=8,type="l")

  ## plot progress tracked over time
  n <- nrow(x$X)
  start <- n - length(track) +1
  plot(start:n, log(track), type="l", xlim=c(1,n),
       ylab=paste("log", method), xlab="time")
  abline(v=start-0.5, col=2, lty=2)

  ## plot active learning statistic
  plotal(x, method=method, xlab="x")

  ## histogram of X-samples
  hist(x$X, main="", xlab="X samples")
}



## plotal:
##
## plot the active learning heuristic --
## ised nu plotprogress

plotal <- function(x, method=c("alm", "alc", "ei"), prec=1,
                   add=FALSE, each=FALSE, root=FALSE, ylim=NULL, ...)
  {
    ## complete the method argument
    method <- match.arg(method)

    ## extract the variance (or ALC or EI) to be plotted
    v <- alcalc(x, method, prec)

    ## shall we plot the sqrt?
    if(root) v <- sqrt(v)
    if(is.null(ylim)) ylim <- range(v)

    ## coerse to a matrix
    if(is.null(dim(v))) v <- matrix(v, ncol=1)

    ## get the order of the XX
    o <- order(x$XX)
    
    if(each) {
      ## extract R and check
      R <- x$R
      if(is.null(R)) stop("no multiple runs in xect")

      ## plot the first variance
      if(add) lines(x$XX[o], v[o,1], pch=20, type="l", ...)
      else plot(x$XX[o], v[o,1], pch=20, ylim=ylim, type="l",
                ylab=method, ...)      

      ## plot the rest
      for(r in 2:R) lines(x$XX[o], v[o,r], col=2)
      
    } else { ## plot the variance
      if(add) lines(x$XX[o], apply(v, 1, mean)[o], ...)
      else plot(x$XX[o], apply(v, 1, mean)[o], type="l",
                ylab=method, ...)
    }
  }



