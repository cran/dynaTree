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


## plot.dynaTree:
##
## plot the contents of a single dynaTree object that
## contains predictive information

plot.dynaTree <-
  function(x, proj=1, add=FALSE, ylim=NULL, col=2, lwd=1,
           ptype=c("each", "mean"), ...)
  {
    ## sanity check
    if(is.null(x$XX)) stop("nothing to plot; call predict first")

    ## match ptype
    ptype = match.arg(ptype)

    ## check proj
    if(length(proj) != 1 || proj <= 0 || proj > ncol(x$X))
      stop("bad proj")
    if(x$icept == "augmented") proj <- proj + 1
    
    ## extract projected X and XX
    X <- x$X[,proj]; XX <- x$XX[,proj]
    
    ## usually plot the data first
    if(add==FALSE) {
      if(is.null(ylim)) {
        if(!is.null(x$q1)) ylim <- range(c(x$q1, x$q2))
        else ylim <- range(x$y)
      }
      plot(X, x$y, pch=20, ylim=ylim, ...)
    }

    ## for plotting in orer
    o <- order(XX)

    ## plot a single set of (3) lines
    if(is.null(x$R)) {
      lines(XX[o], x$mean[o], col=col, lwd=lwd)
      if(!is.null(x$q1)) {
        lines(XX[o], x$q1[o], col=col, lty=2, lwd=lwd)
        lines(XX[o], x$q2[o], col=col, lty=2, lwd=lwd)
      }
    } else { ## or plot R sets

      if(ptype == "each") {
        ## plotting each pair of lines
        for(r in 1:x$R) {
          lines(XX[o], x$mean[o,r], col=col, lwd=lwd)
          if(!is.null(x$q1)) {
            lines(XX[o], x$q1[o,r], col=col, lty=2, lwd=lwd)
            lines(XX[o], x$q2[o,r], col=col, lty=2, lwd=lwd)
          }
        }
      } else if(ptype == "mean") {
        ## plotting the average of the lines
        lines(XX[o], apply(x$mean, 1, mean)[o],
              col=col, lwd=lwd)
        if(!is.null(x$q1)) {
          lines(XX[o], apply(x$q1, 1, mean)[o],
                col=col, lty=2, lwd=lwd)
          lines(XX[o], apply(x$q2, 1, mean)[o],
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
                         track, f1d, prec, ...)
{  
  ## set up to plot
  par(mfrow=c(2,2), mai=c(.8,.7,.2,.1))
  
  ## plot current data and newly chosen point
  plot(x, ylab="y", xlab="x", ...)
  points(xstar, ystar, col=3, pch=20)
  abline(v=xstar, col=3, lty=2)
  
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
  plotal(x, method=method, prec=prec, xlab="x")
  abline(v=xstar, col=3, lty=2)

  ## histogram of X-samples
  hist(x$X, main="", xlab="X samples")
}



## plotal:
##
## plot the active learning heuristic --
## ised nu plotprogress

plotal <- function(x, method=c("alm", "alc", "ei", "ieci", "qEntropy", "qEI"),
                   prec=1, add=FALSE, each=FALSE, root=FALSE, ylim=NULL, ...)
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



