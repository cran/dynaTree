\name{alcX.dynaTree}
\alias{alcX.dynaTree}
\alias{alc.dynaTree}
\alias{entropyX.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the ALC or predictive entropy statistic at the X
  locations, or ALC at new XX predictive locations
}
\description{
  Uses analytic integration (at the leaves) to calculate
  the (regression) ALC statistic, or calculate the predictive (class) entropy
 at the input (X) locations; or calculate ALC at new predictive
 locations either analytically or numerically
}
\usage{
\method{alcX}{dynaTree}(object, rect = NULL, categ = NULL,
     approx = FALSE, verb = 0)
\method{entropyX}{dynaTree}(object, verb = 0)
\method{alc}{dynaTree}(object, XX, rect = NULL, categ = NULL,
     approx = FALSE, Xref = NULL, probs = NULL, verb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
  }
  \item{rect}{
    for \code{alcX}, a \code{matrix} with two columns and
    \code{ncol(object$X)} rows
    describing the bounding rectangle for the ALC integration; the
    default that is used when \code{rect = NULL} is the bounding
    rectangle obtained by applying \code{\link{range}} to each
    column of \code{object$X} (taking care to remove the
    first/intercept column of \code{object$X} if \code{icept =
      "augmented"}; only applies to regression models 
      (\code{object$model != "class"}); for \code{alc}, \code{rect}
      must be a scalar logical: see \code{Xref} below
    }
  \item{categ}{
    A vector of logicals of length \code{ncol(object$X)} indicating
    which, if any, dimensions of the input space should be treated
    as categorical; this input is used to helps with the analytic
    integrals from a \code{rect}-based calculation, which means
    it should not specified along with \code{Xref}; the default \code{categ}
    argument is \code{NULL} meaning that the categorical inputs
    are derived from \code{object$X} in a sensible way
  }
  \item{approx}{ a scalar logical that, when \code{TRUE}, causes
    the number of data points in a node/leaf to be used as a proxy for
    its area in the analytic calculations }
  \item{XX}{
    a design \code{matrix} of predictive locations (where \code{ncol(XX) ==
      ncol(X)}; only used by \code{alc} 
 }
 \item{Xref}{
    \code{Xref} input can be optionally used to
    specify a grid of reference locations for the numerical ALC
    calculation - a \code{matrix} with \code{ncol(X)} columns.  If
    \code{NULL}, the default, then the \code{XX} is taken as both
    candidate and reference locations.
}
\item{probs}{
   weights for the reference locations to be used in a Monte Carlo
  approximation; usually these weights are class probabilities for
  response surfaces under constraints
}
  \item{verb}{
    a positive scalar integer indicating how many predictive locations
    (iterations) after which a progress statement should be
    printed to the console; a (default) value of \code{verb = 0} is quiet
  }
}
\details{
  This function is most useful for selecting \code{object$X}
  locations to remove from the analysis, perhaps in an online inference
  setting. See \code{\link{retire.dynaTree}} for more details.  The
  output is the same as using \code{\link{predict.dynaTree}}
  using \code{XX = object$X}, \code{alc = "rect"}, and \code{Xref =
    rect}

   \code{entropyX} only apples to classification models
   (\code{object$model != "class"}), and \code{alcX} applies (only)
   to the other, regression, models

   The \code{alc} function is more generic and allows ALC calculations
   at new, predictive, \code{XX} locations.  This functionality used
   to be part of the \code{\link{predict.dynaTree}} function, but were
   separated out for computational reasons.  The previous version was
   Monte Carlo-based (using \code{Xref}) whereas the new version
   also allows analytic calculation (now the default, via \code{rect})
}
\value{
  The entire \code{object} is returned with a new entry called
  \code{alcX} containing a vector of length \code{nrow(X)} with
  the ALC values, or \code{entropyX} containing the entropy values,
  or \code{alc} if general ALC calculations at new \code{XX} locations
}
\references{
  Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1636
  
  \url{http://faculty.chicagobooth.edu/robert.gramacy/dynaTree.html}
}
\author{
 Robert B. Gramacy \email{rbgramacy@chicagobooth.edu} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}}

\seealso{
  \code{\link{dynaTree}}, \code{\link{predict.dynaTree}}, and
  \code{\link{retire.dynaTree}}
}
\examples{
## fit the model to the parabola data
n <- 100
Xp <- runif(n,-3,3)
Yp <- Xp + Xp^2 + rnorm(n, 0, .2)
rect <- c(-3,3)
out <- dynaTree(Xp, Yp, model="linear", icept="augmented")

## calculate the alcX
out <- alcX(out, rect=rect)

## to compare to analytic 
out <- alc(out, XX=out$X[,-1], rect=rect)

## plot comparison between alcX and predict-ALC
plot(out$X[,-1], out$alcX)
o <- order(out$X[,2])
lines(out$X[o,-1], out$alc[o], col=2, lty=2)

## now compare to approximate analytic
## (which writes over out$alc)
out <- alc(out, XX=out$X[,-1], rect=rect, approx=TRUE)
lines(out$X[o,-1], out$alc[o], col=3, lty=3)

## clean up
deletecloud(out)

## similarly with entropyX for classification models

## see demo("design") for more iterations and
## design under other active learning heuristics
## like ALC, and EI for optimization; also see
## demo("online") for an online learning example where
## ALC is used for retirement
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
