\name{retire.dynaTree}
\alias{retire.dynaTree}

\title{
  Retire (i.e. remove) data from the a dynaTree model
}
\description{
  Allows the removal (or \dQuote{retireing}
  of \code{X}-\code{y} pairs from a
  \code{"dynaTree"}-class object to facilitate online
  learning; \dQuote{retireed} pairs ar absorbed into
  the leaf prior(s)
}
\usage{
\method{retire}{dynaTree}(object, indices, lambda = 1, verb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
  }
  \item{indices}{
    a vector of positive integers in \code{1:nrow(object$X)} indicating
    which \code{X}-\code{y} pairs to \dQuote{retire}; must
    have \code{length(indices) <= nrow(object$X)}
}
\item{lambda}{
   a scalar proportion (forgetting factor) used to downweight the previous prior
  summary statistics
}
  \item{verb}{
    a nonzero scalar causes info about the \dQuote{retireed} indices,
    i.e., their \code{X}-\code{y} values, to be printed to the
    screen as they are \dQuote{retireed}
}
}
\details{
  Primarily for use in online learning contexts.  After
  \dQuote{retireing} the predictive distribution remains unchanged,
  because the sufficient statistics of the removed pairs enters
  the prior in the leaves of the tree of each particle.  Further
  \code{\link{update.dynaTree}} calls (adding data) may cause
  changes to the posterior predictive as grow moves cannot keep
  the \dQuote{retires}; see a forthcoming paper for more
  details.  In many ways, \code{\link{retire.dynaTree}} is the
  opposite of \code{\link{update.dynaTree}} except that the loss of
  information upon \dQuote{retireing} is not complete.

  Drifting regression or classification relationships may be modeled
  with a forgetting factor \code{lambda < 1}

  The \code{\link{alcX.dynaTree}} provides a good, and computationally
  efficient, heuristic for choosing which points to \dQuote{retire} for 
  regression models, and likewise \code{link{entropyX.dynaTree}} for 
  classification models.

  Note that classification models (\code{model = "class"}) are
  not supported, and implicit intercepts (\code{icept = "implicit"})
  with linear models (\code{model = "linear"}) are not supported
  at this time
}
\note{ In order to use \code{model = "linear"} with
  \code{\link{dynaTree}} and retirement one must also specify
  \code{icept = "augmented"} which automatically augments an
  extra column of ones onto the input \code{X} design matrix/matrices.
  The \code{retire} function only supports this \code{icept} case
}
\value{
  returns a \code{"dynaTree"}-class object with updated attributes
}
\references{
  Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1586 

  \url{http://faculty.chicagobooth.edu/robert.gramacy/dynaTree.html}
}
\author{
Robert B. Gramacy \email{rbgramacy@chicagobooth.edu} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dynaTree}}, \code{\link{alcX.dynaTree}}, 
  \code{\link{entropyX.dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{rejuvenate.dynaTree}}
}
\examples{
n <- 100
Xp <- runif(n,-3,3)
XX <- seq(-3,3, length=200)
Yp <- Xp + Xp^2 + rnorm(n, 0, .2)
rect <- c(-3,3)
out <- dynaTree(Xp, Yp, model="linear", icept="augmented")

## predict and plot
out <- predict(out, XX)
plot(out, main="parabola data", lwd=2)

## randomly remove half of the data points
out <- retire(out, sample(1:n, n/2, replace=FALSE))

## predict and add to plot -- shouldn't change anything
out <- predict(out, XX)
plot(out, add=TRUE, col=3)
points(out$X[,-1], out$y, col=3)

## now illustrating rejuvenation, which should result
## in a change to the predictive surface
out <- rejuvenate(out)
out <- predict(out, XX)
plot(out, add=TRUE, col=4)
legend("top", c("original", "retired", "rejuvenated"),
       col=2:4, lty=1)

## clean up
deletecloud(out)

## see demo("online") for an online learning example
## where ALC is used for retirement
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
