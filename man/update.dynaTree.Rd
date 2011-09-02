\name{update.dynaTree}
\alias{update.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Updating a Dynamic Tree Model With New Data
}
\description{
Updating an already-initialized dynamic tree model with
new input/output pairs, primarily to facilitate sequential
design and optimization applications
}
\usage{
\method{update}{dynaTree}(object, X, y, verb = round(length(y)/10), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{X}{
    an augmenting design \code{matrix} of real-valued predictors
    with \code{ncol(X) = object$m}
}
  \item{y}{
    an augmenting vector of real-valued responses or integer
    categories with \code{length(y) = nrow(X)}
}
  \item{verb}{
    a positive scalar integer indicating how many time steps
    (iterations) after which a progress statement should be
    printed to the console; a value of \code{verb = 0} is quiet
}
    \item{\dots}{
    to comply with the generic \code{\link{predict}} method --
    currently unused
}
}
\details{
  This function updates the \code{\link{dynaTree}} fit with
  new \code{(X,y)} pairs by the Particle Learning (PL)
  algorithm.  The updated fit will be for data combined
  as \code{rbind(object$X, X)} and \code{c(object$y, y)}.

  The primary use of this function is to facilitate sequential
  design by optimization and active learning.  Typically one
  would use \code{\link{predict.dynaTree}} to estimate active
  learning statistics at candidate location.
  These are used to pick new \code{(X,y)}
  locations to add to the design -- the new fit being facilitated
  by this function; see the examples below
}
\value{
  The returned list is the same as \code{\link{dynaTree}} --
  i.e., a \code{"dynaTree"}-class object
}
\references{
 Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1636

  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing}.
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.

  \url{http://faculty.chicagobooth.edu/robert.gramacy/dynaTree.html}
}
\author{
Robert B. Gramacy \email{rbgramacy@chicagobooth.edu} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}
}
\note{
  The object (\code{object}) must contain a pointer to a particle
  cloud (\code{object$num}) which has not been deleted by
  \code{\link{deletecloud}}.  In particular, it cannot be
  an object returned from \code{\link{dynaTrees}}
}

\seealso{
  \code{\link{predict.dynaTree}}, \code{\link{dynaTree}},
  \code{\link{plot.dynaTree}}, \code{\link{deletecloud}},
  \code{\link{getBF}}
}
\examples{
## simple function describing (x,y) data
f1d <- function(x, sd=0.1){
  return( sin(x) - dcauchy(x,1.6,0.15) + rnorm(1,0,sd))
} 

## initial (x,y) data
X <- seq(0, 7, length=30)
y <- f1d(X)

## PL fit to initial data
obj <- dynaTree(X=X, y=y, N=1000, model="linear")

## a predictive grid
XX <- seq(0,7, length=100)
obj <- predict(obj, XX, quants=FALSE)

## follow the ALM algorithm and choose the next
## point with the highest predictive variance
m <- which.max(obj$var)
xstar <- drop(obj$XX[m,])
ystar <- f1d(xstar)

## plot the next chosen point
par(mfrow=c(2,1))
plot(obj, ylab="y", xlab="x", main="fitted surface")
points(xstar, ystar, col=3, pch=20)
plot(obj$XX, sqrt(obj$var), type="l", xlab="x",
     ylab="predictive sd", main="active learning")

## update the fit with (xstar, ystar)
obj <- update(obj, xstar, ystar)

## new predictive surface
obj <- predict(obj, XX, quants=FALSE)

## plotted
plot(obj, ylab="y", xlab="x", main="updated fitted surface")
plot(obj$XX, sqrt(obj$var), type="l", xlab="x",
     ylab="predictive sd", main="active learning")

## delete the cloud to prevent a memory leak
deletecloud(obj); obj$num <- NULL

## see demo("design") for more iterations and
## design under other active learning heuristics
## like ALC, and EI for optimization; also see
## demo("online") for an online learning example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ design }
