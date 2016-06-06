\name{varpropuse}
\alias{varpropuse.dynaTree}
\alias{varproptotal.dynaTree}
\alias{treestats.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the proportion of variables used in tree splits, and
  average summary stats of tree heights and leaf sizes
}
\description{
  Calculates the proportion of particles which use each input
  to make a tree split and the proportion of all splits in trees
  of each particle that correspond to each input variable; also
  provides tree height and leaf size summary information
}
\usage{
\method{varpropuse}{dynaTree}(object)
\method{varproptotal}{dynaTree}(object)
\method{treestats}{dynaTree}(object)
}
\arguments{
  \item{object}{
     a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
   }
}
\details{
  \code{varpropuse} gives the proportion of times a particle
  uses each input variable in a tree split; \code{varproptotal} gives
  the proportion of total uses by the tree in each particle (i.e.,
  averaged over the total number of splits used in the tree).

  Usually, \code{varpropuse} returns a vector of (nearly) all ones
  unless there are variables which are not useful in predicting
  the response.  Using \code{model = "linear"} is not recommended
  for this sort of variable selection.

  \code{treestats} returns the average tree height, and the average
  leaf size, both active and retired
}
\value{
  For \code{varprop*}, a
  vector of proportions of length \code{ncol(object$X))} is returned;
  for \code{treestats} a 1-row, 4-column \code{\link{data.frame}} is
  returned
}
\references{
  Gramacy, R.B., Taddy, M.A., and S. Wild (2011).
 \dQuote{Variable Selection and Sensitivity Analysis via
   Dynamic Trees with an Application to Computer Code Performance Tuning}
  arXiv:1108.4739
  
  \url{http://bobby.gramacy.com/r_packages/dynaTree/}
}
\author{
  Robert B. Gramacy \email{rbgramacy@chicagobooth.edu}, \cr
Matt Taddy \email{taddy@chicagobooth.edu}, and \cr
Christoforos Anagnostopoulos \email{christoforos.anagnostopoulos06@imperial.ac.uk}}

\seealso{ \code{\link{dynaTree}}, \code{\link{sens.dynaTree}},
\code{\link{relevance.dynaTree}} }
\examples{
## ffit a dynaTree model to the Ozone data
X <- airquality[,2:4]
y <- airquality$Ozone
na <- apply(is.na(X), 1, any) | is.na(y)
out <- dynaTree(X=X[!na,], y=y[!na])

## obtain variable usage proportions
varpropuse(out)
varproptotal(out)

## gather relevance statistics which are more meaningful
out <- relevance(out)
boxplot(out$relevance)
abline(h=0, col=2, lty=2)

## obtain tree statistics
treestats(out)

## clean up
deletecloud(out)
}
\keyword{ method }
\keyword{ tree }
