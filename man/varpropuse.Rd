\name{varpropuse}
\alias{varpropuse}
\alias{varproptotal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the proportion of variables used in tree splits
}
\description{
  Calculates the proportion of particles which use each input
  to make a tree split and the proportion of all splits in trees
  of each particle that correspond to each input variable
}
\usage{
varpropuse(object)
varproptotal(object)
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
  for this sort of variable selection
}
\value{
  A vector of proportions of length \code{ncol(object$X))} is returned
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

\seealso{ \code{\link{dynaTree}} }
\examples{
## ffit a dynaTree model to the Ozone data
X <- airquality[,2:4]
y <- airquality$Ozone
na <- apply(is.na(X), 1, any) | is.na(y)
out <- dynaTree(X=X[!na,], y=y[!na])

## obtain variable usage proportions
varpropuse(out)
varproptotal(out)

## clean up
deletecloud(out)
}
\keyword{ method }
\keyword{ tree }
