\name{relevance.dynaTree}
\alias{relevance.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate relevance statistics for input coordinates
}
\description{
  Computes relevance statistics for each input coordinate by
  calculating their particle-averaged mean reduction in variance
  each time that coordinate is used as a splitting variable in
  (an internal node of) the tree(s)
}
\usage{
relevance.dynaTree(object, rect = NULL, categ = NULL,
     approx = FALSE, verb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{rect}{
    an optional \code{matrix} with two columns and
    \code{ncol(object$X)} rows describing the bounding rectangle
    for the ALC integration; the
    default that is used when \code{rect = NULL} is the bounding
    rectangle obtained by applying \code{\link{range}} to each
    column of \code{object$X} (taking care to remove the
    first/intercept column of \code{object$X} if \code{icept =
      "augmented"}}
  \item{categ}{
    A vector of logicals of length \code{ncol(object$X)} indicating
    which, if any, dimensions of the input space should be treated
    as categorical; the default \code{categ}
    argument is \code{NULL} meaning that the categorical inputs
    are derived from \code{object$X} in a sensible way}
  \item{approx}{
    a scalar logical indicating if the count of the number of
    data points in the leaf should be used in place of its area;
    this can help with numerical accuracy in high dimensional input
    spaces
  }
  \item{verb}{
    a positive scalar integer indicating how many particles should
    be processed (iterations) before a progress statement should be
    printed to the console; a (default) value of \code{verb = 0} is quiet
  }
}
\details{
  Each binary split in the tree (in each particle) emits a reduction
  in variance (for regression models) or a reduction in entropy
  (for classification).  This function calculates these reductions
  and attributes them to the variable(s) involved in the split(s).
  Those with the largest relevances are the most useful for prediction.
  A sensible variable selection rule based on these relevances is to
  discard those variables whose median relevance is not positive.  See
  the Gramacy, Taddy, \& Wild (2011) reference below for more details.
  
  The new set of particles is appended to the old set.  However
  after a subsequent \code{\link{update.dynaTree}} call the total
  number of particles reverts to the original amount.

  Note that this does not work well with \code{\link{dynaTree}} objects
  which were built with \code{model="linear"}.  Rather, a full
  sensitivity analysis (\code{\link{sens.dynaTree}}) is needed.  Usually
  it is best to first do \code{model="constant"} and then use
  \code{\link{relevance.dynaTree}}.  Bayes factors (\code{\link{getBF}})
  can be used to back up any variable selections implied by the
  relevance.  Then, if desired, one can re-fit on the new (possibly
  reduced) set of predictors with \code{model="linear"}.

  There are no caveats with \code{model="class"}
}
\value{
  The entire \code{object} is returned with a new entry called
  \code{relevance} containing a \code{matrix} with \code{ncol(X)}
  columns.  Each row contains the sample from the relevance of
  each input, and there is a row for each particle
}
\references{
  Gramacy, R.B., Taddy, M.A., and S. Wild (2011).
 \dQuote{Variable Selection and Sensitivity Analysis via
   Dynamic Trees with an Application to Computer Code Performance Tuning}
  arXiv:1108.4739
  
  \url{http://bobby.gramacy.com/r_packages/dynaTree/}
}
\author{
 Robert B. Gramacy \email{rbg@vt.edu}, \cr
Matt Taddy \email{taddy@chicagobooth.edu}, and \cr
Christoforos Anagnostopoulos \email{christoforos.anagnostopoulos06@imperial.ac.uk}}

\seealso{
  \code{\link{dynaTree}}, \code{\link{sens.dynaTree}},
   \code{\link{predict.dynaTree}}
  \code{\link{varpropuse}}, \code{\link{varproptotal}}
}
\examples{
## see the examples in sens.dynaTree for the relevances;
## Also see varpropuse and the class2d demo via
## demo("class2d")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ htest }
