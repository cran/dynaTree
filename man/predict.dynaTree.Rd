\name{predict.dynaTree}
\alias{predict.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for Dynamic Tree Models
}
\description{
Predicting and calculating sequential design and optimization
statistics at new design points (i.e., active learning heuristics)
for dynamic tree models
}
\usage{
\method{predict}{dynaTree}(object, XX, yy = NULL, quants = TRUE,
                 ei = FALSE, verb = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{XX}{
    a design \code{matrix} of predictive locations (where \code{ncol(XX) ==
      ncol(X)})
}
  \item{yy}{
  an optional vector of \dQuote{true} responses at the \code{XX}
  predictive locations at which the log posterior probability are
  to be reported
  }  	 
 \item{quants}{
    a scalar \code{logical} indicating if predictive quantiles
    are desired (useful for visualization, but less so for active
    learning); calculating predictive quantiles is expensive and should
    be turned off if prediction is not being used for visualization,
    e.g., if used for active learning
}
  \item{ei}{
    a scalar \code{logical} indicating if the expected improvement
    statistic (for optimization) should be calculated and returned
}
  \item{verb}{
    a positive scalar integer indicating how many predictive locations
    (iterations) after which a progress statement should be
    printed to the console; a (default) value of \code{verb = 0} is quiet
}
  \item{\dots}{
    to comply with the generic \code{\link{predict}} method --
    currently unused
}
}
\details{
  Returns predictive summary statistics by averaging over the
  samples from the posterior predictive distribution obtained
  from each of the particles in the cloud pointed to by the
  object (\code{object})

  The value(s) calculated are appended to \code{object}; the new
  fields are described below

  ALC calculations have been moved to the \code{\link{alc.dynaTree}}
  function(s)
}
\value{
  The object returned is of class \code{"dynaTree"}, which includes a
  copy of the list elements from the \code{object} passed in, 
  with the following (predictive)
  additions depending on whether \code{object$model} is for
  regression (\code{"constant"} or \code{"linear"}) or classification
  (\code{"class"}).

  For regression:
    
  \item{mean }{ a vector containing an estimate of the predictive mean
    at the \code{XX} locations }
  \item{var }{ a vector containing an estimate of the predictive
    variance ad the \code{XX} locations }
  \item{q1 }{ a vector containing an estimate of the 5\% quantile of
    the predictive distribution at the \code{XX} locations,
    unless \code{quants = FALSE}}
  \item{q2 }{ a vector containing an estimate of the 95\% quantile of
    the predictive distribution at the \code{XX} locations, unless
    \code{quants = FALSE}}
  \item{yypred }{ if \code{yy != NULL} then this contains the
    predictive probability of the true \code{yy} values at the
    \code{XX} locations}
 \item{ei }{ a vector containing an estimate of the EI statistic,
    unless \code{ei = FALSE}};

  For classification:

  \item{p }{ a \code{nrow(XX)}-by-\code{max(object$y)} \code{matrix} of
    mean class probabilities for each of \code{max(object$y)} classes
    at the predictive data locations }
  \item{entropy }{ a \code{nrow(XX)} vector of predictive entropys
    at the predictive data locations}
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

\seealso{
\code{\link{dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{plot.dynaTree}}, \code{\link{alc.dynaTree}}, 
  \code{\link{entropyX.dynaTree}}}
\examples{
## see the example(s) section(s) of dynaTree and
## update.dynaTree and the demos (demo(package=dynaTree))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
