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
\method{predict}{dynaTree}(object, XX, quants = TRUE, alc = FALSE,
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
  \item{quants}{
    a scalar \code{logical} indicating if predictive quantiles
    are desired (useful for visualization, but less so for active
    learning)
}
  \item{alc}{
    a scalar \code{logical} indicating if the ALC statistic
    (expected reduction in predictive variance) should be
    calculated and returned
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
}
\value{
  The object returned is of class \code{"dynaTree"}, which includes a
  copy of the list elements \code{object} with the following (predictive)
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
  \item{alc }{ a vector containing an estimate of the ALC statistic,
    unless \code{alc = FALSE}}
  \item{ei }{ a vector containing an estimate of the ALC statistic,
    unless \code{ei = FALSE}}

  For classification:

  \item{p }{ a \code{nrow(XX)}-by-\code{max(object$y)} \code{matrix} of
    mean class probabilities for each of \code{max(object$y)} classes
    at the predictive data locations }
  \item{entropy }{ a \code{nrow(XX)} vector of predictive entropys
    at the predictive data locations}
}
\references{

  Taddy, M.A., Gramacy, R.B., and Polson, N. (2009).
  \dQuote{Dynamic trees for learning and design}
  Tech. Rep. arXiv:0912.1636, University of Cambridge.
  
  \url{http://www.statslab.cam.ac.uk/~bobby/dynaTree.html}
}
\author{
Robert B. Gramacy \email{bobby@statslab.cam.ac.uk} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}}

\seealso{
\code{\link{dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{plot.dynaTree}}}
\examples{
## see the example(s) section(s) of dynaTree and
## update.dynaTree and the demos (demo(package=dynaTree))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
