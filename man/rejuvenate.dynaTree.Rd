\name{rejuvenate.dynaTree}
\alias{rejuvenate.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Rejuvenate particles from the dynaTree posterior
}
\description{
  Re-pass the \code{X}-\code{y} pairs in the \code{object}
  in a random (or specified) order to temporarily double the
  size of the particle set
}
\usage{
\method{rejuvenate}{dynaTree}(object, odr = order(runif(length(object$y))),
     verb = round(length(object$y)/10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{odr}{
    an integer vector of \code{length(object$y)} specifying the order
    in which the \code{object$X}-\code{object$y} paris should be
    processed for the rejuvenated particles
}
  \item{verb}{
    a positive scalar integer indicating how many time steps
   (iterations) should pass before a progress statement is
   printed to the console; a value of \code{verb = 0} is quiet
}
}
\details{
  The \code{rejuvenate} function causes the particle set to
  temporarily double, to have size \code{2 * object$N}.  The new
  \code{object$N} particles represent a discrete approximation
  to the \code{\link{dynaTree}} posterior under the ordering specified
  by \code{odr}, which may be random.  Subsequent calls to
  \code{\link{update.dynaTree}} cause the particle set to revert back
  to \code{object$N} particles as only that many are obtained from
  the particle learning resample step.

  This function can be particularly useful in online learning contexts,
  where \code{\link{retire.dynaTree}} is used to retain information
  on discarded data, especially when the data is discarded historically
  to deal with drifting concepts.  Since the new, rejuvenated, particles
  are based only on the active data, \code{object$X}-\code{object$y}
  pairs (and not the retired data via informative leaf priors),
  subsequent \code{\link{update.dynaTree}} steps allow the data
  to dictate if old (informative prior) or new (default prior) particles
  are best for the new concept
}
\value{
   The returned list is the same as \code{\link{dynaTree}} --
   i.e., a \code{"dynaTree"}-class object but with \code{2 * object$N}
   particles.  Note that \code{object$N} is not updated to reflect this
   fact, but the C-side object will indeed have a double particle set.
   Repeated calls to \code{rejuvenate} will cause the particle set to
   double again.
}
\references{
  Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1586

  Anagnostopoulos, C., Gramacy, R.B. (2013) \dQuote{Information-Theoretic 
  Data Discarding for Dynamic Trees on Data Streams.} Entropy, 15(12), 
  5510-5535; arXiv:1201.5568

  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing}.
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.

  \url{https://bobby.gramacy.com/r_packages/dynaTree/}
}
\author{
Robert B. Gramacy \email{rbg@vt.edu}, \cr
Matt Taddy and Christoforos Anagnostopoulos
}
\note{
 The object (\code{object}) must contain a pointer to a particle
  cloud (\code{object$num}) which has not been deleted by
  \code{\link{deletecloud}}.  In particular, it cannot be
  an object returned from \code{\link{dynaTrees}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dynaTree}}, \code{\link{alcX.dynaTree}}, 
  \code{\link{entropyX.dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{retire.dynaTree}}
}
\examples{
## see retire.dynaTree for a combined example
## illustrating rejuvenation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
