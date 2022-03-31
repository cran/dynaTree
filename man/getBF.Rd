\name{getBF}
\alias{getBF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract a Path of (log) Bayes Factors
}
\description{
  Extract a path (log) Bayes factors (BFs) from the log marginal posterior
  probabilities of two \code{"dynaTree"}-class objects
}
\usage{
getBF(obj1, obj2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{obj2}{
    another \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
}
\details{
  Simply calculates a difference in log marginal posterior
  probabilities, setting BFs to zero for initial elements of the
  path where one of the objects has more zero marginal probabilities
  than the other.  The BF is for the model in \code{obj1} over
  \code{obj2}. If the objects are the output of repeated
  fits as obtained from \code{\link{dynaTrees}}, then multiple
  traces are returned
}
\value{
  Returns a vector or \code{matrix} of a trace(s) of Bayes factors that
  can be plotted; see examples below
}
\references{
 Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1586

  Gramacy, R.B., Taddy, M.A., and S. Wild (2011).
   \dQuote{Variable Selection and Sensitivity Analysis via
   Dynamic Trees with an Application to Computer Code Performance Tuning}
  arXiv:1108.4739

  \url{https://bobby.gramacy.com/r_packages/dynaTree/}
}
\author{
Robert B. Gramacy \email{rbg@vt.edu}, \cr
Matt Taddy and Christoforos Anagnostopoulos
}

\seealso{
  \code{\link{dynaTree}}, \code{\link{update.dynaTree}},
  \code{link{logpost}}
}
\examples{
## parabola data
n <- 100
Xp <- sort(runif(n,-3,3))
Yp <- Xp + Xp^2 + rnorm(n, 0, .2)
XXp <- seq(-3,3,length=100)

## comparison by log Bayes Factor
R <- 2 ## use R >= 10 for better results
## small R is for faster CRAN checks
o <- apply(matrix(runif(n*(R-1)), ncol=R-1), 2, order)
lpc.p <- dynaTrees(Xp, Yp, R=R, rorder=o, verb=0)
lpl.p <- dynaTrees(Xp, Yp, model="linear", R=R, rorder=o, verb=0)
bf.p <- getBF(lpl.p, lpc.p)

## plot the log Bayes factors
matplot(bf.p, type="l", lty=1, col="gray", main="parabola",
        xlab="time", ylab="log Bayes factor")

## see demo("reg1d") for further examples
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
