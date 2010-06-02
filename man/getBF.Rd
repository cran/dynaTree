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
  Taddy, M.A., Gramacy, R.B., and Polson, N. (2009).
  \dQuote{Dynamic trees for learning and design}
  Tech. Rep. arXiv:0912.1636, University of Cambridge.
}
\author{
Robert B. Gramacy \email{bobby@statslab.cam.ac.uk} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}
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
o <- apply(matrix(runif(n*9), ncol=9), 2, order)
lpc.p <- dynaTrees(Xp, Yp, rorder=o, verb=0)
lpl.p <- dynaTrees(Xp, Yp, model="linear", rorder=o, verb=0)
bf.p <- getBF(lpl.p, lpc.p)

## plot the log Bayes factors
matplot(bf.p, type="l", lty=1, col="gray", main="parabola",
        xlab="time", ylab="log Bayes factor")

## see demo("reg1d") for further examples
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
