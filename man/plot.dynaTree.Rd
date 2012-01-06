\name{plot.dynaTree}
\alias{plot.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Predictive Distributions of Dynamic Tree models
}
\description{
Plotting predictive distributions constructed from dynamic
tree (regression) models for 1-d data -- provided primarily
for use in our 1-d examples and for illustrative purposes
}
\usage{
\method{plot}{dynaTree}(x, proj = 1, add = FALSE, ylim = NULL,
              col = 2, lwd = 1, ptype = c("each", "mean"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{add}{
    a scalar \code{logical} indicating if the lines/points should
    be \dQuote{added} to an existing plot
  }
  \item{proj}{
    when \code{ncol(x$X) > 1} this argument can be used to plot
    a 1-d projection by specifying which column of \code{x$X}
    should be used to make the plot
  }
  \item{ylim}{
    user-specified y-axis limits values; see \code{\link{plot}}
}
\item{col}{
  user-specified color value; see \code{\link{plot}}
}
  \item{lwd}{
    user-specified line-width value; see \code{\link{plot}}
}
  \item{ptype}{
    type of plot used to visualize several predictive
    samples obtained from \code{\link{dynaTrees}}: \code{"each"}
    shows each surface with its own set of three lines, and
    \code{"mean"} shows the three lines obtained by averaging
}
  \item{\dots}{
    other arguments to the generic \code{\link{plot}} method
}
}
\details{
  This plotting function only handles the predictive output from
  1-dimensional regression \code{\link{dynaTree}} models as obtained by
  first calling \code{\link{dynaTree}} and then
  \code{\link{predict.dynaTree}} on the resulting output at new
  \code{XX} locations.  It is provided to help make the illustration
  of our 1-d examples easier and to serve as an aid in a user's
  development of custom plotting functions in higher dimensions
}
\value{
The only output of this function is a pretty plot
}
\references{
 Taddy, M.A., Gramacy, R.B., and Polson, N. (2011).
  \dQuote{Dynamic trees for learning and design}
  Journal of the American Statistical Association, 106(493), pp. 109-123;
  arXiv:0912.1586

  \url{http://faculty.chicagobooth.edu/robert.gramacy/dynaTree.html}}
\author{
Robert B. Gramacy \email{rbgramacy@chicagobooth.edu} and\cr
Matt Taddy \email{taddy@chicagobooth.edu}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{predict.dynaTree}}, \code{\link{dynaTree}},
  \code{\link{update.dynaTree}}
}
\examples{
## see dynaTree, dynaTrees and update.dynaTree for examples
## which use this plot function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ hplot }
