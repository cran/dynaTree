\name{dynaTrees}
\alias{dynaTrees}
\alias{dynaTree}
\title{
Fitting Dynamic Tree Models
}
\description{
A function to initialize and fit dynamic tree models
to regression and classification data by the sequential
Monte Carlo (SMC) method of particle learning (PL)
}
\usage{
dynaTree(X, y, N = 1000, model = c("constant", "linear", "class"),
         alpha = 0.95, beta = 2, minp = NULL, verb = 10)
dynaTrees(X, y, N = 1000,  R = 10,
          model = c("constant", "linear", "class"),
          alpha = 0.95, beta = 2, minp = NULL, XX = NULL,
          plotit = FALSE, rorder = TRUE, verb = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    A design \code{matrix} of real-valued predictors
}
  \item{y}{
    A vector of length \code{nrow(X)} containing real-valued
    responses (for regression) or positive
    integer-valued class labels (for classification)
}
\item{N}{
  a positive scalar integer indicating the number
  of particles to be used
}
  \item{R}{
    a scalar integer \code{>= 2} indicating the number of
    \dQuote{repeats} or passes through the data,
    as facilitated by \code{dynaTrees}; see details below
}
  \item{model}{
    indicates the type of model to be used at the leaves of the tree;
    \code{"constant"} and \code{"linear"} apply to regression,
    and \code{"class"} to multinomial classification
}
  \item{alpha}{
    tree prior parameter; see details below
}
  \item{beta}{
    tree prior parameter; see details below
}
   \item{minp}{
    A positive scalar integer describing the smallest allowable
    region in the treed partition; if \code{NULL} (default) a
    suitable minimum is calculated based on \code{dim(X)} and
    the type of \code{model} being fit
}
  \item{XX}{
    a design matrix of predictive locations (where \code{ncol(XX) ==
      ncol(X)} for \code{dynaTrees}; also see
    \code{\link{predict.dynaTree}} and some explanation in the
    details below
}
  \item{plotit}{
    a scalar \code{logical} indicating if the fit should be plotted
    after each of the \code{R} repeats; only applies to 1-d data
    and \code{dynaTrees}
}
  \item{rorder}{
    a scalar \code{logical} indicating if
    the rows of \code{X} (and corresponding components of
    \code{y}) shoould be randomly re-ordered for repeats
    \code{2:R} in order to assess the how the time-ordering
    of the SMC effects the Monte Carlo
    error; only applies to \code{dynaTrees}.  Alternatively,
    one can specify an \code{nrow(X)}-by-\code{(R-1)} matrix
    of orderings (permutations of \code{1:nrow(X)})
}
\item{verb}{
    a positive scalar integer indicating how many time steps
    (iterations) should pass before a progress statement is
    printed to the console; a value of \code{verb = 0} is quiet
  }
  \item{\dots}{
    extra arguments to \code{\link{predict.dynaTree}} passed
    from \code{dynaTrees}
    }
}
\details{
  The \code{dynaTree} function processes the \code{X} and \code{y}
  pairs serially via PL.  It builds up a particle cloud
  which is stored as an object in \code{C}.  A \dQuote{pointer} to that
  object is the primary return value. The \code{dynaTrees} function
  fits several (\code{R}) different dynamic tree models on different
  time-oderings of the data indices and also
  obtains samples from the posterior predictive distribution at
  new \code{XX} locations.  These predictions can be averaged
  over each repeat, or used to assess the Monte Carlo predictive
  error.

  Three different leaf \code{model}s are supported: two for
  regression and one for classification.  If \code{model == "class"}
  then the \code{y} values must contain representatives from
  every class (\code{1:max(y)}).  For details of these models and
  the complete description of their use at the leaves of the dynamic
  trees, see the Taddy, et al., (2009) reference, below.
  
  The tree prior is specified by \code{alpha}, \code{beta} and
  \code{minp}.  It was originally described by Chipman et al., (1998, 2002)
  \deqn{p_{\mbox{\tiny split}}(\eta, \mathcal{T}) =
       \alpha*(1+\eta)^\beta}{p(split leaf eta) = alpha*(1+depth(eta))^(-beta)}
     and subsequently augmented to enforce a minimum number of points
     (\code{minp}) in each region.

  Once a \code{"dynaTree"}-class object has been built (by
  \code{dynaTree}, predictions and estimates of sequential design and
  optimization criteria can be obtained via
  \code{\link{predict.dynaTree}}, a generic prediction method.
  These values can be used to augment the design, and the
  \code{\link{update.dynaTree}} function can be used to quickly
  update the fit with the augmenting data
}
\value{
  Both functions return an object of class \code{"dynaTree"},
  which is a list containing the following fields

  \item{m }{ \code{ncol(X)} } 
  \item{T }{ \code{nrow(X)} }
  \item{X }{ a copy of the design matrix \code{X} }
  \item{y }{ a copy of the responses \code{y} }
  \item{model}{ a copy of the specified leaf model }
  \item{params}{ a 3-vector containing \code{c(alpha, beta, minp)}}
  \item{verb }{ a copy of the verbosity argument }
  \item{lpred }{ a vector of \code{log} posterior probabilities
    for of each observation, conditional on the ones previous,
    for all time \code{(2*minp):T}; see \code{getBF} for calculating
  Bayes factors from these}
  \item{time }{ the total computing time used to build the
    particle cloud}
  \item{num }{ a \dQuote{pointer} to the \code{C}-side particle
    cloud; see the note below }
      
  The \code{dynaTrees} function also does prediction (via
  \code{\link{predict.dynaTree}}) at each of the \code{R} repeats.
  Therefore, the \code{"dynaTree"} object returned contains
  extra fields collecting these predictive samples, primarily
  comprising of \code{R} columns of information for each
  of the fields returned by  \code{\link{predict.dynaTree}};
  see that function for more details.  The \code{dynaTrees}
  function does not return \code{num} since it does not leave
  any allocated particle clouds on the \code{C}-side
}
\references{
   
  Taddy, M.A., Gramacy, R.B., and Polson, N. (2009).
  \dQuote{Dynamic trees for learning and design}
  Tech. Rep. arXiv:0912.1636, University of Cambridge.

  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing}.
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.

  Chipman, H., George, E., \& McCulloch, R. (1998).
  \emph{Bayesian CART model search (with discussion).}
  Journal of the American Statistical Association, \bold{93},
  935--960.

  Chipman, H., George, E., \& McCulloch, R. (2002).
  \emph{Bayesian treed models.}
  Machine Learning, \bold{48}, 303--324.

  \url{http://www.statslab.cam.ac.uk/~bobby/dynaTree.html}

}
\author{
  Robert B. Gramacy \email{bobby@statslab.cam.ac.uk} and\cr
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\note{
  As mentioned in the details section, above, the
  \code{dynaTree} function returns a pointer to a particle
  cloud allocated in \code{C}.  This is pointer is used
  for prediction, via \code{\link{predict.dynaTree}} and for
  later updating/augmentation of data, via
  \code{\link{update.dynaTree}}.
  This information will not be \dQuote{freed} unless
  the user specifically calls \code{\link{deletecloud}(num)}
  or \code{\link{deleteclouds}()}.  Failing to call one
  of these functions (when done with the corresponding
  object(s)) could result in a memory leak;
  see their documentation for more details.

  The \code{C}-side memory cannot be saved in the workspace,
  so they cannot persist across \code{R} sessions
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{predict.dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{plot.dynaTree}}, \code{\link{deletecloud}},
  \code{\link{getBF}}
}
\examples{
## simple parabolic data
n <- 100
Xp <- sort(runif(n,-3,3))
Yp <- Xp + Xp^2 + rnorm(n, 0, .2)

## fit a piece-wise linear model
parab.fit <- dynaTree(Xp, Yp, model="linear")

## obtain predictions at a new set of locations
## and plot
parab.fit <- predict(parab.fit, XX=seq(-3, 3, length=100))
plot(parab.fit)

## must delete the cloud or memory may leak
deletecloud(parab.fit); parab.fit$num <- NULL

## for more examples of dynaTree see update.dynaTree

## Motorcycle accident data
library(MASS)
data(mcycle)
Xm <- mcycle[,1]
Ym <- mcycle[,2]
XXm <- seq(min(mcycle[,1]), max(mcycle[,1]), length=100)

## fit constant model with 10 repeats and predictions
moto.fit <- dynaTrees(Xm, Ym, XX=XXm, plotit=TRUE)

## plot the averages
plot(moto.fit, ptype="mean")

## clouds automatically deleted by dynaTrees

## 2-d/3-class classification data
library(plgp)
library(tgp)
xx <- seq(-2, 2, length=20)
XX <- expand.grid(xx, xx)
X <- dopt.gp(125, Xcand=XX)$XX
C <- exp2d.C(X)

## fit a classification model with 3 repeats
class.fit <- dynaTrees(X, C, XX=XX, R=3, model="class")

## for plot the output (no generic plotting available)
cols <- c(gray(0.85), gray(0.625), gray(0.4))
par(mfrow=c(1,2))
library(akima)

## plot R-averaged predicted class
mclass <- apply(class.fit$p, 1, which.max)
image(interp(XX[,1], XX[,2], mclass), col=cols,
      xlab="x1", ylab="x2", main="repeated class mean")
points(X)
## plot R-averaged entropy
ment <-  apply(class.fit$entropy, 1, mean)
image(interp(XX[,1], XX[,2], ment),
      xlab="x1", ylab="x2", main="repeated entropy mean")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ tree }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ classif }
