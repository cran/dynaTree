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
dynaTree(X, y, N = 1000, model = c("constant", "linear", "class", "prior"),
         nu0s20 = c(0,0), ab = c(0.95, 2), minp = NULL, sb = NULL, 
	 nstart = minp, icept = c("implicit", "augmented", "none"), 
         rprop = c("luvar", "luall", "reject"), verb = round(length(y)/10))
dynaTrees(X, y, N = 1000,  R = 10, sub = length(y),
          model = c("constant", "linear", "class", "prior"), nu0s20 = c(0,0),
          ab=c(0.95, 2), minp = NULL, sb = NULL, nstart = minp,
          icept =  c("implicit", "augmented", "none"),
          rprop = c("luvar", "luall", "reject"), XX = NULL, yy = NULL,
	  varstats = FALSE, lhs = NULL, plotit = FALSE, proj = 1,
          rorder = TRUE, verb = round(sub/10),  pverb=round(N/10), ...)
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
 \item{sub}{
   Optional argument allowing only a subset of the \code{length(y)}
   \code{X}-\code{y} pairs to be used in each repeat of 
  \code{dynaTrees}; each repeat will use a different random subset
  of size \code{sub}
 }
  \item{model}{
    indicates the type of model to be used at the leaves of the tree;
    \code{"constant"} and \code{"linear"} apply to regression,
    and \code{"class"} to multinomial classification; finally \code{"prior"}
    was recently added to explore sampled without data
}
 \item{nu0s20}{
   a two-vector indicating Inverse Gamma prior parameters
   \code{c(nu0, sigma20)} for the variance in each leaf node,
   \eqn{\sigma^2}{sigma2}.  A \code{c(0,0)} setting indicates 
   a default, scale-invariant, prior; does not apply to the
   \code{"class"} model 
}
\item{ab}{
   tree prior parameter \code{c(alpha, beta)}; see details below
}
\item{minp}{
    a positive scalar integer describing the smallest allowable
    region in the treed partition; if \code{NULL} (default) a
    suitable minimum is calculated based on \code{dim(X)} and
    the type of \code{model} being fit
  }
  \item{sb}{
    an optional two-vector of positive integers indicating
    \code{c(splitmin, basemax)} for the \code{"linear"} model.
    It gives the first column of
    \code{X} on which treed partitioning is allowed, and the last
    column of \code{X} to use as covariates in the linear model
    at the leaves, respectively
  }
  \item{nstart}{
    a positive scalar integer \code{>= minp} indicating
    the time index at which treed partitioning is allowed to start
  }
  \item{icept}{
    indicates the type of intertcept term used (only applies to
    \code{model="linear"}).  The default, \code{"implicit"} 
    causes the inputs \code{X} to be centered so the intercept can be
    implied as an afterthought; \code{"augmented"} causes the inputs 
    \code{X} to automatically gain a leading column of ones in a way
    that is transparent to the user; and \code{"none"} assumes that
    no intercept is being used, or that the user has pre-treated
    \code{X} to have a column of ones.  The main advantage of
    \code{"implicit"} over \code{"augmented"} is that the former can
    default to a constant model fit if leaf design matrices become
    rank deficient.  The latter defaults to the zero-model in such
    cases
  }
  \item{XX}{
    a design \code{matrix} of predictive locations (where \code{ncol(XX) ==
      ncol(X)}) for \code{dynaTrees}; also see
    \code{\link{predict.dynaTree}} and some explanation in the
    details below
  } 	 
  \item{yy}{
  an optional vector of \dQuote{true} responses at the \code{XX}
  predictive locations at which the log posterior probability are
  to be reported
  }  	 
  \item{varstats}{
  if \code{TRUE} causes the \code{\link{varpropuse}},
  \code{\link{varproptotal}}, and \code{\link{relevance.dynaTree}}
  functions to be called on after each
  repetition to collect the usage proportions of each input variable
 (column of \code{X}); see those documentation files for more details
}
\item{lhs}{
  an optional \code{lhs} argument to
    \code{\link{sens.dynaTree}} if a sensitivity analysis step is
    desired after each restart (\code{XX="sens"})
  }
  \item{plotit}{
    a scalar \code{logical} indicating if the fit should be plotted
    after each of the \code{R} repeats; only applies to 1-d data
    and \code{dynaTrees}
  }
  \item{proj}{
    when \code{ncol(x$X) > 1} and \code{plotit = TRUE}
    this argument is passed to \code{\link{plot.dynaTree}} to make
    a 1-d projection using \code{x$X[,proj]}
  }
  \item{rorder}{
    a scalar \code{logical} indicating if
    the rows of \code{X} (and corresponding components of
    \code{y}) should be randomly re-ordered for repeats
    \code{2:R} in order to assess the how the time-ordering
    of the SMC effects the Monte Carlo
    error; only applies to \code{dynaTrees}.  Alternatively,
    one can specify an \code{nrow(X)}-by-\code{(R-1)} matrix
    of orderings (permutations of \code{1:nrow(X)})
}
\item{rprop}{
    indicates the scheme used to construct a grow proposal.
    The best setting, \code{"luall"} uses the lower (L) and upper (U)
    rectangle method based on \code{minp} (above) as described 
    in the seminal reference in which
    the growing location and dimension is sampled uniformly.  It can
    be computationally intensive for large input spaces.  A thriftier
    option (the default) in this case is \code{"luvar"} which uniformly chooses the
    splitting variable first and then uses the LU method marginally.
    Thriftier still is \code{"reject"} which just proposes uniformly
    in the bounding leaf rectangle and rejects subsequent grows that
    lead to partitions with too few data points; (see the \code{minp}
    argument)
}
\item{verb}{
   a positive scalar integer indicating how many time steps
   (iterations) should pass before a progress statement is
   printed to the console; a value of \code{verb = 0} is quiet
 }
\item{pverb}{
   a positive scalar integer indicating after many particles
   should be processed for prediction before a progress statement is
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
  time-orderings of the data indices and also
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
  
  The tree prior is specified by \code{ab=c(alpha, beta)}
  via the and \code{minp}.
  It was originally described by Chipman et al., (1998, 2002)
  \deqn{p_{\mbox{\tiny split}}(\eta, \mathcal{T}) =
       \alpha*(1+\eta)^\beta}{p(split leaf eta) = alpha*(1+depth(eta))^(-beta)}
     and subsequently augmented to enforce a minimum number of points
     (\code{minp}) in each region.

  Once a \code{"dynaTree"}-class object has been built (by
  \code{dynaTree}), predictions and estimates of sequential design and
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
  \item{N }{ the number of particles used }
  \item{X }{ a copy of the design matrix \code{X} }
  \item{y }{ a copy of the responses \code{y} }
  \item{model}{ a copy of the specified leaf model }
  \item{params}{ a vector containing \code{c(nu0s20, alpha, beta, minp,
      sb, icept, rprop)}, where the latter two are in integer form }
  \item{verb }{ a copy of the verbosity argument }
  \item{lpred }{ a vector of \code{log} posterior probabilities
    for each observation, conditional on the ones previous,
    for all time \code{(2*minp):T}; see \code{getBF} for calculating
  Bayes factors from these}
  \item{icept}{ a copy of the intercept argument }
  \item{time }{ the total computing time used to build the
    particle cloud}
  \item{num }{ a \dQuote{pointer} to the \code{C}-side particle
    cloud; see the note below }

  -\cr 
  The \code{dynaTrees} function can obtain predictive samples
  (via \code{\link{predict.dynaTree}}) at each of the \code{R}
  repeats.  Therefore, the \code{"dynaTree"} object returned contains
  extra fields collecting these predictive samples, primarily
  comprising of \code{R} columns of information for each of the fields
  returned by \code{\link{predict.dynaTree}}; see that function for
  more details.  Likewise, when \code{varstats = TRUE} the returned
  object also contains \code{vpu}, \code{vpt} and \code{parde[} fields
  whose columns contain the \code{\link{varpropuse}} and
  \code{\link{varproptotal}} outputs.

  Likewise, \code{dynaTrees}, can provide variable usage summaries
  if \code{varstats = TRUE}, in which case the output includes
  \code{vpu} and \code{vpt} fields; See \code{\link{varpropuse}}
  and \code{\link{varproptotal}} for more details

  The \code{dynaTrees} function does not return \code{num} since 
  it does not leave  any allocated particle clouds on the \code{C}-side
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

  \url{http://bobby.gramacy.com/r_packages/dynaTree/}

}
\author{
  Robert B. Gramacy \email{rbg@vt.edu}, \cr
  Matt Taddy \email{taddy@chicagobooth.edu}, and \cr
  Christoforos Anagnostopoulos \email{christoforos.anagnostopoulos06@imperial.ac.uk}
}

\note{
  As mentioned in the details section, above, the
  \code{dynaTree} function returns a pointer to a particle
  cloud allocated in \code{C}.  This pointer is used
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

  To copy a \code{"dynaTree"}-class object, use
  \code{\link{copy.dynaTree}}, which will also copy the \code{C}-side
  memory allocated to the object
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{predict.dynaTree}}, \code{\link{update.dynaTree}},
  \code{\link{plot.dynaTree}}, \code{\link{deletecloud}},
  \code{\link{copy.dynaTree}}, \code{\link{getBF}}, 
  \code{\link{varpropuse}}, \code{\link{varproptotal}},
  \code{\link{sens.dynaTree}}, \code{\link{relevance.dynaTree}}
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

## try duplicating the object
parab.fit.copy <- copy(parab.fit)

## must delete the cloud or memory may leak
deletecloud(parab.fit); parab.fit$num <- NULL
## to delete all clouds, do:
deleteclouds()

## for more examples of dynaTree see update.dynaTree

## Motorcycle accident data
library(MASS)
data(mcycle)
Xm <- mcycle[,1]
Ym <- mcycle[,2]
XXm <- seq(min(mcycle[,1]), max(mcycle[,1]), length=100)

R <- 2 ## use R >= 10 for better results
## small R is for faster CRAN checks
## fit constant model with R=2 repeats and predictions
moto.fit <- dynaTrees(Xm, Ym, XX=XXm, R=R, plotit=TRUE)

## plot the averages
plot(moto.fit, ptype="mean")

## clouds automatically deleted by dynaTrees

\dontrun{
## 2-d/3-class classification data
library(plgp)
library(tgp)
xx <- seq(-2, 2, length=20)
XX <- expand.grid(xx, xx)
X <- dopt.gp(125, Xcand=XX)$XX
C <- exp2d.C(X)

## fit a classification model with R=10 repeats, 
class.fit <- dynaTrees(X, C, XX=XX, model="class")

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ tree }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ classif }
