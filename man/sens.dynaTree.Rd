\name{sens.dynaTree}
\alias{sens.dynaTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Monte Carlo Sensitivity Analysis for dynaTree Models
}
\description{
  A Monte Carlo sensitivity analysis using random Latin
  hypercube samples (LHSs) or bootstrap resamples for each particle
  to estimate main effects as well as 1st order and total
  sensitivity indices
}
\usage{
\method{sens}{dynaTree}(object, class = NULL, nns = 1000, nME = 100,
              span = 0.3, method = c("lhs", "boot"),
              lhs = NULL, categ = NULL, verb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
}
  \item{class}{
    only valid for \code{object$model = "class"}, allows the
    user to specify the subset of class labels in
    \code{unique(object$y)} for which sensitivity indices are
    calculated.  The implementation loops over the vector
    of labels provided.  The default of \code{NULL} results in
    \code{class = unique(object$y)}
}
\item{nns}{
  A positive integer scalar indicating the size of each
  LHS or bootstrap drawn for use in the Monte Carlo
  integration scheme underlying the sensitivity analysis;
  the total number of locations is \code{nn.lhs*(ncol(X)+2)}
}
\item{nME}{
  A positive integer scalar indicating number of grid points,
  in each input dimension, upon which
  main effects will be estimated
}
\item{span}{
  A positive real-valued scalar giving the smoothing parameter
  for main effects integration: the fraction of \code{nns} points
  that will be included in a moving average window that is used
  to estimate main effects at the \code{nME} locations in each
  input dimension
}
  \item{method}{
    indicates whether LHS or bootstrap should be used
}
  \item{lhs}{
    if \code{method = "lhs"} then this argument should be
    a \code{\link{list}} with entries \code{rect}, \code{shape}
    and \code{mode} describing the marginal distributions of
    the Latin Hypercube; specify \code{NULL}
    for a default specification for for
    \code{method = "boot"}.  The fields should have the following
    format(s):
    \itemize{
      \item \code{rect}: Optional rectangle describing the domain of
      the uncertainty
    distribution with respect to which the sensitivity is to be
    determined.  This defines the domain from which the LH sample
    is to be taken.  The rectangle should be a \code{matrix} or
    \code{data.frame} with \code{ncol(rect) = 2}, and number of rows
    equal to the dimension of the domain.  For 1-d data, a vector of
    length 2 is allowed.  The default is the input
    data range of each column of (\code{object$X}).

    \item \code{shape}: Optional vector of shape parameters for Beta
    marginal distributions having length \code{ncol(object$X)} and
    elements \code{> 1}, i.e., concave Beta distributions. If
    specified, the uncertainty distribution (i.e. the LHS) is
    proportional to a joint pdf formed by independent Beta
    distributions in each dimension of the domain, scaled and shifted
    to have support defined by \code{rect}.  If unspecified, the
    uncertainty distribution is uniform over \code{rect}. The
    specification \code{shape[i]=0} instructs \code{sens} to treat the
    i'th dimension as a binary variable. In this case, \code{mode[i]}
    is the probability parameter for a bernoulli uncertainty
    distribution, and we must also have \code{rect[i,]=c(0,1)}.

    \item \code{mode}: Optional vector of mode values for the Beta
    uncertainty distribution. Vector of length equal to the dimension
    of the domain, with elements within the support defined by
    \code{rect}.  If \code{shape} is specified, but this is not,
    then the scaled Beta distributions will be symmetric.
  }
}
  \item{categ}{
    A vector of logicals of length \code{ncol(object$X)} indicating
    which, if any, dimensions of the input space should be treated
    as categorical; this input is used to help set the default
    \code{lhs$shape} argument if not specified; the default \code{categ}
    argument is \code{NULL} meaning that the categorical inputs
    are derived from \code{object$X} in a sensible way
}
\item{verb}{
    a positive scalar integer indicating how many predictive locations
    (iterations) after which a progress statement should be
    printed to the console; a (default) value of \code{verb = 0} is quiet
}}
\details{
  Saltelli (2002) describes a Latin Hypercube sampling based method for
  estimation of the 'Sobol' sensitivity indices:
  
  1st Order for input \eqn{i}{i},
  \deqn{S(i) = \mbox{Var}(E[f|x_i])/\mbox{Var}(f),}{S(i) = var(E[f|x[i]])/var(f),}
  where \eqn{x_i}{x[i]} is the  \eqn{i}{i}-th input.
 
  Total Effect for input \eqn{i}{i},
  \deqn{T(i) = E[\mbox{Var}(f|x_{-i})]/\mbox{Var}(f),}{T(i) = E[var(f|x[-i])]/var(f),}
  where \eqn{x_{-i}}{x[-i]} is all inputs except for the \eqn{i}{i}-th.
  
  All moments are with respect to the appropriate marginals of the
  uncertainty distribution \eqn{U}{U} -- that is, the probability
  distribution on the inputs with respect to which sensitivity is being
  investigated.
  Under this approach, the integrals involved are approximated through
  averages over properly chosen samples based on two LH samples
  proportional to U. If \code{nns} is the sample size for the
  Monte Carlo estimate, this scheme requires \code{nns*(ncol(X)+2)}
  function evaluations.

  The \code{sens.dynaTree} function implements the method for unknown functions
  \eqn{f}, through prediction via one of the \pkg{tgp} regression
  models conditional on an observed set of \code{X} locations.
  For each particle, treated as sample from the \code{\link{dynaTree}}
  model posterior,
  the \code{nns*(ncol(X)+2)} locations are drawn randomly from the
  LHS scheme and realizations of the sensitivity indices are
  calculated. Thus we obtain a posterior sample of the indices,
  incorporating variability from both the Monte Carlo estimation and
  uncertainty about the function output. Since a subset of the
  predictive locations are actually an LHS proportional to the
  uncertainty distribution, we can also estimate the main effects
  through simple non-parametric regression (a moving average).

  See the Gramacy, Taddy, \& Wild (2011) reference below for more details.
  
  If \code{method = "boot"} is used then simply replace LHS above
  with a bootstrap resample of the \code{object$X} locations.

  As with prediction, the \code{\link{dynaTrees}} function enables
  repeated calls to \code{sens.dynaTree}
}
\value{

  The object returned is of class \code{"dynaTree"}, which includes a
  copy of the list elements from the \code{object} passed in, 
  with the following (sensitivity-analysis specific)
  additions.

  \item{MEgrid }{ An \code{nME}-by-\code{ncol(object$X)} matrix
    containing the main effects predictive grid at which the following
    \code{MEmean}, \code{MEq1}, and \code{MEq2} quantities were obtained }
  \item{MEmean }{ A \code{\link{matrix}} with \code{ncol(object$X)}
    columns and \code{nME} rows containing the mean main effects
    for each input dimension }
  \item{MEq1 }{ same as \code{MEmean} but containing the 5\% quantiles }
  \item{MEq2 }{ same as \code{MEmean} but containing the 95\% quantiles }
  \item{S }{ An \code{object$N}-row and \code{ncol(object$X)}
    \code{\link{matrix}} containing the posterior (samples) of the
    1st Order Sobol sensitivity indices }
  \item{T }{ same as \code{S} but containing the Total Effect indices }

  In the case of \code{object$model = "class"} the entries
  listed above will themselves be lists with an entry for each
  \code{class} specified on input, or all classes as is the
  default
}
\references{
 
  Saltelli, A. (2002) 
  \emph{Making best use of model evaluations to compute sensitivity indices.}
  Computer Physics Communications, 145, 280-297.

  Gramacy, R.B., Taddy, M.A., and S. Wild (2011).
 \dQuote{Variable Selection and Sensitivity Analysis via
   Dynamic Trees with an Application to Computer Code Performance Tuning}
  arXiv:1108.4739
  
  \url{http://faculty.chicagobooth.edu/robert.gramacy/dynaTree.html}
}
\author{
  Robert B. Gramacy \email{rbgramacy@chicagobooth.edu} and\cr
  Matt Taddy \email{taddy@chicagobooth.edu}}

\note{
  The quality of sensitivity analysis is dependent on the size of
  the LHSs used for integral approximation; as with any Monte
  Carlo integration scheme, the sample size (\code{nns}) must
  increase with the dimensionality of the problem.  The total
  sensitivity indices \eqn{T}{T} are forced non-negative,
  and if negative values occur it is necessary to increase
  \code{nnd}. Postprocessing replaces negative values with \code{NA}
}

\seealso{
  \code{\link{dynaTree}}, \code{\link{predict.dynaTree}},
  \code{\link{relevance.dynaTree}}, 
  \code{\link{varpropuse}}, \code{\link{varproptotal}}
}
\examples{
## friedman data
library(tgp)
f <- friedman.1.data(1000)
X <- f[,1:6]
Z <- f$Y

## fit the model and do the sensitivity analysis
out <- dynaTree(X=X, y=Z, N=1000, ab=c(0.01,2))
## also try with model="linear"

## gather relevance statistics
out <- relevance(out)
boxplot(out$relevance)
abline(h=0, col=2, lty=2)
## relevance stats are not as useful when model="linear"
## since it will appear that x4 and x5 not helpful; these
## interact linearly with the response

## full simulation-based sensitivity analysis
out <- sens.dynaTree(out, verb=100)

## plot the main effects
r <- range(rbind(c(out$MEmean, out$MEq1, out$MEq2)))
par(mfrow=c(1,ncol(out$X)), mar=c(5,3,2,1))
plot(out$MEgrid[,1], out$MEmean[,1], type="l", ylim=r, lwd=2,
     ylab="", xlab=colnames(out$MEmean)[1])
lines(out$MEgrid[,1], out$MEq1[,1], lty=2, lwd=2)
lines(out$MEgrid[,1], out$MEq2[,1], lty=2, lwd=2)
if(ncol(out$X) > 1) {
  for(d in 2:ncol(out$X)) {
    plot(out$MEgrid[,d], out$MEmean[,d], col=d, type="l", ylim=r,
         lwd=2, xlab=colnames(out$MEmean)[d], ylab="")
    lines(out$MEgrid[,d], out$MEq1[,d], col=d, lty=2)
    lines(out$MEgrid[,d], out$MEq2[,d], col=d, lty=2)
  }
}

## Sobol indices
par(mfrow=c(1,2), mar=c(5,4,4,2))
boxplot(out$S, main="first order indices", xlab="inputs")
boxplot(out$T, main="total indices", xlab="inputs")
## these look better when model="linear"

## clean up
deletecloud(out)

## for a classification example using the sensitivity hooks
## in the dynaTrees function,  see the class2d demo
## i.e., demo("class2d")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ htest }
