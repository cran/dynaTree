\name{elec2}
\alias{elec2}
\docType{data}
\title{
 The ELEC2 Data Set
}
\description{
  Electricity Pricing Data Set Exhibiting Concept Drift
}
\usage{data(elec2)}
\format{
  A data frame with 27552 observations on the following 5 variables.
  \describe{
    \item{\code{x1}}{a numeric vector}
    \item{\code{x2}}{a numeric vector}
    \item{\code{x3}}{a numeric vector}
    \item{\code{x4}}{a numeric vector}
    \item{\code{y}}{class label}
  }
}
\details{
  This data has become a benchmark of sorts in streaming
  classification.  It was first described by Harries (1999)
  and used thereafter for several performance comparisons
  [e.g., Baena-Garcia et al. (2006); Kuncheva and Plumpton, (2008)].
  It holds information for the Australian New South Wales (NSW)
  Electricity Market, containing 27552 records dated from
  May 1996 to December 1998, each referring to a period of 30 minutes
  subsampled as the completely observed portion of 45312 total records
  with missing values.
  These records have seven fields: a binary class label,
  two time stamp indicators (day of week, time),
  and four covariates capturing aspects of electricity demand and
  supply.

  An appealing property of this dataset is that it is expected to
  contain drifting data distributions since, during the recording
  period,  the electricity market was expanded to include adjacent
  areas. This allowed for the production surplus of one region
  to be sold in the adjacent region, which in turn dampened price levels.
}
\source{
  M. Harries. \dQuote{Splice-2 Comparative Evaluation: Electricity
    Pricing}. University of New South Wales, School of Computer Science
    and Engineering technical report (1999)
}
\references{

  M. Baena-Garcia, J. del Campo-Avila, R., Fidalgo, A. Bifet,
  R. Gavalda and R. Morales-Bueno. \dQuote{Early drift detection
    method}. \emph{ECML PKDD 2006 Workshop on Knowledge Discovery from
    Data Streams}, pp. 77-86 (2006)

  L.I. Kuncheva C.O. and Plumpton. \dQuote{Adaptive Learning Rate for
    Online Linear Discriminant Classifiers}. \emph{SSPR and SPR 2008,
    Lecture Notes in Computer Science (LNCS)}, 5342, pp. 510-519 (2008)
}
\examples{
## this is a snipet from the "elec2" demo; see that demo
## for a full comparison to dynaTree models which can
## cope with drifting concepts

## set up data
data(elec2)
X <- elec2[,1:4]
y <- drop(elec2[,5])

## predictive likelihood for repated trials
T <- 200 ## use nrow(X) for a longer version,
## short T is for faster CRAN checks
hits <-  rep(NA, T)

## fit the initial model
n <- 25; N <- 1000
fit <- dynaTree(X[1:n,], y[1:n], N=N, model="class")

w <- 1
for(t in (n+1):T) {

  ## predict the next data point
  ## full model
  fit <- predict(fit, XX=X[t,], yy=y[t])
  hits[t] <- which.max(fit$p) == y[t]

  ## sanity check retiring index
  if(any(fit$X[w,] != X[t-n,])) stop("bad retiring")

  ## retire
  fit <- retire(fit, w)
  ## update retiring index
  w <- w + 1; if(w >= n) w <- 1

  ## update with new point
  fit <- update(fit, X[t,], y[t], verb=100)
}

## free C-side memory
deleteclouds()

## plotting a moving window of hit rates over time
rhits <- rep(0, length(hits))
for(i in (n+1):length(hits)) {
  rhits[i] <- 0.05*as.numeric(hits[i]) + 0.95*rhits[i-1]
}

## plot moving window of hit rates
plot(rhits, type="l", main="moving window of hit rates",
     ylab="hit rates", xlab="t")
}
\keyword{datasets}
