\name{deletecloud}
\alias{deletecloud}
\alias{deleteclouds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delete C-side Clouds in dynaTree Objects
}
\description{
  Function to free the memory used to represent particle clouds
  stored on the \code{C}-side of a \code{"dynaTree"}-class object --
  essential to avoid memory leaks
}
\usage{
deletecloud(obj)
deleteclouds()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    a \code{"dynaTree"}-class object built by \code{\link{dynaTree}}
    with \code{obj$num != NULL}
}
}
\details{
  \code{deletecloud} frees the particular \code{C}-side memory
  associated to a particular object, and \code{deleteclouds} frees
  all \code{C}-side clouds (i.e., all memory) 
  
  An object with \code{obj$num != NULL} has already had its
  \code{C}-side cloud freed
}
\value{
  These functions do not return any value
}
\author{
  Robert B. Gramacy \email{bobby@statslab.cam.ac.uk} and\cr
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
 \code{\link{dynaTree}}, \code{\link{update.dynaTree}}
}
\examples{
deleteclouds()
## see the dynaTree and update.dynaTree help files
## and the demos (demo(package="dynaTree")) for
## more examples
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ internal }
