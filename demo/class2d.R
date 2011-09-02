## Classification (3-class) example with data from the
## plgp package; predictive entropy could be used for
## sequential design, i.e., sequential exploration of
## the boundaries between classes of different labels

## load the necessary libraries
library(dynaTree)
library(plgp)
library(tgp)
library(akima)

## close down old graphics windows
graphics.off()

## predictive locations and true classes
## for comparisons
xx <- seq(-2, 2, length=20)
XX <- expand.grid(xx, xx)
CC <- exp2d.C(XX)

## create the design and data in a bounding rectangle
end <- 125
X <- dopt.gp(end, Xcand=XX)$XX
C <- exp2d.C(X)

## now compare to repetition on a fixed design
objR <- dynaTrees(X, C, XX=XX, model="class", yy=CC)

## plotting colors and grid
cols <- c(gray(0.85), gray(0.625), gray(0.4))
par(mfrow=c(1,2))

## plot R-averaged predicted class
CCpR <- apply(objR$p, 1, which.max)
image(interp(XX[,1], XX[,2], CCpR), col=cols,
      xlab="x1", ylab="x2", main="repeated class mean")
missR <- CCpR != CC
points(X); points(XX[missR,], pch=18, col=2)
## plot R-averaged entropy
image(interp(XX[,1], XX[,2], apply(objR$entropy, 1, mean)), 
      xlab="x1", ylab="x2", main="repeated entropy mean")


## predictive locations and true classes
## for comparisons
library(plgp)
xx <- seq(-2, 2, length=20)
XX <- expand.grid(xx, xx)
X <- dopt.gp(125, Xcand=XX)$XX
C <- exp2d.C(X)

## sensitivity analysis can be done with XX="sens",
## and relevance stats can be obtained with varstats=TRUE
objR <- dynaTrees(X, C, XX="sens", model="class", R=3, varstats=TRUE)

## first look at the relevance statistics, which are gathered
## along with varpropuse and varproptotal when varstats=TRUE
par(mfrow=c(1,1))
boxplot(objR$relevance)
abline(h=0, col=2, lty=2)

## plot main effects and Sobol indices for all classes
Cs <- sort(unique(C))
par(mfrow=c(length(Cs), 3))
for(cls in Cs) {
  plot(objR$MEgrid[,1], objR$MEmean[[cls]][,1], type="l",
       main=paste("class", cls, "main effects"), ylab="main effect",
       xlab="x", ylim=c(0,1))
  lines(objR$MEgrid[,1], objR$MEq1[[cls]][,1], type="l", lty=2)
  lines(objR$MEgrid[,1], objR$MEq2[[cls]][,1], type="l", lty=2)
  lines(objR$MEgrid[,2], objR$MEmean[[cls]][,2], type="l", col=2)
  lines(objR$MEgrid[,2], objR$MEq1[[cls]][,2], type="l", col=2, lty=2)
  lines(objR$MEgrid[,2], objR$MEq2[[cls]][,2], type="l", col=2, lty=2)
  boxplot(objR$S[[cls]], main=paste("class", cls, "S indices"), ylim=c(0,1))
  boxplot(objR$T[[cls]], main=paste("class", cls, "T indices"), ylim=c(0,1))
}
