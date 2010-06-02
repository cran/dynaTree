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
objR <- dynaTrees(X, C, XX=XX, model="class")

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
