## Simple 1-d regrssion demos, plotting predictive surfaces
## and the evolution of Bayes factors over time

## load the package
library(dynaTree)

## parabola data
n <- 100
Xp <- runif(n,-3,3)
Xp <- sort(Xp)
Yp <- Xp + Xp^2 + rnorm(n, 0, .2)
XXp <- seq(-3,3,length=100)

## fits and plots
const.p <- dynaTrees(Xp, Yp, XX=XXp, plotit=TRUE)
plot(const.p, ptype="mean", xlab="x", ylab="y",
     main="(mean) parabola with constant model")
linear.p <- dynaTrees(Xp, Yp, XX=XXp, model="linear", plotit=TRUE)
plot(linear.p, ptype="mean", xlab="x", ylab="y",
     main="(mean) parabola with linear model")

## comparison by log Bayes Factor
o <- apply(matrix(runif(n*9), ncol=9), 2, order)
tree.c.p <- dynaTrees(Xp, Yp, rorder=o, verb=0)
tree.l.p <- dynaTrees(Xp, Yp, model="linear", rorder=o, verb=0)
bf.p <- getBF(tree.l.p, tree.c.p)
matplot(bf.p, type="l", lty=1, col="gray", main="parabola",
        xlab="observation", ylab="log Bayes factor")
lines(apply(bf.p,1,mean))
legend("topleft", c("individual runs", "mean"), col=c(8,1), lwd=2, bty="n")

## Motorcycle accident data
library(MASS)
data(mcycle)
Xm <- mcycle[,1]
Ym <- mcycle[,2]
XXm <- seq(min(mcycle[,1]), max(mcycle[,1]), length=100)

## fits and plots
const.m <- dynaTrees(Xm, Ym, XX=XXm, plotit=TRUE)
plot(const.m, ptype="mean", xlab="time", ylab="accel",
     main="(mean) motorcycle with constant model")
linear.m <- dynaTrees(Xm, Ym, XX=XXm, model="linear", plotit=TRUE)
plot(linear.m, ptype="mean", xlab="time", ylab="accel",
     main="(mean) motorcycle with constant model")

## comparison by log Bayes factor
o <- apply(matrix(runif(length(Xm)*9), ncol=9), 2, order)
tree.c.m <- dynaTrees(Xm, Ym, rorder=o, verb=0)
tree.l.m <- dynaTrees(Xm, Ym, model="linear", rorder=o, verb=0)
bf.m <- getBF(tree.l.m, tree.c.m)
matplot(bf.m, type="l", lty=1, col="gray", main="motorcycle",
        xlab="observation", ylab="log Bayes factor")
lines(apply(bf.m,1,mean))
legend("topleft", c("individual runs", "mean"), col=c(8,1), lwd=2, bty="n")
