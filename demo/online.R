## load the libraries
library(tgp)
library(dynaTree)

## designed predictive locations
NN <- 1000
rect <- matrix(rep(c(0,1), 5), ncol=2, byrow=TRUE)
ff <- friedman.1.data(NN)
XX <- ff[,1:5]
YY <- ff$Ytrue

## storage for predictive performance stats
yypred <- data.frame(orig=NA, or=NA, oalc=NA, full=NA)
time <- rmse <- yypred

## random data
n <- 100
f <- friedman.1.data(n)
Xp <- f[,1:5]
Yp <- f$Y

## accumlate
Xall <- Xp
Yall <- Yp
  
## fit a piece-wise linear model; must use "augmented"
## in order to retire data
N <- 1000
rr <- dynaTree(Xp, Yp, N=N, model="linear", icept="augmented")
## comparing four estimators; one using only the initial
## subset of the data; one with random retirement (rr),
## one retirement by ALC; and one on the full data;
## this is one iteration of the experiment from the online
## paper
ralc <- copy(rr) 
  
## predict and record RMSE, predictive probability and
## computational time for the original subset estimator
rr <- predict(rr, XX, yy=YY, quants=FALSE)
yypred$orig <- exp(mean(log(rr$yypred)))
time$orig <- rr$time
rmse$orig <- sqrt(mean((YY-rr$mean)^2))

## now, streaming data aquisition of 1900 more pairs
more <- 1900
for(i in 1:more) {
    
  ## retire by ALC, approx=TRUE not necessary here by just
  ## used for illustrative purposes
  ralc <- alcX(ralc, rect=rect, approx=TRUE)
  rem <- which.min(ralc$alcX)
  ralc <- retire(ralc, rem)    

  ## random removal
  rem <- sample(1:nrow(rr$X), 1)
  rr <- retire(rr, rem)

  ## update with 1 more
  f <- friedman.1.data(1)
  Xp <- f[,1:5]
  Yp <- f$Y

  ## accumulate with new point
  Xall <- rbind(Xall, Xp)
  Yall <- c(Yall, Yp)

  ## update random and ALC estimators with new point
  rr <- update(rr, Xp, Yp, verb=100)
  ralc <- update(ralc, Xp, Yp, verb=100)
}

## final predictive probability calculation for
## random discarding online version
rr <- predict(rr, XX=XX, yy=YY, quants=FALSE)
yypred$or <- exp(mean(log(rr$yypred)))
time$or <- rr$time
rmse$or <- sqrt(mean((YY-rr$mean)^2))

## final predictive probability calculation for
## ALC retiring online version
ralc <- predict(ralc, XX=XX, yy=YY, quants=FALSE)
yypred$oalc <- exp(mean(log(ralc$yypred)))
time$oalc <- ralc$time
rmse$oalc <- sqrt(mean((YY-ralc$mean)^2))
  
## final estimator: fit with all data
full <- dynaTree(Xall, Yall, N=N, model="linear", icept="augmented")

## final predictive probability calculation for full version
full <- predict(full, XX=XX, yy=YY, quants=FALSE)
yypred$full <- exp(mean(log(full$yypred)))
time$full <- full$time
rmse$full <- sqrt(mean((YY-full$mean)^2))

## print the results of the comparison
cat("yypred"); print(yypred)
cat("rmse"); print(rmse)
cat("time"); print(time)

## free C-side memory
