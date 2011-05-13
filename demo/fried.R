## load the libraries
library(tgp)
library(dynaTree)

## allocate results vector
rmse <- rep(NA, 5)

## generate data for training 
f <- friedman.1.data(200)    
X <- f[,1:5]
Z <- f$Y
## and testing
ff <- friedman.1.data(1000)
XX <- ff[,1:5]

## MCMC treed constant model
fit.bcart <- bcart(X=X, Z=Z, XX=XX, pred.n=FALSE, R=10,
                   BTE=c(2000,12000,2))
rmse[1] <- sqrt(mean((fit.bcart$ZZ.mean - ff$Ytrue)^2))

## dynamic treed constant model
fit.dtc <- dynaTrees(X=X, y=Z, XX=XX, plotit=FALSE)
rmse[2] <- sqrt(mean((apply(fit.dtc$mean, 1, mean) - ff$Ytrue)^2))
  
## MCMC treed linear model
fit.btlm <- btlm(X=X, Z=Z, XX=XX, pred.n=FALSE, R=10,
                 BTE=c(2000,12000,2))
rmse[3] <- sqrt(mean((fit.btlm$ZZ.mean - ff$Ytrue)^2))

## dynamic treed linear model
fit.dtl <- dynaTrees(X=X, y=Z, XX=XX, model="linear", plotit=FALSE)
rmse[4] <- sqrt(mean((apply(fit.dtl$mean, 1, mean) - ff$Ytrue)^2))

## MCMC GP
fit.bgp <- bgp(X=X, Z=Z, XX=XX, pred.n=FALSE, BTE=c(2000,12000,2))
rmse[5] <- sqrt(mean((fit.bgp$ZZ.mean - ff$Ytrue)^2))

## compare rmses
rmse <- matrix(rmse, nrow=1)
colnames(rmse) <- c("tc", "dtc", "tl", "dtl", "gp")
print(rmse)
