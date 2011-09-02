library(dynaTree)

## read in the data
data(elec2)
X <- elec2[,1:4]
y <- drop(elec2[,5])

## predictive likelihood for repated trials
T <- nrow(X)
hits <- rep(NA, T)
hits <- data.frame(f=hits, h=hits, hd=hits, hr=hits, hdr=hits)

## random data
n <- 25

## fit the initial model
N <- 1000
ffit <- dynaTree(X[1:n,], y[1:n], N=N, model="class")
hfit <- copy(ffit)
hdfit <- copy(ffit)
hrfit <- copy(ffit)
hdrfit <- copy(ffit)

w <- 1
for(t in (n+1):T) {

  ## predict the next data point
  ## full model
  ffit <- predict(ffit, XX=X[t,], yy=y[t])
  hits$f[t] <- which.max(ffit$p) == y[t]

  ## historical retiring model
  hfit <- predict(hfit, XX=X[t,], yy=y[t])
  hits$h[t] <- which.max(hfit$p) == y[t]

  ## historical retiring model, with drift
  hdfit <- predict(hdfit, XX=X[t,], yy=y[t])
  hits$hd[t] <- which.max(hdfit$p) == y[t]
  
  ## historical retiring model, with rejuvenation
  hrfit <- predict(hrfit, XX=X[t,], yy=y[t])
  hits$hr[t] <- which.max(hrfit$p) == y[t]

  ## historical retiring model, with drift and rejuvenation
  hdrfit <- predict(hdrfit, XX=X[t,], yy=y[t])
  hits$hdr[t] <- which.max(hdrfit$p) == y[t]
  
  ## sanity check retiring index
  if(any(hfit$X[w,] != X[t-n,])) stop("bad retiring in h")
  if(any(hdfit$X[w,] != X[t-n,])) stop("bad retiring in hd") 
  if(any(hrfit$X[w,] != X[t-n,])) stop("bad retiring in hr")
  if(any(hdrfit$X[w,] != X[t-n,])) stop("bad retiring in hdr") 
  
  ## retire
  hfit <- retire(hfit, w)
  hdfit <- retire(hdfit, w, lambda=0.9)
  hrfit <- retire(hrfit, w)
  hdrfit <- retire(hdrfit, w, lambda=0.9)
  
  ## update retiring index
  w <- w + 1; if(w >= n) w <- 1
  
  ## rejuvenate  every 100
  if(t %% n == 0) {
    stop()
    hrfit <- rejuvenate(hrfit, odr=1:(n-1), verb=0)
    hdrfit <- rejuvenate(hdrfit, odr=1:(n-1), verb=0)
  }
    
  ## update with new point
  ffit <- update(ffit, X[t,], y[t], verb=100)
  hfit <- update(hfit, X[t,], y[t], verb=100)
  hdfit <- update(hdfit, X[t,], y[t], verb=100)
  hrfit <- update(hrfit, X[t,], y[t], verb=100)
  hdrfit <- update(hdrfit, X[t,], y[t], verb=100)
}

## free C-side memory
deleteclouds()

## plotting a moving window of hit rates over time
apply(hits, 2, mean, na.rm=TRUE)
n <- 25
rhits <- matrix(0, nrow=nrow(hits), ncol=ncol(hits))
for(i in (n+1):nrow(hits)) {
  rhits[i,] <- 0.05*as.numeric(hits[i,]) + 0.95*rhits[i-1,]
}
matplot(rhits, type="l")
