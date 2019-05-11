fmri.detrend <- function(data,degree=1,nuisance=NULL,accoef=0) {
  if (!class(data) == "fmridata") {
    warning("fmri.lm: data not of class <fmridata>. Try to proceed but strange things may happen")
  }
  cat("Start trend removal \n")
  ttt <- extractData(data)
  dimttt <- dim(ttt)
  mask <- data$mask
  if (length(dimttt) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  n <- dimttt[4]
  z <- rep(1,n)
  if(degree>0) z <- cbind(rep(1,n),poly(1:n,degree))
  if(!is.null(nuisance)){
    if(dim(nuisance)[1]!=n) nuisance <- t(nuisance)
    if(dim(nuisance)[1]!=n) stop("incompatible dimensions of nuisance ts")
    z <- cbind(z,nuisance)
  }
  u <- svd(z,nv=0)$u
  dim(ttt) <- c(prod(dimttt[1:3]),dimttt[4])
  ttt[mask,] <- ttt[mask,] - ttt[mask,]%*%u%*%t(u)
  dim(ttt) <- dimttt
  cat("Finished trend removal \n")
  if(accoef>0){
  cat("Start prewhitening \n")
     rho0 <- 1/sqrt(1-accoef^2)
     rho1 <- accoef*rho0
     dim(ttt) <- c(prod(dimttt[1:3]),dimttt[4])
     ttt[mask,-1] <- rho0*ttt[mask,-1] - rho1*ttt[mask,-n]
     dim(ttt) <- dimttt
  cat("Finished prewhitening  \n")
  }
  data$ttt <- writeBin(as.numeric(ttt),raw(),4)
  invisible(data)
}

smooth.fmridata <- function(data,bw=0,unit=c("SD","FWHM"),what="spatial"){
  if (!class(data) == "fmridata") {
    warning("smooth.fmridata: data not of class <fmridata>. Try to proceed but strange things may happen")
  }
  ttt <- extractData(data)
  dimttt <- dim(ttt)
  if (length(dimttt) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  if(what=="spatial"){
    cat("Start spatial smoothing \n")
  for(i in 1:dimttt[4]){
     ttt[,,,i] <- aws::kernsm(ttt[,,,i],h=bw,unit=unit)@yhat
  }
  }
  if(what=="temporal"){
    cat("Start temporal smoothing \n")
    n <- dimttt[4]
    if(unit=="FWHM") bw=bw2fwhm(bw)
    span <- as.integer(5*bw)
    wghts <- dnorm(0:span,0,bw)
    tttn <- ttt*wghts[1]
    for(i in 1:span){
       tttn[,,,-(1:i)] <- tttn[,,,-(1:i)] + tttn[,,,-(n+1-1:i)]*wghts[i]
       tttn[,,,-(n+1-1:i)] <- tttn[,,,-(n+1-1:i)] + tttn[,,,-(1:i)]*wghts[i]
    }
    ttt <- tttn/(2*sum(wghts)-wghts[1])
}
  cat("Finished smoothing \n")
  data$ttt <- writeBin(as.numeric(ttt),raw(),4)
  invisible(data)
}
