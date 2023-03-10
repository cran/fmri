fmri.stimulus <- function(scans = 1,
                          onsets = c(1),
                          durations = c(1),
                          TR = 2,
                          times = FALSE,
                          sliceorder = NULL,
                          type = c("canonical", "gamma", "boxcar", "user"),
                          par = NULL,
                          scale = 10,
                          hrf = NULL,
                          verbose = FALSE) {

  ## match the type of HRF function to the defaults
  type <- match.arg(type)
  if ((type == "user") && (!inherits(hrf, "function")))
    stop("HRF type is user, but specified hrf is not a function!")

  ## is information for slice timing present ?

  ## re-calculate design spec. from Scans to Time
  if (!times) {
    onsets <- onsets*TR
    durations <- durations*TR
  }
  slicetiming <- !is.null(sliceorder)
  if(slicetiming) {
     nslices <- length(sliceorder)
     scale <- max(scale,nslices)
     # resulution should at least reflect slice times
     slicetimes <- (1:nslices-1)[sliceorder]/TR*scale
  }
  ## consider microtime
  if(scale<1/TR) scale <- 1/TR
  onsets <- onsets * scale
  durations <- durations * scale
  scans <- scans * TR * scale
  TR <- TR/scale
  slicetiming <- !is.null(sliceorder)
  if(slicetiming) {
     nslices <- length(sliceorder)
     slicetimes <- ceiling((1:nslices-1)[sliceorder]/nslices*scale)
  }
  ## normalization constant for the user-defined to make stimuli comparable
  if (type == "user") shrf <- sum(hrf(0:(ceiling(scans)-1)/scale))

  ## basic consistency checks for design spec
  no <- length(onsets)
  if (length(durations) == 1) {
    durations <- rep(durations, no)
  } else if (length(durations) != no)  {
    stop("Length of duration vector does not match the number of onsets!")
  }

  ## create boxcar function
  if(slicetiming){
    stimulus <- matrix(0,ceiling(scans),nslices)
    for (j in 1:nslices) for(i in 1:no)
       stimulus[pmax(1,onsets[i]:(onsets[i]+durations[i]-1)-slicetimes[j]),j] <- 1
  } else {
    stimulus <- rep(0, ceiling(scans))
    for (i in 1:no) stimulus[onsets[i]:(onsets[i]+durations[i]-1)] <- 1
  }
  ## define some HRF
  ## t: time in seconds
  .canonicalHRF <- function(t, par = NULL) {
    ttpr <- par[1] * par[3]
    ttpu <- par[2] * par[4]
    (t/ttpr)^par[1] * exp(-(t-ttpr)/par[3]) - par[5] * (t/ttpu)^par[2] * exp(-(t-ttpu)/par[4])
  }

  ## t: time in seconds
  .gammaHRF <- function(t, par = NULL) {
    th <- 0.242 * par[1]
    1/(th*factorial(3)) * (t/th)^3 * exp(-t/th)
  }

  ## prepare parameters for HRF
  if (type == "canonical") {
    if (is.null(par)) par <- c(6, 12, 0.9, 0.9, 0.35)
    if (!is.numeric(par[1:5]) || any(is.na(par[1:5]))) {
      warning("parameter vector c(", paste(par, collapse=", "),
              ") for canonical HRF is not numeric or has unsufficient length (<5)!\nUsing default parameters!",
              paste(par <- c(6, 12, 0.9, 0.9, 0.35), collapse=", "))
    }
  } else if (type =="gamma") {
    if (is.null(par)) par <- 4
    if (!is.numeric(par[1])) {
      warning("parameter ", par[1],
              " for gamma HRF is not numeric!\nUsing default parameter!", par <- 4)
    }
  }

  ## convolve with chosen HRF
  if (verbose) cat("fmriStimulus: Using", type, "HRF for stimulus creation\n")
  y <- switch(type,
              canonical = .canonicalHRF(0:(20*scale)/scale, par)/2.885802,
              gamma = .gammaHRF(0:(28*scale)/scale, par),
              boxcar = scale,
              user = hrf(0:(ceiling(scans)-1)/scale)/shrf)
  if(slicetiming) {
     for(j in 1:nslices) {
      stimulus[,j] <-
        convolve(stimulus[,j], rev(y), type="open")[1:dim(stimulus)[1]]
      }
      ## final operations to get BOLD at scan time
    ind <- unique((trunc(min(scale*TR,1)*scale):scans)%/%(scale^2*TR))*scale^2*TR
    stimulus <- stimulus[ind,]/(scale^2*TR)
  } else {
     stimulus <- convolve(stimulus, rev(y), type="open")
    ## final operations to get BOLD at scan time
     ind <- unique((trunc(min(scale*TR,1)*scale):scans)%/%(scale^2*TR))*scale^2*TR
     stimulus <- stimulus[ind]/(scale^2*TR)
  }
  ## return mean corrected stimulus function
  if(slicetiming) sweep(stimulus,2,apply(stimulus,2,mean),"-") else stimulus - mean(stimulus)
}

fmri.design <- function(stimulus,
                        order = 2,
                        cef = NULL,
                        verbose = FALSE) {

  ## create matrices and make consistency checks
  if(is.list(stimulus)){
     nstimulus <- length(stimulus)
     dims <- dim(stimulus[[1]])
     if(!is.null(dims)){
        slicetiming <- TRUE
        nslices <- dims[2]
        scans <- dims[1]
        stims <- array(0,c(scans,nstimulus,nslices))
        for(j in 1:nstimulus){
           if(!all(dim(stimulus[[j]])==dims)) stop("Inconsistent dimensions in stimulus list")
           stims[,j,] <- as.matrix(stimulus[[j]])
        }
     }
  } else {
     slicetiming <- FALSE
     stims <- as.matrix(stimulus)
     dims <- dim(stims)
     nstimulus <- dims[2]
     scans <- dims[1]
     nslices <- 1
     dim(stims) <- c(dims,nslices)
  }
  if (!is.null(cef)) {
    cef <- as.matrix(cef)
    if (dim(stims)[1] != dim(cef)[1])
      stop("Length of stimuli ", dim(stimulus)[1], " does not match the length of confounding effects ", dim(cef)[1])
    neffects <- dim(cef)[2]
  } else {
    neffects <- 0
  }


  ## create empty design matrix and fill first columns with Stimulus

    dz <- c(scans, nstimulus + neffects + order + 1, nslices)
    z <- array(0, dz)
    if (verbose) cat("fmriDesign: Adding stimuli to design matrix\n")
    z[, 1:nstimulus,] <- stims

  ## this is the mean effect
  if (verbose) cat("fmriDesign: Adding mean effect to design matrix\n")
  z[, neffects + nstimulus + 1,] <- 1
  ## now confounding effects
  if (neffects != 0) {
    if (verbose) cat("fmriDesign: Adding", neffects, "confounding effect(s) to design matrix\n")
    z[, (nstimulus + 1):(nstimulus + neffects),] <- cef
}

  ## now the polynomial trend (make orthogonal to stimuli)
  if (order != 0) {
    if (verbose) cat("fmriDesign: Adding polynomial trend of order", order, "to design matrix\n")
    for(k in 1:nslices) {
      ortho <- t(stims[,,k]) %*% stims[,,k]
       hz <- numeric(nstimulus)
       for (i in (neffects + nstimulus + 2):(neffects + nstimulus + order + 1)) {
          z[, i, k] <- (1:scans)^(i - nstimulus - neffects - 1)
       z[, i, k] <- z[, i, k]/mean(z[,i,k])
  }
  }
}
if(dim(z)[3]==1) dim(z) <- dim(z)[1:2]
  ## thats it!
  z
}


fmri.lm <- function(ds,
         z,
         mask = NULL,
         actype = c("smooth", "noac", "ac", "accalc"),
         contrast = c(1),
         verbose = FALSE) {
  ##
  ##  should be faster and more memory efficient
  ##  if length(dim(z)) == 3 the third dimension is assumed to correspond to number of slices
  ##  and code separade design matrices for each slice (slice timing)
  ## get function call as character
  call <- as.character(as.expression(sys.call(-1)))

  if (verbose) cat("fmri.lm: entering function at", format(Sys.time()), "\n")

  ## some settings
  actype <- match.arg(actype)

  ## extract the compressed data
  if(!is.null(mask)) ds <- condensefMRI(ds, mask=mask)
  ttt <- extractData(ds, maskOnly=TRUE)
  ## test dimensionality of object and design matrix
  dy <- ds$dim
  if (length(dy) != 4)
    stop("fmri.lm: Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!")
  dz <- dim(z)
  ldz <- length(dz)
  if(!(ldz%in%c(2,3)))
    stop("fmri.lm: Something is wrong with the design matrix")
  slicetiming <- ldz==3
  if(slicetiming){
     if(dz[3] != dy[3]){
       stop("fmri.lm: Slice timing requested with inapproriate number of slices in design")
       }
  }
  mask <- ds$mask
  dm <- dim(mask)
  mask <- array(as.logical(mask),dm)
  nvoxel <- sum(mask)
## mask needs to be logical for index operations
  if (length(dm) != 3)
    stop("fmri.lm: mask is not three-dimensional array!")
  if (any(dy[1:3] != dm))
    stop("fmri.lm: mask dimensionality does not match functional dataset")
  ## first consider contrast vector! NO test whether it is real contrast!!
  length(contrast) <- dim(z)[2]
  contrast[is.na(contrast)] <- 0
  if(slicetiming){
     slices <- (1:dy[3])[apply(mask,3,sum)>0]
# these slices contain active voxel
     nslices <- length(slices)
     xtx <- array(0,c(dz[2],dz[2],nslices))
     lambda1 <- array(0,c(dz[2],nslices))
     zinv <- array(0,c(dz[1:2],nslices))
     Minv <- array(0,c(2,2,nslices))
     for(i in 1:nslices){
       svdresult <- svd(z[,,slices[i]])
       u <- svdresult$u # remember this for the determination of df
       v <- svdresult$v # remember this for the determination of df
       vt <- t(v)       # remember this for the determination of df
       lambda1[,i] <- 1/svdresult$d
       xtx[,,i] <- svdresult$v %*% diag(lambda1[,i]^2) %*% t(svdresult$v)
       R <- diag(1, dy[4]) - svdresult$u %*% t(svdresult$u)
       zinv[,,i] <- svdresult$u %*% diag(lambda1[,i]) %*% t(svdresult$v)
       ## calculate matrix R for bias correction in correlation coefficient
       ## estimate (see Worsley 2005)
       m00 <- dy[4] - dim(z)[2]
       m01 <- sum(R[,-1]*R[,-dy[4]])
       m11 <- sum(R[-1, -dy[4]] * R[-dy[4], -1] + R[-1, -1] * R[-dy[4], -dy[4]])
       Minv[,,i] <- matrix(c(m11, -m01, -m01, m00), 2, 2)/(m00 * m11 - m01 * m01) # inverse
       anzslicemask <- apply(mask,3,sum)[slices]
       inslice <- rep(1:nslices,anzslicemask)

     }
  } else {
  ## (X^T X)^-1
     svdresult <- svd(z)
     u <- svdresult$u # remember this for the determination of df
     v <- svdresult$v # remember this for the determination of df
     vt <- t(v)       # remember this for the determination of df
     lambda1 <- 1/svdresult$d
     xtx <- svdresult$v %*% diag(lambda1^2) %*% t(svdresult$v)
     R <- diag(1, dy[4]) - svdresult$u %*% t(svdresult$u)
     zinv <- svdresult$u %*% diag(lambda1) %*% t(svdresult$v)
     ## calculate matrix R for bias correction in correlation coefficient
     ## estimate (see Worsley 2005)
     m00 <- dy[4] - dim(z)[2]
     m01 <- sum(R[,-1]*R[,-dy[4]])
     m11 <- sum(R[-1, -dy[4]] * R[-dy[4], -1] + R[-1, -1] * R[-dy[4], -dy[4]])
     Minv <- matrix(c(m11, -m01, -m01, m00), 2, 2)/(m00 * m11 - m01 * m01) # inverse
  }
  rm(R)
  ## now we have z = svdresult$u diag(lambda1^(-1)) t(svdresult$v)

  ## define some variables and make object a matrix
  voxelcount <- prod(dy[1:3])
  if(nvoxel==voxelcount) ttt <- ttt[mask,]
  dim(ttt) <- c(nvoxel, dy[4]) ## ttt only contains voxel in mask
  arfactor <- numeric(nvoxel)
  variance <- numeric(nvoxel)
  if(slicetiming){
    beta <- matrix(0,nvoxel,dim(zinv)[2])
     cxtx <- matrix(0,nvoxel)
     residuals <- ttt
     for(i in 1:nslices){
        cxtx[inslice==i] <- t(contrast) %*% xtx[,,i] %*% contrast
        beta[inslice==i,] <- ttt[inslice==i,] %*% zinv[,,i]
        residuals[inslice==i,] <- ttt[inslice==i,] - beta[inslice==i,] %*% t(z[,,slices[i]])
     }
  } else {
     cxtx <- rep.int(t(contrast) %*% xtx %*% contrast, nvoxel)
     beta <- ttt %*% zinv
     residuals <- ttt - beta %*% t(z)
  }
  ## calculated the paramters and residuals for all voxels

  ## actype == "smooth" ... calc AC, smooth AC, calc prewhitened model
  ## actype == "accalc" ... calc AC, calc prewhitened model
  ## actype == "ac"     ... calc AC only
  ## "accalc" is actually a special case of "smooth" (hmax=1)
  if (actype %in% c("smooth", "accalc", "ac")) {
    if (verbose) {
      cat("fmri.lm: calculating AR(1) model\n")
      #      pb <- txtProgressBar(0, nvoxel, style = 3)
    }
    #    for (i in (1:nvoxel)) {
    #     if (verbose) setTxtProgressBar(pb, i)
    #
    ## calculate the coefficients of ACR(1) time series model
    #    a0 <- residuals[i, ] %*% residuals[i, ]
    #    a1 <- residuals[i, -1] %*% residuals[i, -dim(z)[1]]
    #    an <- Minv %*% c(a0, a1)
    #    if (an[1] != 0) arfactor[i] <- an[2]/an[1]
    #  }
    a0 <- as.vector((residuals^2)%*%rep(1,dy[4]))
    a1 <- as.vector((residuals[, -1]*residuals[, -dim(z)[1]])%*%rep(1,dy[4]-1))
    if(slicetiming){
       an <- matrix(0,2,nvoxel)
       for(i in 1:nslices) an[,inslice==i] <- Minv[,,i] %*% rbind(a0[inslice==i], a1[inslice==i])
    } else {
       an <- Minv %*% rbind(a0, a1)
  }
    arfactor <- an[2,]/an[1,]
    arfactor[is.na(arfactor)] <- 0
    arfactor[arfactor >= 1] <- 0.999
    #    if (verbose) close(pb)

    if (actype == "smooth") {
      har <- 3.52
      if (verbose) cat("fmri.lm: smoothing AR(1) parameters with (har):", har, "... ")
      ## now smooth arfactor (if actype is such) with AWS
      ## arfactor only contains voxel within mask
      arfactor <- aws::smooth3D(arfactor, lkern = "Gaussian", h = har, wghts = ds$weights, mask = mask)
      if (verbose) cat("done\n")
    }
    ##
    ##  We now have estimated AR-coeffs for all voxel
    ##
    if (actype %in% c("smooth", "accalc")) {
      ## re- calculated the linear model with prewithened object
      ## NOTE: sort arfactors and bin them! this combines calculation in
      ## NOTE: voxels with similar arfactor, see Worsley
      step <- 0.01
      arlist <- seq(min(arfactor) - step/2, max(arfactor) + step/2,  length = diff(range(arfactor)) / step + 2)
      if (verbose) {
        cat("fmri.lm: re-calculating linear model with prewhitened object\n")
        pb <- txtProgressBar(0, length(arlist) - 1, style = 3)
      }
      for (i in 1:(length(arlist)-1)) {
        if (verbose) setTxtProgressBar(pb, i)
        if(slicetiming){
          for(j in 1:nslices){
            indar <- (arfactor > arlist[i]) & (arfactor <= arlist[i+1]) & (inslice==j)
          if (sum(indar) > 0) {
            ## create prewhitening matrix
            rho <- mean(arlist[i:(i+1)])
            rho0 <- 1/sqrt(1 - rho^2)
            a <- diag(c(1, rep(rho0, dy[4] - 1)) , dy[4])
            a[col(a) == row(a) - 1] <- -rho * rho0
            zprime <- a %*% z[,,j]
            ## calc SVD of prewhitened design
            svdresult <- svd(zprime)
            ## xtx * <estimate of variance> of prewhitened noise is variance of parameter estimate
            xtx <- svdresult$v %*% diag(1/svdresult$d^2) %*% t(svdresult$v)
            cxtx[indar] <- t(contrast) %*% xtx %*% contrast
            tttprime <- ttt[indar, ] %*% t(a)
            ## estimate parameter
            beta[indar,] <- tttprime %*% svdresult$u %*% diag(1/svdresult$d) %*% t(svdresult$v)
            ## calculate residuals
            residuals[indar,] <- tttprime - beta[indar,] %*% t(zprime)
          }
          }
        } else {
        indar <- (arfactor > arlist[i]) & (arfactor <= arlist[i+1])
        if (sum(indar) > 0) {
          ## create prewhitening matrix
          rho <- mean(arlist[i:(i+1)])
          rho0 <- 1/sqrt(1 - rho^2)
          a <- diag(c(1, rep(rho0, dy[4] - 1)) , dy[4])
          a[col(a) == row(a) - 1] <- -rho * rho0
          zprime <- a %*% z
          ## calc SVD of prewhitened design
          svdresult <- svd(zprime)
          ## xtx * <estimate of variance> of prewhitened noise is variance of parameter estimate
          xtx <- svdresult$v %*% diag(1/svdresult$d^2) %*% t(svdresult$v)
          cxtx[indar] <- t(contrast) %*% xtx %*% contrast
          tttprime <- ttt[indar, ] %*% t(a)
          ## estimate parameter
          beta[indar,] <- tttprime %*% svdresult$u %*% diag(1/svdresult$d) %*% t(svdresult$v)
          ## calculate residuals
          residuals[indar,] <- tttprime - beta[indar,] %*% t(zprime)
        }
        }
      }
      rm(ttt,tttprime)
      if (verbose) close(pb)
      ##  prewhitened residuals don't have zero mean, therefore sweep mean over time from them
      meanres <- residuals%*%rep(1/dy[4],dy[4])
      residuals <- residuals - as.vector(meanres)
    }
  }

  ##
  ##  now calculate variances in prewhitened model
  ##
  variance <- apply(residuals^2, 1, sum) / (dy[4] - dim(z)[2]) * cxtx
  residuals <- t(residuals)
  variance[variance == 0] <- 1e20

  ## calculate contrast of parameters
  cbeta <- beta %*% contrast
  ##
  ##    we now need the full arrays
  ##

  ## re-arrange residual dimensions,
    if (verbose) cat("fmri.lm: calculating spatial correlation ... ")

  lags <- c(5, 5, 3)
  corr <- aws::residualSpatialCorr(residuals,mask,lags=lags,compact=TRUE)
  #
  #    "compress" the residuals
  #
  scale <- max(abs(range(residuals)))/32767
  residuals <- writeBin(as.integer(residuals/scale), raw(), 2)
  ## residuals (condensed to mask) will be returned with results
  gc()

  ## determined local smoothness
  bw <- optim(c(2, 2, 2),
              corrrisk,
              method = "L-BFGS-B",
              lower = c(.59, .59, .59),
              upper = c(10, 10, 10),
              lag = lags,
              data = corr)$par
  bw[bw < .6] <- 0
  dim(bw) <- c(1, 3)
  if ((max(bw) > 4 ) || (corr[lags[1], 1, 1] + corr[1, lags[2], 1] + corr[1, 1, lags[3]] > 0.5))
    warning(paste("fmri.lm: Local smoothness characterized by large bandwidth ", bw[1], bw[2], bw[3], " check residuals for structure", collapse = ","))
  rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, bw[3]))
  dim(rxyz) <- c(1, 3)

  if (verbose) cat("done\n")
  ## re-arrange dimensions
  beta0 <- beta
  beta <- matrix(0,voxelcount,dim(z)[2])
  beta[mask,] <- beta0
  rm(beta0)
  dim(beta) <- c(dy[1:3], dim(z)[2]) # vx * vy * vz * p
  cbeta0 <- cbeta
  cbeta <- array(0,dy[1:3])
  cbeta[mask] <- cbeta0
  rm(cbeta0)
  variance0 <- variance
  variance  <- array(0,dy[1:3])
  variance[mask] <- variance0
  arfactor0 <- arfactor
  arfactor  <- array(0,dy[1:3])
  arfactor[mask] <- arfactor0

  if (verbose) cat("fmri.lm: determining df ... ")
  white <- switch(actype,
                  "smooth" = 1L,
                  "accalc" = 2L,
                  3L)
  if(slicetiming) lambda1 <- apply(lambda1,1,mean)
  cx <- u %*% diag(lambda1) %*% vt %*% contrast
  tau1 <- sum(cx[-1] * cx[-length(cx)]) / sum(cx * cx)
  df <- switch(white,
               abs(diff(dim(z)[1:2])) / (1 + 2*(1 + 2 * prod(har/bw)^0.667)^(-1.5) * tau1^2),
               abs(diff(dim(z)[1:2])) / (1 + 2* tau1^2),
               abs(diff(dim(z)[1:2])))
  if (verbose) cat(df, "done\n")

  ## create the result and leave the function
  result <- list()
  result$beta <- beta
  result$cbeta <- cbeta
  result$var <- variance
  result$mask <- mask
  result$residuals <- residuals
  result$resscale <- scale
  result$maskOnly <- TRUE
  result$arfactor <- arfactor
  result$rxyz <- rxyz
  result$scorr <- corr
  result$weights <- ds$weights
  result$dim <- ds$dim
  if(slicetiming){
     hrf <- matrix(0,dim(z)[1],nslices)
     for(i in 1:nslices) hrf[,i] <- z[,,i] %*% contrast
     result$hrf <- hrf
  } else result$hrf <- z %*% contrast
  result$bw <- bw
  result$df <- df
  result$call <- args
  result$roixa <- ds$roixa
  result$roixe <- ds$roixe
  result$roiya <- ds$roiya
  result$roiye <- ds$roiye
  result$roiza <- ds$roiza
  result$roize <- ds$roize
  result$roit <- ds$roit
  result$header <- ds$header
  result$format <- ds$format
  class(result) <- c("fmridata","fmrispm")

  attr(result, "file") <- attr(ds, "file")
  attr(result, "design") <- z
  attr(result, "white") <- white
  attr(result, "residuals") <- !is.null(scale)

  if (verbose) cat("fmri.lm: exiting function at", format(Sys.time()), "\n")

  invisible(result)
}

sincfilter <- function(t,x,wr=8){
   .Fortran(C_sincfilter,
            as.double(t),
            as.integer(length(t)),
            as.double(x),
            as.integer(length(x)),
            ft = double(length(t)),
            as.integer(wr))$ft
}

slicetiming <- function(fmridataobj, sliceorder=NULL){
#
#  performs sinc interpolation for slicetiming
#
   dy <- fmridataobj$dim
   mask <- fmridataobj$mask
   if(is.null(mask)) mask <- array(TRUE,dy[1:3])
   nvoxel <- sum(mask)
   data <- extractData(fmridataobj, maskOnly=FALSE)
   data <- aperm(data,c(4,1:3))
   if(is.null(sliceorder)) sliceorder <- 1:dim(data)[4]
   if(length(sliceorder)!=dim(data)[4])
      return(warning("Inproper length of sliceorder"))
   sliceorder <- pmax(1,pmin(dim(data)[4],as.integer(sliceorder)))
   newdata <- .Fortran(C_slicetim,
                       as.double(data),
                       as.integer(dim(data)[1]),
                       as.integer(dim(data)[2]),
                       as.integer(dim(data)[3]),
                       as.integer(dim(data)[4]),
                       slicetimed=double(prod(dim(data))),
                       double(dim(data)[1]),
                       as.integer(sliceorder))$slicetimed
    dim(newdata) <- dim(data)
    newdata <- aperm(newdata,c(2:4,1))
    dim(newdata) <- c(prod(dy[1:3]),dy[4])
    datascale <- max(abs(range(newdata)))/32767
    fmridataobj$ttt <- writeBin(as.integer(newdata[mask,]/datascale),raw(),2)
    fmridataobj$datascale <- datascale
    fmridataobj$maskOnly <- TRUE
    fmridataobj
}
