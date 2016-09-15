### create expected BOLD convoluted with HRF
fmri.stimulus <- function(scans = 1,
                          onsets = c(1),
                          durations = c(1),
                          TR = 2,
                          times = FALSE,
                          type = c("canonical", "gamma", "boxcar", "user"),
                          par = NULL,
                          scale = 10,
                          hrf = NULL,
                          verbose = FALSE) {
  
  ## match the type of HRF function to the defaults
  type <- match.arg(type)
  if ((type == "user") && (class(hrf) != "function"))
    stop("HRF type is user, but specified hrf is not a function!")
  
  ## re-calculate design spec. from Scans to Time
  if (!times) {
    onsets <- onsets*TR
    durations <- durations*TR
  }
  
  ## consider microtime
  onsets <- onsets * scale
  durations <- durations * scale
  scans <- scans * TR * scale
  TR <- TR/scale
  
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
  stimulus <- rep(0, ceiling(scans))
  for (i in 1:no) stimulus[onsets[i]:(onsets[i]+durations[i]-1)] <- 1
  
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
  stimulus <- convolve(stimulus, rev(y), type="open")
  
  ## final operations to get BOLD at scan time
  stimulus <- stimulus[unique((scale:scans)%/%(scale^2*TR))*scale^2*TR]/(scale^2*TR)  
  
  ## return mean corrected stimulus function
  stimulus - mean(stimulus)
}

## OLD VERSION OF THE FUNCTION!
fmri.stimulus.old <- function(scans=1 ,onsets=c(1) ,durations=c(1),
                              rt=3, times=NULL, mean=TRUE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35) {
  
  mygamma <- function(x, a1, a2, b1, b2, c) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- ( x/d1 )^a1
    c2 <- c * ( x/d2 )^a2
    res <- c1 * exp(-(x-d1)/b1) - c2 * exp(-(x-d2)/b2)
    res
  }
  
  
  
  if (is.null(times)) {
    scale <- 1
  } else {
    scale <- 100
    onsets <- times/rt*scale
    durations <- durations/rt*scale
    rt <- rt/scale
    scans <- scans*scale
  }
  numberofonsets <- length(onsets)
  
  if (length(durations) == 1) {
    durations <- rep(durations,numberofonsets)
  } else if (length(durations) != numberofonsets)  {
    stop("Length of duration vector does not match the number of onsets!")
  }
  stimulus <- rep(0, scans)
  
  for (i in 1:numberofonsets) {
    for (j in onsets[i]:(onsets[i]+durations[i]-1)) {
      stimulus[j] <- 1
    }
  }
  stimulus <- c(rep(0,20*scale),stimulus,rep(0,20*scale))
  #  just fill with zeros to avoid bounding effects in convolve
  hrf <- convolve(stimulus,mygamma(((40*scale)+scans):1, a1, a2, b1/rt, b2/rt, cc))/scale
  hrf <- hrf[-(1:(20*scale))][1:scans]
  hrf <- hrf[unique((scale:scans)%/%scale)*scale]
  
  dim(hrf) <- c(scans/scale,1)
  
  if (mean) {
    hrf - mean(hrf)
  } else {
    hrf
  }
}




### create fmriDesign object (design-matrix etc.) 
fmri.design <- function(stimulus,
                        order = 2,
                        cef = NULL,
                        verbose = FALSE) {
  
  ## create matrices and make consistency checks
  stimulus <- as.matrix(stimulus)
  if (!is.null(cef)) {
    cef <- as.matrix(cef)
    if (dim(stimulus)[1] != dim(cef)[1])
      stop("Length of stimuli ", dim(stimulus)[1], " does not match the length of confounding effects ", dim(cef)[1])   
    neffects <- dim(cef)[2]
  } else {
    neffects <- 0
  }
  
  scans <- dim(stimulus)[1]
  nstimulus <- dim(stimulus)[2]
  
  ## create empty design matrix and fill first columns with Stimulus 
  z <- matrix(0, scans, nstimulus + neffects + order + 1)
  if (verbose) cat("fmriDesign: Adding stimuli to design matrix\n")
  z[, 1:nstimulus] <- stimulus
  
  ## needed to make all effects orthogonal to stimulus
  ortho <- t(stimulus) %*% stimulus
  
  ## this is the mean effect
  if (verbose) cat("fmriDesign: Adding mean effect to design matrix\n")
  z[, neffects + nstimulus + 1] <- 1
  
  ## now confounding effects (make orthogonal to stimuli)
  if (neffects != 0) {
    if (verbose) cat("fmriDesign: Adding", neffects, "confounding effect(s) to design matrix\n")
    z[, (nstimulus + 1):(nstimulus + neffects)] <- cef
    hz <- numeric(nstimulus)
    for (i in (nstimulus + 1):(nstimulus + neffects)) {
      for (j in 1:nstimulus) hz[j] <- z[, j] %*% z[, i]
      z[, i] <- z[, i] + as.vector(as.matrix(stimulus) %*% as.vector(lm(-hz~ortho-1)$coeff))
    }
  }
  
  ## now the polynomial trend (make orthogonal to stimuli)
  if (order != 0) {
    if (verbose) cat("fmriDesign: Adding polynomial trend of order", order, "to design matrix\n")
    hz <- numeric(nstimulus)
    for (i in (neffects + nstimulus + 2):(neffects + nstimulus + order + 1)) {
      z[, i] <- (1:scans)^(i - nstimulus - neffects - 1)
      z[, i] <- z[, i]/mean(z[,i])
      for (j in 1:nstimulus) hz[j] <- z[,j]%*%z[,i]
      z[, i] <- z[, i] + as.vector(as.matrix(stimulus) %*% as.vector(lm(-hz~ortho-1)$coeff))
    }
  }
  
  ## thats it!
  z
}

## OLD version of the function
fmri.design.old <- function(hrf, order=2) {
  stimuli <- dim(hrf)[2]
  scans <- dim(hrf)[1]
  
  z <- matrix(0, scans, stimuli+order+1)
  
  for (i in 1:stimuli) {
    z[,i] <- hrf[,i]
  }
  
  ortho <- matrix(0, stimuli, stimuli)
  for (i in 1:stimuli) {
    for (j in 1:stimuli) {
      ortho[i,j] <- z[,i]%*%z[,j]
    }
  }
  
  z[,stimuli+1] <- 1
  
  if (order != 0) {
    for (i in (stimuli+2):(stimuli+order+1)) {
      z[,i] <- (1:scans)^(i-stimuli-1)
      z[,i] <- z[,i]/mean(z[,i])
      hz <- numeric(stimuli)
      for (j in 1:stimuli) {
        hz[j] <- z[,j]%*%z[,i]
      }
      tmp <- lm(-hz~ortho-1)
      z[,i] <- z[,i] + as.vector(as.matrix(hrf) %*% as.vector(tmp$coeff))
    }
  }
  
  z
}

fmri.lm <- function(ds,
                    z,
                    mask = NULL,
                    actype = c("smooth", "noac", "ac", "accalc"),
                    contrast = c(1),
                    verbose = FALSE) {
  
  ## get function call as character
  call <- as.character(as.expression(sys.call(-1)))
  
  if (verbose) cat("fmri.lm: entering function at", format(Sys.time()), "\n")
  
  ## some settings
  actype <- match.arg(actype)
  
  ## extract the compressed data
  ttt <- extract.data(ds)
  
  ## test dimensionality of object and design matrix
  dy <- dim(ttt)
  if (length(dy) != 4)
    stop("fmri.lm: Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!")
  if (!is.matrix(z))
    stop("fmri.lm: Something is wrong with the design matrix")
  if (is.null(mask)) mask <- ds$mask
  dm <- dim(mask)
  if (length(dm) != 3)
    stop("fmri.lm: mask is not three-dimensional array!")
  if (any(dy[1:3] != dm))
    stop("fmri.lm: mask dimensionality does not match functional dataset")
  
  ## first consider contrast vector! NO test whether it is real contrast!!
  length(contrast) <- dim(z)[2]
  contrast[is.na(contrast)] <- 0
  
  ## (X^T X)^-1
  svdresult <- svd(z)
  u <- svdresult$u # remember this for the determination of df
  v <- svdresult$v # remember this for the determination of df
  vt <- t(v)       # remember this for the determination of df
  lambda1 <- diag(1/svdresult$d)
  xtx <- svdresult$v %*% diag(1/svdresult$d^2) %*% t(svdresult$v)
  ## now we have z = svdresult$u lambda1^(-1) t(svdresult$v)
  
  ## define some variables and make object a matrix
  voxelcount <- prod(dy[1:3])
  dim(ttt) <- c(voxelcount, dy[4])
  arfactor <- numeric(voxelcount)
  variance <- numeric(voxelcount)
  cxtx <- rep.int(t(contrast) %*% xtx %*% contrast, voxelcount)
  
  ## calculate matrix R for bias correction in correlation coefficient
  ## estimate (see Worsley 2005)
  R <- diag(1, dy[4]) - svdresult$u %*% t(svdresult$u)
  m00 <- dy[4] - dim(z)[2]
  m01 <- 0
  for (k in 1:dy[4]) 
    m01 <- m01 + sum(R[k, -1] * R[k, -dy[4]])
  m10 <- m01
  m11 <- 0
  for (k in 1:(dy[4]-1))
    m11 <- m11 + sum(R[k+1, -dy[4]] * R[k, -1] + R[k+1, -1] * R[k, -dy[4]])
  Minv <- matrix(c(m11, -m01, -m10, m00), 2, 2)/(m00 * m11 - m01 * m10) # inverse
  
  ## calculate the paramters and residuals for all voxels
  beta <- ttt %*% svdresult$u %*% lambda1 %*% t(svdresult$v)
  residuals <- ttt - beta %*% t(z)
  
  ## actype == "smooth" ... calc AC, smooth AC, calc prewhitened model
  ## actype == "accalc" ... calc AC, calc prewhitened model
  ## actype == "ac"     ... calc AC only
  ## "accalc" is actually a special case of "smooth" (hmax=1)
  if (actype %in% c("smooth", "accalc", "ac")) { 
    if (verbose) {
      cat("fmri.lm: calculating AR(1) model\n")
      pb <- txtProgressBar(0, sum(mask), style = 3)
    }
    for (i in (1:voxelcount)[mask]) {
      if (verbose) setTxtProgressBar(pb, i)
      
      ## calculate the coefficients of ACR(1) time series model
      a0 <- residuals[i, ] %*% residuals[i, ]
      a1 <- residuals[i, -1] %*% residuals[i, -dim(z)[1]]
      an <- Minv %*% c(a0, a1)
      if (an[1] != 0) arfactor[i] <- an[2]/an[1]
    }
    arfactor[arfactor >= 1] <- 0.999
    if (verbose) close(pb)
    
    if (actype == "smooth") {
      hmax <- 3.52
      if (verbose) cat("fmri.lm: smoothing AR(1) parameters with (hmax):", hmax, "... ")
      dim(arfactor) <- dy[1:3]
      ## now smooth (if actype is such) with AWS
      arfactor <- smooth3D(arfactor, lkern = "Gaussian", hmax = hmax, wghts = ds$weights, mask = mask)
      dim(arfactor) <- voxelcount
      if (verbose) cat("done\n")
    }
    
    if (actype %in% c("smooth", "accalc")) {
      ## re- calculated the linear model with prewithened object
      ## NOTE: sort arfactors and bin them! this combines calculation in
      ## NOTE: voxels with similar arfactor, see Worsley
      step <- 0.01
      arlist <- seq(min(arfactor) - step/2, max(arfactor) + step/2,  length = diff(range(arfactor)) / step + 1)
      if (verbose) {
        cat("fmri.lm: re-calculating linear model with prewithened object\n")
        pb <- txtProgressBar(0, length(arlist) - 1, style = 3)
      }
      for (i in 1:(length(arlist)-1)) {
        if (verbose) setTxtProgressBar(pb, i)
        
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
      if (verbose) close(pb)
      ##  prewhitened residuals don't have zero mean, therefore sweep mean over time from them
      residuals <- .Fortran("sweepm",
                            residuals = as.double(t(residuals)),
                            as.logical(mask),
                            as.integer(dy[1]),
                            as.integer(dy[2]),
                            as.integer(dy[3]),
                            as.integer(dy[4]),
                            PACKAGE = "fmri")$residuals
      dim(residuals) <- c(dy[4], prod(dy[1:3]))
      variance <- apply(residuals^2, 2, sum) / (dy[4] - dim(z)[2]) * cxtx
    } else {
      variance <- apply(residuals^2, 1, sum) / (dy[4] - dim(z)[2]) * cxtx
      residuals <- t(residuals)
    }
  } else { # actype == "noac"
    variance <- apply(residuals^2, 1, sum) / (dy[4] - dim(z)[2]) * cxtx
    residuals <- t(residuals)
  }
  variance[variance == 0] <- 1e20
  
  ## calculate contrast of parameters
  cbeta <- beta %*% contrast
  
  ## re-arrange dimensions
  dim(beta) <- c(dy[1:3], dim(z)[2]) # vx * vy * vz * p
  dim(cbeta) <- dy[1:3]              # vx * vy * vz
  dim(variance) <- dy[1:3]           # vx * vy * vz
  dim(arfactor) <- dy[1:3]           # vx * vy * vz
  dim(residuals) <- dy[c(4, 1:3)]    # n * vx * vy * vz
  
  if (verbose) cat("fmri.lm: calculating spatial correlation ... ")
  
  lags <- c(5, 5, 3)
  corr <- .Fortran("mcorr",
                   as.double(residuals),
                   as.logical(mask),
                   as.integer(dy[1]),
                   as.integer(dy[2]),
                   as.integer(dy[3]),
                   as.integer(dy[4]),
                   scorr = double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE = "fmri")$scorr
  dim(corr) <- lags                     
  
  ## "compress" the residuals
  scale <- max(abs(range(residuals)))/32767
  residuals <- writeBin(as.integer(residuals/scale), raw(), 2)
  
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
  
  if (verbose) cat("fmri.lm: determining df ... ")
  white <- switch(actype,
                  "smooth" = 1L,
                  "accalc" = 2L,
                  3L)
  cx <- u %*% lambda1 %*% vt %*% contrast
  tau1 <- sum(cx[-1] * cx[-length(cx)]) / sum(cx * cx)
  df <- switch(white,
               abs(diff(dim(z))) / (1 + 2*(1 + 2 * prod(hmax/bw)^0.667)^(-1.5) * tau1^2),
               abs(diff(dim(z))) / (1 + 2* tau1^2),
               abs(diff(dim(z))))
  if (verbose) cat(df, "done\n")
  
  ## create the result and leave the function
  result <- list()
  result$beta <- beta 
  result$cbeta <- cbeta 
  result$var <- variance 
  result$mask <- mask 
  result$res <- residuals 
  result$resscale <- scale
  result$arfactor <- arfactor 
  result$rxyz <- rxyz 
  result$scorr <- corr 
  result$weights <- ds$weights 
  result$dim <- ds$dim 
  result$hrf <- z %*% contrast 
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
  result$dim0 <- ds$dim0
  class(result) <- c("fmridata","fmrispm")
  
  attr(result, "file") <- attr(ds, "file")
  attr(result, "design") <- z
  attr(result, "white") <- white
  attr(result, "residuals") <- !is.null(scale)
  
  if (verbose) cat("fmri.lm: exiting function at", format(Sys.time()), "\n")
  
  invisible(result)
}


fmri.lm.old <- function(data,z,actype="smooth",vtype="var",step=0.01,contrast=c(1),vvector=c(1),keep="all") {
  cat("fmri.lm: entering function\n")
  
  hmax <- 3.52
  if (!class(data) == "fmridata") {
    warning("fmri.lm: data not of class <fmridata>. Try to proceed but strange things may happen")
  }
  args <- sys.call()
  
  ttt <- extract.data(data)
  
  if (length(dim(ttt)) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  
  create.arcorrection <- function(scans, rho=0) {
    rho0 <- 1/sqrt(1-rho^2)
    a <- numeric(scans*scans)
    
    a[1] <- 1
    ind <- (2:scans) *(scans+1) - 2*scans
    a[ind] <- -rho*rho0
    a[ind+scans] <- rho0
    dim(a) <- c(scans,scans)
    
    a
  }
  
  # first consider contrast vector! NO test whether it is real contrast!!
  if (length(contrast) <= dim(z)[2]) contrast <- c(contrast,rep(0,dim(z)[2]-length(contrast)))
  length(contrast) <- dim(z)[2]
  if ((length(vvector) < dim(z)[2]) && (length(vvector) > 1)) vvector <- c(vvector,rep(0,dim(z)[2]-length(vvector)))
  
  # first get the SVD for the design matrix
  svdresult <- svd(z)
  u <- svdresult$u
  v <- svdresult$v
  vt <- t(v)
  lambda1 <- diag(1/svdresult$d)
  lambda2 <- diag(1/svdresult$d^2)
  xtx <- v %*% lambda2 %*% vt
  # now we have z = u lambda1^(-1) vt
  
  # define some variables and make ttt a matrix
  dy <- dim(ttt)
  voxelcount <- prod(dy[1:3])
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  arfactor <- rep(0,length=prod(dy[1:3]))
  variance <- rep(0,length=prod(dy[1:3]))
  if (length(vvector) > 1) variancem <- array(0,c(dy[1:3],length(vvector)^2))
  
  # calculate matrix R for bias correction in correlation coefficient
  # estimate
  R <- diag(1,dy[4]) - u %*% t(u)
  m00 <- dy[4] - dim(z)[2]
  m01 <- 0
  for (k in 1:dy[4]) {m01 <- m01 + sum(R[k,-1]*R[k,-dy[4]])}
  m10 <- m01
  m11 <- 0
  for (k in 1:(dy[4]-1)) {m11 <- m11 + sum(R[k+1,-dy[4]]*R[k,-1] + R[k+1,-1]*R[k,-dy[4]]) }
  Minv <- matrix(c(m11,-m01,-m10,m00),2,2)/(m00*m11-m01*m10) # inverse
  
  # calculate the paramters and residuals for all voxels
  beta <- ttt %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  
  # actype == "smooth" ... calc AC, smooth AC, calc prewhitened model
  # actype == "accalc" ... calc AC, calc prewhitened model
  # actype == "ac"     ... calc AC only
  # "accalc" is actually a special case of "smooth" (hmax=1), but we
  # leave it here for clearer function interface, so parameter set is
  # only needed for "smooth"
  if ((actype == "smooth") || (actype == "accalc") || (actype == "ac")) { 
    progress = 0
    cat("fmri.lm: calculating AR(1) model\n")
    for (i in (1:voxelcount)[data$mask]) {
      if (i > progress/100*voxelcount) {
        cat(progress,"% . ",sep="")
        progress = progress + 10
      }
      # calculate the Koeff of ACR(1) time series model
      a0 <- residuals[i,] %*% residuals[i,]
      a1 <- residuals[i,-1] %*% residuals[i,-dim(z)[1]]
      an <- Minv %*% c(a0,a1)
      if (an[1] != 0) {
        arfactor[i] <- an[2]/an[1]
      } else {
        arfactor[i] <- 0
      }
      
      if (arfactor[i] >= 1) arfactor[i] <- 0.999
      
      ### this method does have a bias!
      #      if (sum(abs(residuals[i,]))) {
      #        arfactor[i] <- residuals[i,2:dim(z)[1]] %*% residuals[i,1:(dim(z)[1]-1)] / residuals[i,] %*% residuals[i,]
      #      } else {
      #        arfactor[i] <- 0 # avoid NaN
      #      }
      ### leave it for historical reasons
      
    }
    cat("\n")
    
    if (actype == "smooth") {
      cat("fmri.lm: smoothing with (hmax):",hmax,"\n")
      dim(arfactor) <- dy[1:3]
      # now smooth (if actype is such) with AWS
      hinit <- 1
      #      arfactor <- gkernsm(arfactor,rep(hmax,3)*0.42445)$gkernsm
      arfactor <- smooth3D(arfactor,lkern="Gaussian",hmax=hmax,wghts=data$weights,mask=data$mask)
      dim(arfactor) <- voxelcount
      cat("fmri.lm: finished\n")
    }
    
    if ((actype == "smooth") || (actype == "accalc")) {
      progress = 0
      cat("fmri.lm: re-calculating linear model with prewithened data\n")
      # re- calculated the linear model with prewithened data
      # NOTE: sort arfactors and bin them! this combines calculation in
      # NOTE: voxels with similar arfactor, see Worsley
      variancepart <- rep(0,length=prod(dy[1:3]))
      if (length(vvector) > 1) variancepartm <- array(0,c(sum(as.logical(vvector))^2,prod(dy[1:3])))
      arlist <- seq(range(arfactor)[1]-step/2,range(arfactor)[2]+step/2,length=diff(range(arfactor))/step+1)
      for (i in 1:(length(arlist)-1)) {
        if (i > progress/100*length(arlist)) {
          cat(progress,"% . ",sep="")
          progress = progress + 10
        }
        indar <- as.logical((arfactor > arlist[i]) * (arfactor <=
                                                        arlist[i+1]))
        if(sum(indar)>0){
          a <- create.arcorrection(dy[4],mean(arlist[i:(i+1)])) # create prewhitening matrix
          zprime <- a %*% as.matrix(z)
          svdresult <- svd(zprime) # calc SVD of prewhitened design
          v <- svdresult$v
          vtt <- t(v)
          xtx <- v %*% diag(1/svdresult$d^2) %*% vtt # xtx * <estimate of varince> of prewhitened noise is variance of parameter estimate
          tttprime <- ttt[indar,] %*% t(a)
          beta[indar,] <- tttprime %*% svdresult$u %*% diag(1/svdresult$d) %*% vtt # estimate parameter
          residuals[indar,] <- tttprime - beta[indar,] %*% t(zprime) # calculate residuals
          variancepart[indar] <- t(contrast) %*% xtx %*% contrast # variance estimate
          if (length(vvector) > 1) variancepartm[,indar] <- as.vector(xtx[as.logical(vvector),as.logical(vvector)])
        }
        gc()
      }
      #  prewhitened residuals don't have zero mean, therefore sweep mean over time from them
      residuals <- .Fortran("sweepm",residuals=as.double(t(residuals)),
                            as.logical(data$mask),
                            as.integer(dy[1]),
                            as.integer(dy[2]),
                            as.integer(dy[3]),
                            as.integer(dy[4]),
                            PACKAGE="fmri")$residuals
      dim(residuals) <- c(dy[4],prod(dy[1:3]))
      residuals <- t(residuals)
      b <- rep(1/dy[4],length=dy[4])
      variance <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2])) * variancepart
      if (length(vvector) > 1) variancem <- as.vector((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2])) * t(variancepartm)
      cat("\n")
    } else {
      b <- rep(1/dy[4],length=dy[4])
      cxtx <- t(contrast) %*% xtx %*% contrast
      if (length(vvector) > 1) variancem <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2])) %*% as.vector(xtx[as.logical(vvector),as.logical(vvector)])
      variance <- ((residuals^2 %*% b) * dy[4] / (dy[4]-dim(z)[2])) %*% cxtx
    }
  } else { # actype == "noac"
    # estimate variance, add more paramters if needed
    b <- rep(1/dy[4],length=dy[4])
    cxtx <- t(contrast) %*% xtx %*% contrast
    if (length(vvector) > 1) variancem <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2])) %*% as.vector(xtx[as.logical(vvector),as.logical(vvector)])
    variance <- ((residuals^2 %*% b) * dy[4] / (dy[4]-dim(z)[2])) %*% cxtx
  }
  cbeta <- beta %*% contrast
  # re-arrange dimensions
  dim(beta) <- c(dy[1:3],dim(z)[2])
  dim(cbeta) <- dy[1:3]
  dim(variance) <- dy[1:3]
  if (length(vvector) > 1) dim(variancem) <- c(dy[1:3],sum(as.logical(vvector)),sum(as.logical(vvector)))
  dim(arfactor) <- dy[1:3]
  residuals <- t(residuals)
  
  cat("fmri.lm: calculating spatial correlation\n")
  
  lags <- c(5,5,3)
  corr <- .Fortran("mcorr",as.double(residuals),
                   as.logical(data$mask),
                   as.integer(dy[1]),
                   as.integer(dy[2]),
                   as.integer(dy[3]),
                   as.integer(dy[4]),
                   scorr=double(prod(lags)),
                   as.integer(lags[1]),
                   as.integer(lags[2]),
                   as.integer(lags[3]),
                   PACKAGE="fmri")$scorr
  dim(corr) <- lags                     
  scale <- NULL
  if(keep=="all"){
    qscale <- range(residuals)
    scale <- max(abs(qscale))/32767
    residuals <- writeBin(as.integer(residuals/scale),raw(),2)
  }
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.59,.59,.59),upper=c(10,10,10),lag=lags,data=corr)$par  
  bw[bw<=.6] <- 0
  if( (max(bw) > 2.5 ) || (corr[lags[1],1,1]+corr[1,lags[2],1]+corr[1,1,lags[3]] >0.5) ) warning(paste("Local smoothness characterized by large bandwidth ",bw," check residuals for structure",collapse=","))
  rxyz <- c(resel(1,bw[1]), resel(1,bw[2]), resel(1,bw[3]))
  dim(rxyz) <- c(1,3)
  
  variance[variance == 0] <- 1e20
  if (length(vvector) > 1) variancem[variancem == 0] <- 1e20
  
  if (length(vvector) > 1) {
    cbeta <- beta[,,,1:sum(vvector>0)]
    vwghts <- numeric(sum(vvector>0))
    for (i in 1:sum(vvector>0))  vwghts[i] <- median(variancem[,,,i,i]/variancem[,,,1,1])
  } else {
    vwghts <- 1
  }
  
  cat("fmri.lm: determining df: ")
  if (actype == "smooth") {
    white <- 1
  } else if (actype == "accalc") {
    white <- 2
  } else {
    white <- 3
  }
  cx <- u %*% lambda1 %*% vt %*% contrast
  tau1 <- sum(cx[-1] * cx[-length(cx)]) / sum(cx * cx)
  df <- switch(white,abs(diff(dim(z))) / (1 + 2*(1 + 2 * prod(hmax/bw)^0.667)^(-1.5) * tau1^2) ,
               abs(diff(dim(z))) / (1 + 2* tau1^2)  ,
               abs(diff(dim(z))) )
  cat(df,"\n")
  
  cat("fmri.lm: exiting function\n")
  
  if (keep == "all") {
    result <- list(beta = beta, cbeta = cbeta, var = variance, res =
                     residuals, arfactor = arfactor, rxyz = rxyz, scorr = corr, weights =
                     data$weights, vwghts = vwghts, mask=data$mask, dim =
                     data$dim, hrf = z %*% contrast, resscale=scale, bw=bw, df=df, call=args)
  } else {
    result <- list(cbeta = cbeta, var = variance, rxyz = rxyz, scorr = corr, weights =
                     data$weights, vwghts = vwghts, mask=data$mask, dim = data$dim, 
                   hrf = z %*% contrast, res=NULL, resscale=NULL, bw=bw, df=df, call=args)
  }
  
  if (length(vvector) > 1) {
    result$varm <- variancem
  }
  
  result$roixa <- data$roixa
  result$roixe <- data$roixe
  result$roiya <- data$roiya
  result$roiye <- data$roiye
  result$roiza <- data$roiza
  result$roize <- data$roize
  result$roit <- data$roit
  result$header <- data$header
  result$format <- data$format
  result$dim0 <- data$dim0
  class(result) <- c("fmridata","fmrispm")
  
  attr(result, "file") <- attr(data, "file")
  attr(result, "design") <- z
  attr(result, "white") <- white
  attr(result, "residuals") <- !is.null(scale)
  
  invisible(result)
}



