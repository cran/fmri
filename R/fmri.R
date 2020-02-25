fmri.smooth <- function (spm, hmax = 4, adaptation = "aws",
    lkern = "Gaussian", skern = "Plateau", weighted = TRUE, ...)
{
    cat("fmri.smooth: entering function\n")
    args <- sys.call()
    args <- c(spm$call,args)
    if (!(tolower(adaptation) %in% c("none", "aws", "fullaws", "segment"))) {
        adaptation <- "aws"
    }
    ladjust <- if ("ladjust" %in% names(list(...)))
        list(...)[["ladjust"]]
    else 1
    delta <- if ("delta" %in% names(list(...)))
        list(...)[["delta"]]
    else 0
    alpha <- if ("alpha" %in% names(list(...)))
        list(...)[["alpha"]]
    else 0.05
    propagation <- if ("propagation" %in% names(list(...)))
        list(...)[["propagation"]]
    else FALSE
    restricted <- if ("restricted" %in% names(list(...)))
        as.logical(list(...)[["restricted"]])
    else FALSE
    if (!("fmrispm" %in% class(spm))) {
        warning("fmri.smooth: data not of class <fmrispm>. Try to proceed but strange things may happen")
    }
    if (!is.null(attr(spm, "smooth"))) {
        warning("fmri.smooth: Parametric Map seems to be smoothed already!")
    }
    variance <- spm$var
    variance[variance == 0] <- 1e+20
    if (is.null(spm$weights)) {
        weights <- c(1, 1, 1)
    }
    else {
        weights <- spm$weights
    }
    if (is.null(spm$bw)) {
        bw <- rep(0, 3)
    }
    else {
        bw <- spm$bw
    }
    cat("fmri.smooth: smoothing the Statistical Parametric Map\n")
#
#   only employ voxel within mask
#
    mask <- spm$mask
    if(is.null(mask)) mask <- array(TRUE,spm$dim)
# reduce data entries to contain only voxel within mask
    if(is.null(spm$maskOnly)|!spm$maskOnly) spm <- condensefMRI(spm)
    cbeta <- spm$cbeta[mask]
    variance <- variance[mask]
    res <- extractData(spm, what="residuals", maskOnly=TRUE)
    if(is.null(res)) {
        return(warning("Please specify keep=''all'' when calling fmri.lm"))
    }
    ttthat <- switch(tolower(adaptation),
                     aws = aws3D(y = cbeta,
                                  sigma2 = variance,
                                  hmax = hmax,
                                  mask = spm$mask,
                                  wghts = weights,
                                  h0 = bw,
                                  lkern = lkern,
                                  skern = skern,
                                  weighted = weighted,
                                  res = res,
                                  ladjust = ladjust,
                                  testprop = propagation),
                     fullaws = aws3Dfull(y = cbeta,
                                          sigma2 = variance,
                                          hmax = hmax,
                                          mask = spm$mask,
                                          wghts = weights,
                                          lkern = lkern,
                                          skern = skern,
                                          weighted = weighted,
                                          res = res,
                                          ladjust = ladjust,
                                          testprop = propagation),
                     none = aws3D(y = cbeta,
                                   sigma2 = variance,
                                   hmax = hmax,
                                   mask = spm$mask,
                                   qlambda = 1,
                                   wghts = weights,
                                   h0 = bw,
                                   lkern = lkern,
                                   skern = skern,
                                   weighted = weighted,
                                   res = res,
                                   ladjust = ladjust),
                     segment = segm3D(y = cbeta,
                                      sigma2 = variance,
                                      hmax = hmax,
                                      mask = spm$mask,
                                      wghts = weights,
                                      df = spm$df,
                                      h0 = bw,
                                      weighted = weighted,
                                      residuals = res,
                                      ladjust = ladjust,
                                      delta = delta,
                                      alpha = alpha,
                                      restricted = restricted))
    scale <- max(abs(range(ttthat$res)))/32767
    ttthat$res <- writeBin(as.integer(ttthat$res/scale), raw(), 2)
    ttthat$resscale <- scale

    cat("\n")
    cat("fmri.smooth: determine local smoothness\n")
    if (is.null(ttthat$scorr)) {
        bw <- get3Dh.gauss(ttthat$vred, weights)
    }
    else {
        bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B",
            lower = c(0.25, 0.25, 0.25), upper = c(20, 20, 20),
            lag = c(5,5,3), data = ttthat$scorr)$par
        bw[bw <= 0.25] <- 0
        dim(bw) <- c(1, 3)
    }
    rxyz <- c(resel(1, bw[, 1]), resel(1, bw[, 2]), resel(1,
        bw[, 3]))
    dim(rxyz) <- c(dim(bw)[1], 3)
    bw0 <- get3Dh.gauss(ttthat$vred0, weights)
    rxyz0 <- c(resel(1, bw0[, 1]), resel(1, bw0[, 2]), resel(1,
        bw0[, 3]))
    dim(rxyz0) <- c(dim(bw0)[1], 3)
    cat("fmri.smooth: exiting function\n")
    cbeta <- variance <- array(0, dim(mask))
    cbeta[mask] <- ttthat$theta
    variance[mask] <- ttthat$var
    if(!is.null(ttthat$segm)){
        segm <- array(0, dim(mask))
        segm[mask] <- ttthat$segm
    } else {
        segm <- NULL
    }
    z <- list(cbeta = cbeta, var = variance, rxyz = rxyz,
              rxyz0 = rxyz0, scorr = spm$scorr, weights = spm$weights,
              bw = bw, hmax = ttthat$hmax, dim = spm$dim, hrf = spm$hrf,
              segm = segm, mask = mask,
              resscale=ttthat$resscale, res=ttthat$residuals,
              maskOnly=TRUE, call = args)
    if (adaptation == "segment") {
      class(z) <- c( "fmrisegment", "fmridata")
      z$alpha <- alpha
      z$delta <- delta
    } else {
      class(z) <- c("fmrispm", "fmridata")
    }
    z$roixa <- spm$roixa
    z$roixe <- spm$roixe
    z$roiya <- spm$roiya
    z$roiye <- spm$roiye
    z$roiza <- spm$roiza
    z$roize <- spm$roize
    z$roit <- spm$roit
    z$header <- spm$header
    z$format <- spm$format
    z$scorr <- ttthat$scorr
    z$call <- args
    attr(z, "file") <- attr(spm, "file")
    attr(z, "white") <- attr(spm, "white")
    attr(z, "design") <- attr(spm, "design")
    attr(z, "residuals") <- attr(spm, "residuals")
    if (!is.null(attr(spm, "smooth"))) {
        attr(z, "smooth") <- paste("Already smoothed before:\n",
            attr(spm, "smooth"), "\nnow with:\n  adaptation  :",
            adaptation, "\n  bandwidth :", signif(hmax,
                3), "\n  lkern     :", lkern, "\n  skern     :",
            skern, "\n")
    }
    else {
        attr(z, "smooth") <- paste("Smoothed with:\n  adaptation  :",
            adaptation, "\n  bandwidth :", signif(hmax,
                3), "\n  lkern     :", lkern, "\n  skern     :",
            skern, "\n")
    }
    invisible(z)
}

fmri.pvalue <- function(spm, mode="basic", na.rm=FALSE, minimum.signal=0, alpha=0.05 ) {
    args <- sys.call()
    args <- c(spm$call,args)
  cat("fmri.pvalue: entering function\n")

  if (!("fmrispm" %in% class(spm)) ) {
    warning("fmri.pvalue: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (!is.null(attr(spm, "smooth"))) {
    if (!is.null(attr(spm, "residuals"))) {
      type <- "t"
      df <- spm$df
      if(is.null(df)) df <- abs(diff(dim(attr(spm, "design"))))
    } else {
      type <- "norm"
      df <- 1000 # this is actually not needed, placeholder
    }
  } else {
    type <- "t"
    df <- spm$df
  }

    stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "voxelwise"){
      pv <- switch(type, "norm" = pnorm(-stat), "t" = pt(-stat,df))
      thresh <- switch(type, "norm" = qnorm(1-alpha), "t" = qt(1-alpha,df))
    } else if (mode == "Bonferroni"){
      pv <- switch(type, "norm" = pnorm(-stat), "t" = pt(-stat,df))
      alpha <- alpha/sum(spm$mask)
      thresh <- switch(type, "norm" = qnorm(1-alpha), "t" = qt(1-alpha,df))
    } else if (mode == "local") {
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
     } else if (mode == "FDR") {
      pv <- 1-switch(type,"norm"=pnorm(stat),"t"=pt(stat,df))
      ind <- fdr(pv[spm$mask],alpha)
      # restrict to voxel in brain mask
      thresh <- min(stat[spm$mask][ind])
      alpha <- max(pv[spm$mask][ind])
      # thats what needed for scale info in plot.fmripvalue
   } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
    }

  cat("fmri.pvalue: thresholding\n")
  mask <- rep(TRUE,length=prod(spm$dim[1:3]))
  mask[stat < thresh] <- FALSE
  mask <- mask&spm$mask
  pv[!mask] <- NA
  dim(pv) <- spm$dim[1:3]

  if (na.rm) {
    pv[spm$var > 9e19] <- NA
  }
  pv[pv<1e-10] <- 1e-10
  # avoid extremely small values
  cat("fmri.pvalue: exiting function\n")

  z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim, hrf = spm$hrf, mask = spm$mask)

  class(z) <- c("fmripvalue")

  z$roixa <- spm$roixa
  z$roixe <- spm$roixe
  z$roiya <- spm$roiya
  z$roiye <- spm$roiye
  z$roiza <- spm$roiza
  z$roize <- spm$roize
  z$roit <- spm$roit
  z$header <- spm$header
  z$format <- spm$format
  z$call <- args
  z$alpha <- alpha
  z$thresh <- thresh
  attr(z, "file") <- attr(spm, "file")
  attr(z, "white") <- attr(spm, "white")
  attr(z, "design") <- attr(spm, "design")
  if (is.null(attr(spm, "smooth"))) {
    attr(z, "smooth") <- "Not smoothed"
  } else {
    attr(z, "smooth") <- attr(spm, "smooth")
  }
  attr(z, "mode") <- paste("Threshold mode:",mode,"\n")

  z
}
