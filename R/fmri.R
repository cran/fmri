fmri.smooth <- function(spm,hmax=4,adaptive=TRUE,adaptation="aws",lkern="Gaussian",skern="Plateau",weighted=TRUE,...) {
  cat("fmri.smooth: entering function\n")
  if(!adaptive) adaptation <- "none"
  if(!(tolower(adaptation)%in%c("none","aws","segment"))) {
      adaptation
  }
  ladjust <- if("ladjust" %in% names(list(...))) list(...)[["ladjust"]] else 1
  fov <- if("fov" %in% names(list(...))) list(...)[["fov"]] else NULL
  thresh <- if("thresh" %in% names(list(...))) list(...)[["thresh"]] else 3.5
  delta <- if("delta" %in% names(list(...))) list(...)[["delta"]] else 0
  propagation <- if("propagation" %in% names(list(...))) list(...)[["propagation"]] else FALSE

  if (!("fmrispm" %in% class(spm))) {
    warning("fmri.smooth: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (!is.null(attr(spm,"smooth"))) {
    warning("fmri.smooth: Parametric Map seems to be smoothed already!")
  }
  
  variance <- spm$var
#  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)
  variance[variance == 0] <- 1e20
  
  if (is.null(spm$weights)) {
    weights <- c(1,1,1)
  } else {
    weights <- spm$weights
  }
  if (is.null(spm$bw)) {
    bw <- rep(0,3)
  } else {
    bw <- spm$bw
  }

  cat("fmri.smooth: smoothing the Statistical Parametric Map\n")
  ttthat <- switch(tolower(adaptation),
                   "aws"=vaws3D(y=spm$cbeta, sigma2=variance, hmax=hmax, mask=spm$mask,
                         wghts=weights, h0=bw, vwghts = spm$vwghts,
                         lkern=lkern,skern=skern,weighted=weighted,res=spm$res,
                         resscale=spm$resscale, ddim=spm$dim,ladjust=ladjust,
                         testprop=propagation),
                   "fullaws"=vaws3Dfull(y=spm$cbeta, sigma2=variance, hmax=hmax,
                         mask=spm$mask,wghts=weights, vwghts = spm$vwghts,
                         lkern=lkern,skern=skern,weighted=weighted,res=spm$res,
                         resscale=spm$resscale, ddim=spm$dim,ladjust=ladjust,
                         testprop=propagation),
                   "none"=vaws3D(y=spm$cbeta, sigma2=variance, hmax=hmax, mask=spm$mask,
                         qlambda = 1, wghts=weights, h0=bw,
                         vwghts = spm$vwghts,lkern=lkern,skern=skern,weighted=weighted,res=spm$res,
                         resscale=spm$resscale, ddim=spm$dim,ladjust=ladjust),
                   "segment"=segm3D(y=spm$cbeta, sigma2=variance, hmax=hmax, mask=spm$mask,
                         wghts=weights, h0=bw,lkern=lkern,weighted=weighted,res=spm$res,
                         resscale=spm$resscale, ddim=spm$dim,ladjust=ladjust,delta=delta,
                         thresh=thresh,fov=fov))
  cat("\n")
  
  cat("fmri.smooth: determine local smoothness\n")
  if(is.null(ttthat$scorr)){
     bw <- get3Dh.gauss(ttthat$vred,weights)
  } else {
     bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.25,.25,.25),upper=c(6,6,6),lag=c(5,5,3),data=ttthat$scorr)$par  
     bw[bw<=.25] <- 0
     dim(bw) <- c(1,3)
  } 
  rxyz <- c(resel(1,bw[,1]), resel(1,bw[,2]), resel(1,bw[,3]))
  dim(rxyz) <- c(dim(bw)[1],3)
  bw0 <- get3Dh.gauss(ttthat$vred0,weights)
  rxyz0 <- c(resel(1,bw0[,1]), resel(1,bw0[,2]), resel(1,bw0[,3]))
  dim(rxyz0) <- c(dim(bw0)[1],3)
  cat("fmri.smooth: exiting function\n")
  if(length(dim(ttthat$theta))==3) dim(ttthat$theta) <- c(dim(ttthat$theta),1)
  if (dim(ttthat$theta)[4] == 1) {
    z <- list(cbeta = ttthat$theta[,,,1], var = ttthat$var, rxyz =
              rxyz, rxyz0 = rxyz0, scorr = spm$scorr, weights =
              spm$weights, vwghts = spm$vwghts, bw=bw, 
              hmax = ttthat$hmax, dim = spm$dim, hrf = spm$hrf, segm = ttthat$segm)
  } else {
    z <- list(cbeta = ttthat$theta, var = ttthat$var, rxyz = rxyz, rxyz0 = rxyz0, 
              scorr = spm$scorr, weights = spm$weights, vwghts = spm$vwghts, bw=bw,
              hmax = ttthat$hmax, dim = spm$dim, hrf = spm$hrf)
  }

  class(z) <- c("fmridata","fmrispm")

  z$roixa <- spm$roixa
  z$roixe <- spm$roixe
  z$roiya <- spm$roiya
  z$roiye <- spm$roiye
  z$roiza <- spm$roiza
  z$roize <- spm$roize
  z$roit <- spm$roit
  z$header <- spm$header
  z$format <- spm$format
  z$dim0 <- spm$dim0
  z$scorr <- ttthat$scorr

  attr(z, "file") <- attr(spm, "file")
  attr(z, "white") <- attr(spm, "white")
  attr(z, "design") <- attr(spm, "design")
  attr(z, "residuals") <- attr(spm, "residuals")

  if (!is.null(attr(spm, "smooth"))) {
    attr(z, "smooth") <-
      paste("Already smoothed before:\n",attr(spm, "smooth"),
            "\nnow with:\n  adaptive  :",as.character(adaptive),
            "\n  bandwidth :",signif(hmax,3),
            "\n  lkern     :",lkern,
            "\n  skern     :",skern,"\n")
  } else {
    attr(z, "smooth") <-
      paste("Smoothed with:\n  adaptive  :",as.character(adaptive),
            "\n  bandwidth :",signif(hmax,3),
            "\n  lkern     :",lkern,
            "\n  skern     :",skern,"\n")      
  }
  z
}

fmri.pvalue <- function(spm, mode="basic", delta=NULL, na.rm=FALSE, minimum.signal=0 ) {
  cat("fmri.pvalue: entering function\n")

  if (!("fmrispm" %in% class(spm)) ) {
    warning("fmri.pvalue: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (!is.null(attr(spm, "smooth"))) {
    if (!is.null(attr(spm, "residuals"))) {
      type <- "t"
      df <- abs(diff(dim(attr(spm, "design"))))
    } else {
      type <- "norm"
      df <- abs(diff(dim(attr(spm, "design")))) # this is actually not needed, placeholder
    }
  } else {
    type <- "t"
    df <- spm$df
  }
  if (df>171) type <- "norm"

  if (length(dim(spm$cbeta)) < 4) {

    stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }        
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
    }

  } else if (!is.null(delta)) {

    l1 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[1]
    l2 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[2]
    theta1 <- atan(l1)
    theta2 <- atan(l2)
    t1 <- spm$cbeta[,,,1]/sqrt(spm$var * spm$vwghts[1])
    t2 <- spm$cbeta[,,,2]/sqrt(spm$var * spm$vwghts[2])
    ratio <- t2/t1
    ratio[t1==0] <- l2 + 1
    w1 <- (t1 + t2 * l1) / sqrt(1+l1^2)
    w2 <- (t1 + t2 * l2) / sqrt(1+l2^2)
    w3 <- (t1 > 0) * (l1 <= ratio) * (ratio <= l2) * sqrt(t1^2 + t2^2)
    stat <- pmax(w1,w2,w3)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm",cone=theta2-theta1)
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm",cone=theta2-theta1)
    }

  } else {

    stat <- spm$cbeta[,,,1]^2/spm$var + spm$cbeta[,,,2]^2/spm$var/spm$vwghts[2]  # Wert der Statistik
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="chisq",df=2)
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="chisq",df=2)
    }

  }
  cat("fmri.pvalue: thresholding\n")
  mask <- rep(TRUE,length=prod(spm$dim[1:3]))
  mask[stat < thresh] <- FALSE
  pv[!mask] <- 1
  dim(pv) <- spm$dim[1:3]

  if (na.rm) {
    pv[spm$var > 9e19] <- 1
  }
  
  cat("fmri.pvalue: exiting function\n")

  z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim, hrf = spm$hrf)
  
  class(z) <- c("fmridata","fmripvalue")

  z$roixa <- spm$roixa
  z$roixe <- spm$roixe
  z$roiya <- spm$roiya
  z$roiye <- spm$roiye
  z$roiza <- spm$roiza
  z$roize <- spm$roize
  z$roit <- spm$roit
  z$header <- spm$header
  z$format <- spm$format
  z$dim0 <- spm$dim0

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




plot.fmridata <- function(x, anatomic = NULL , maxpvalue = 0.05, spm = TRUE,
                            pos = c(-1,-1,-1), type="slice",
                            device="X11", file="plot.png", slice =  1, view = "axial" ,zlim.u =
                            NULL, zlim.o = NULL, col.o = heat.colors(256), col.u = grey(0:255/255), ...) {
  mri.colors <- function (n1, n2, factor=n1/(n1+n2), from=0, to=.2) {
    colors1 <- gray((0:n1)/(n1+n2))
    colors2 <- hsv(h = seq(from,to,length=n2),
                   s = seq(from = n2/(n2+factor*n1) - 1/(2 * (n2+factor*n1)), to =
                     1/(2 * (n2+factor*n1)), length = n2),
                   v = 1,
                   gamma=1)
    list(all=c(colors1,colors2),gray=colors1,col=colors2)
  }

  if ("fmripvalue" %in% class(x)) {

    if ("fmridata" %in% class(anatomic)) {

      img <- show.slice(x, anatomic, maxpvalue = maxpvalue, slice =  slice, view = view, col.u, col.o, zlim.u, zlim.o)
      hex <- c(0:9, LETTERS[1:6])
      hex <- paste(hex[(0:255)%/%16+1],hex[(0:255)%%16+1],sep="")
      color <- paste("#",hex[img[,,1]%/%256+1],hex[img[,,2]%/%256+1],hex[img[,,3]%/%256+1],sep="")

      xxx <- seq(1,dim(img)[1],length=dim(img)[1])
      yyy <- seq(1,dim(img)[2],length=dim(img)[2])
      zzz <- matrix(1:prod(dim(img)[1:2]),nrow = dim(img)[1],ncol = dim(img)[2])
      # display the image
      image(xxx, yyy, zzz, col = color, asp = 1, xlab="",ylab="", ...)
      return(invisible(img))

    } else {

      signal <- x$pvalue
      signal[signal > maxpvalue] <- 1
      signal[signal < 1e-10] <- 1e-10

      signal <- -log(signal)

      if (is.null(anatomic)) anatomic <- array(0,dim=dim(x$pvalue))

      # re-scale anatomic to 0 ... 0.5
      if (diff(range(anatomic)) !=0) {
        anatomic <- 0.5 * (anatomic - range(anatomic,finite=TRUE)[1]) / diff(range(anatomic,finite=TRUE))
      }
      # re-scale signal to 0.5 ... 1
      scale <- range(signal,finite=TRUE)
      if (diff(scale) != 0) {
        signal <- 0.5 + 0.5 * (signal - scale[1]) / diff(scale)
      } else if (scale[1] == 0) {
        signal <- 0.5
      } else {
        signal <- 1
      }
      # create an overlay
      anatomic[signal > 0.5] <- signal[signal > 0.5]
      anatomic[is.na(anatomic)] <- 0
      anatomic[is.infinite(anatomic)] <- 0
    
      if (type == "3d") {
        tt <- fmri.view3d(anatomic,col=mri.colors(255,255)$all,
                          weights=x$weights, scale=scale,scalecol=mri.colors(255,255)$col,
                          type= "pvalue",maxpvalue=maxpvalue,pos=pos)
      } else {
        fmri.view2d(anatomic, device, file, mri.colors(255,255)$all, scale=scale,scalecol=mri.colors(255,255)$col,type="pvalue",maxpvalue=maxpvalue,pos=pos)
      }
    }

  } else if ("fmrispm" %in% class(x)) {

    signal <- if (spm) x$cbeta/sqrt(x$var) else x$cbeta
    
    # re-scale signal to 0 ... 1
    scale <- range(signal,finite=TRUE)
    if (diff(scale) != 0) {
      signal <-  (signal - scale[1]) / diff(scale)
    } else {
      signal <- 0
    }
    signal[is.na(signal)] <- 0
    signal[is.infinite(signal)] <- 0

    ## check !!!!
    quant <- if (!is.null(attr(spm, "smooth"))) qnorm(1-maxpvalue) else qt(1-maxpvalue,length(x$hrf)) 
    
    if (type == "3d") {
      if (spm) {
        tt <- fmri.view3d(signal,col=mri.colors(255,0)$gray,
                          weights=x$weights,
                          scale=scale,scalecol=mri.colors(255,0)$gray,
                          type="spm",pos=pos)
      } else {
        tt <- fmri.view3d(signal,sigma=sqrt(x$var),col=mri.colors(255,0)$gray,
                          weights=x$weights,
                          scale=scale,scalecol=mri.colors(255,0)$gray, type="spm",hrf=x$hrf, quant = quant,pos=pos)
      } 
    } else {
      fmri.view2d(signal, device, file, mri.colors(255,0)$gray, scale=scale,scalecol=mri.colors(255,0)$gray,
                          type="spm",pos=pos)
    }
  } else if ("fmridata" %in% class(x)) {
    signal <- extract.data(x)
    
    # re-scale signal to 0 ... 1
    scale <- range(signal,finite=TRUE)
    if (diff(scale) != 0) {    
      signal <-  (signal - scale[1]) / diff(scale)
    } else {
      signal <- 0
    }
    signal[is.na(signal)] <- 0
    signal[is.infinite(signal)] <- 0
    
    if (type == "3d") {
      tt <- fmri.view3d(signal,col=mri.colors(255,0)$gray,
                        weights=x$weights, scale=scale,scalecol=mri.colors(255,0)$gray, type="data",pos=pos)
    } else {
      fmri.view2d(signal, device, file, mri.colors(255,0)$gray, scale=scale,scalecol=mri.colors(255,0)$gray, type="data",pos=pos)
    }

  } else {
    cat("sorry. plot for this class not implemented\nFalling back to generic function, but this may fail!")
    plot(x)
  }

  if (exists("tt")) invisible(tt)
}

fmri.view2d <- function(ttt, device, file,  col=grey(0:255/255),
                        pos=c(-1,-1,-1), scale=c(0,1), scalecol = col,
                        type = "data", maxpvalue = 0.05) {

  # some basic data properties
#  zlim <- range(ttt)
  zlim <- c(0,1) # requires rescaled data
  dt <- dim(ttt)

  # determine the number of images in x- and y-direction
  partitionx <- partitiony <- ceiling(sqrt(dt[3]+1))
  while ((partitiony-1)*partitionx >= dt[3]+1) partitiony <- partitiony-1

  # choose the output device
  switch (device,
          "png" = png(filename=file, width = 200*partitionx, height = 200*partitiony, pointsize=12, bg="transparent", res=NA),
          "jpeg" = jpeg(filename=file, width = 200*partitionx, height = 200*partitiony,
            quality = 100, pointsize = 12, res=NA),
          "ppm" = bitmap(file,type="ppm",height=2*partitiony,width=2*partitionx,res=64,pointsize=12),
          X11(width=1.2*partitionx,height=1.2*partitiony))
  
  oldpar <- par(mar=c(0.25,0.25,0.25,.25))
  layout(matrix(c(1:dt[3],rep(dt[3]+1,partitionx*partitiony-dt[3])),partitiony,partitionx,byrow=TRUE),height=rep(dt[2],partitiony),width=rep(dt[1],partitionx))
  
  # all the data
  for (i in 1:dt[3]) {
    if (length(dt) == 4) {
      image(1:dt[1],1:dt[2],ttt[,,i,1], xaxt="n", yaxt="n", zlim=zlim, col=col)
    } else {
      image(1:dt[1],1:dt[2],ttt[,,i], xaxt="n", yaxt="n", zlim=zlim, col=col)
    }
    # mark a voxel at pos[1:3]
    if (i == pos[3]) {
      lines(c(0,dt[1])+0.5, c(pos[2],pos[2]), col=2)
      lines(c(pos[1],pos[1]), c(0,dt[2])+0.5, col=2)
    }     
  }

  # add a scale
  par(mgp=c(2,1,0), mar=c(2,0.25,2,0.25))

  if (diff(scale) != 0) {
    if (type == "pvalue") {
      image(seq(-log(maxpvalue),scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
            matrix(rep(seq(-log(maxpvalue),scale[2],length=100),10),100,10),
            yaxt="n",xaxt="n",xlab="", ylab="",zlim=c(-log(maxpvalue),scale[2]), col=scalecol)
      
      lines(c(-log(0.01),-log(0.01)),scale,col=2)
      text(-log(0.01),scale[1]+0.01*diff(scale),pos=4,"1e-2")
      lines(c(-log(0.001),-log(0.001)),scale,col=2)
      text(-log(0.001),scale[1]+0.01*diff(scale),pos=4,"1e-3")
      lines(c(-log(0.0001),-log(0.0001)),scale,col=2)
      text(-log(0.0001),scale[1]+0.01*diff(scale),pos=4,"1e-4")
      lines(c(-log(0.00001),-log(0.00001)),scale,col=2)
      text(-log(0.00001),scale[1]+0.01*diff(scale),pos=4,"1e-5")
      lines(c(-log(0.000001),-log(0.000001)),scale,col=2)
      text(-log(0.000001),scale[1]+0.01*diff(scale),pos=4,"1e-6")
      lines(c(-log(0.0000001),-log(0.0000001)),scale,col=2)
      text(-log(0.0000001),scale[1]+0.01*diff(scale),pos=4,"1e-7")
      lines(c(-log(0.00000001),-log(0.00000001)),scale,col=2)
      text(-log(0.00000001),scale[1]+0.01*diff(scale),pos=4,"1e-8")
      lines(c(-log(0.000000001),-log(0.000000001)),scale,col=2)
      text(-log(0.000000001),scale[1]+0.01*diff(scale),pos=4,"1e-9")
    } else {
      image(seq(scale[1],scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
            matrix(rep(seq(scale[1],scale[2],length=100),10),100,10),
            yaxt="n",xlab="", ylab="",zlim=scale, col=scalecol)
    }    
  }

  # close the device
  par(oldpar)
  switch (device,
          "png" = dev.off(),
          "jpeg" = dev.off(),
          "ppm" = dev.off())
  
}
  

  
fmri.view3d <- function(ttt, sigma=NULL,type = "data", col = grey(0:255/255), ext = 1, weights =
                        c(1,1,1), scale=c(0,1), scalecol = col,
                        hrf=rep(0,100), quant =3, maxpvalue = 0.05,pos=c(-1,-1,-1)) {
  # check wether Tk/Tcl environment is present
  if (!require(tkrplot))
    stop("required package tkrplot not found. Please install from cran.r-project.org")

  # some basic data properties
  dt <- dim(ttt)
  zlim <- range(ttt, na.rm = TRUE)
  label <- c("x", "y", "z", "t", "signal cut-off")

  # center position with Tcl objects
  if (pos[1] == -1) {
    pos <- c(round(dt[1:3])/2, 1, scale[1])
  } else {
    pos <- c(pos,1,scale[1])
  }
   	
  helpFunc <- function(a){
  	a <- tclVar()	
  }	

  posv <- lapply(pos, helpFunc) 


  fmri.image <- function(which, factor) {
    switch(which, x = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
        # plot image
        if (length(dt) == 4) {
          slice <- ttt[pos[1],,,pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        } else {
          slice <- ttt[pos[1],,]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[2],1:dt[3],slice, col=col, zlim=zlim)
        }
        # mark position
        lines(c(pos[2],pos[2]), c(0,dt[3])+0.5, col=2)
        lines(c(0,dt[2])+0.5, c(pos[3],pos[3]), col=2)
      }
    }, y = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
        # plot image
        if (length(dt) == 4) {
          slice <- ttt[,pos[2],,pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        } else {
          slice <- ttt[,pos[2],]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[3],slice, col=col, zlim=zlim)
        }
        # mark position
        lines(c(pos[1],pos[1]), c(0,dt[3])+0.5, col=2)
        lines(c(0,dt[1])+0.5, c(pos[3],pos[3]), col=2)
      }
    }, z = {
      f <- function() {
        oldpar <- par(mar=c(0,0,0,0))
        on.exit(par(oldpar))
        if (type == "spm") thresh <- (as.numeric(tclvalue(posv[[5]])) - scale[1])/diff(scale)  
        # plot image
        if (length(dt) == 4) {
          slice <- ttt[,dt[2]:1,pos[3],pos[4]]
          if (type == "spm") slice[slice<thresh] <- 0
          image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
        } else {
          slice <- ttt[,dt[2]:1,pos[3]]
          if (type == "spm") slice[slice<thresh] <- 0          
          image(1:dt[1],1:dt[2],slice, col=col, zlim=zlim)
        }
         # mark position
        lines(c(pos[1],pos[1]), c(0,dt[2])+0.5, col=2)
        lines(c(0,dt[1])+0.5, c(dt[2]-pos[2]+1,dt[2]-pos[2]+1), col=2)
      }
    })      
    # create the Tk-widget
    tkrplot(tt, f, hscale=ext, vscale=factor*ext)
  }

  fmri.slider <- function(i) {
    f <- function(...) {
      current <- as.numeric(tclvalue(posv[[i]]))
      if (current != pos[i]) {
        pos[i] <<- current
        tkrreplot(img[[1]])
        tkrreplot(img[[2]])
        tkrreplot(img[[3]])
        tkrreplot(img[[4]])
        if (i == 4) tkrreplot(img[[4]])
        tkconfigure(label2, text=pos[i])
      }
    }
    fr <- tkframe(tt)
    s <- tkscale(fr, command=f, from=1, to=dt[i], resolution=1, 
                 variable=posv[[i]], showvalue=FALSE, orient="horiz")
    label1 <- tklabel(fr, text=label[i])
    label2 <- tklabel(fr, text=pos[i])
    tkgrid(label1, s, label2)
    fr
  }

  fmri.threshold <- function(i) {
    f <- function(...) {
      current <- as.numeric(tclvalue(posv[[i]]))
      if (current != pos[i]) {
        pos[i] <<- current
        tkrreplot(img[[1]])
        tkrreplot(img[[2]])
        tkrreplot(img[[3]])
        tkconfigure(label2, text=pos[i])
      }
    }
    fr <- tkframe(tt)
    s <- tkscale(fr, command=f, from=scale[1], to=scale[2], resolution=diff(scale)/100, 
                 variable=posv[[i]], showvalue=FALSE, orient="horiz")
    label1 <- tklabel(fr, text=label[i])
    label2 <- tklabel(fr, text=pos[i])
    tkgrid(label1, s, label2)
    fr
  }

  fmri.scale <- function(which,scale=scale, scalecol=scalecol) {
    switch(which, "data" = {
      f <- function() {
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),width=c(200),height=c(160,40))
        on.exit(par(oldpar))
        # plot timeseries
        plot(ttt[pos[1],pos[2],pos[3],], xlab="Scan", ylab="BOLD signal")
        # mark scan number position
        lines(c(pos[4], pos[4]),range(ttt[pos[1],pos[2],pos[3],]),col=2)
        # draw scale
        image(seq(scale[1],scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
              matrix(rep(seq(scale[1],scale[2],length=100),10),100,10),
              yaxt="n",xlab="", ylab="",zlim=scale, col=scalecol)
      }
      # create the Tk-widget
      tkrplot(tt, f, hscale=ext, vscale=ext)
    }, "spm" = {
      f <- function() {
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),width=c(200),height=c(160,40))
        on.exit(par(oldpar))
        # draw something
        if (!is.null(sigma)) {
          value <- scale[1]+ttt[pos[1],pos[2],pos[3]]*diff(scale)          
          plot(c(1,length(hrf)),range(c(value*hrf,(value-3*sigma[pos[1],pos[2],pos[3]])*hrf,(value+3*sigma[pos[1],pos[2],pos[3]])*hrf)),type="n",xlab="Scan",ylab="Paramter estimate")
          xx <- c(1:length(hrf),length(hrf):1)
          yy <- c((value-quant*sigma[pos[1],pos[2],pos[3]])*hrf,rev((value+quant*sigma[pos[1],pos[2],pos[3]])*hrf))
          polygon(xx,yy,col="gray",lty=1)
          lines(value*hrf)
          lines(c(1,length(hrf)),c(0,0))
#          text(0.1,0.5,paste("Parameter:",signif(ttt[pos[1],pos[2],pos[3]],3)),pos=4,cex=1.5)
        } else {
          value <- scale[1]+ttt[pos[1],pos[2],pos[3]]*diff(scale)          
          plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
          text(0.1,0.5,paste("t-value:",signif(value,3)),pos=4,cex=1.5)
        }
        # draw scale
        image(seq(scale[1],scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
              matrix(rep(seq(scale[1],scale[2],length=100),10),100,10),
              yaxt="n",xlab="", ylab="", zlim=scale, col=scalecol)
        lines(c(value,value),scale,col="white")
      }
      # create the Tk-widget
      tkrplot(tt, f, hscale=ext, vscale=ext)
    }, "pvalue" = {
      f <- function() {
        if (ttt[pos[1],pos[2],pos[3]] <= 0.5) {
          value <- -log(maxpvalue)
        } else {
          value <- scale[1]+2*(ttt[pos[1],pos[2],pos[3]]-0.5)*diff(scale)          
        }
        oldpar <- par(mar=c(3,3,0.25,0.25), mgp=c(2,1,0))
        layout(matrix(1:2,2,1,byrow=TRUE),width=c(200),height=c(160,40))
        on.exit(par(oldpar))
        # draw something
        plot(c(0,1),c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
        if (value == -log(maxpvalue)) {
          text(0.2,0.5,paste("p-value: >",signif(exp(-value),3)),pos=4,cex=1.5)
        } else if (value == scale[2]) {
          text(0.2,0.5,paste("p-value: <",signif(exp(-value),3)),pos=4,cex=1.5)
        } else {
          text(0.2,0.5,paste("p-value:",signif(exp(-value),3)),pos=4,cex=1.5)
        }
        # draw scale
        image(seq(-log(maxpvalue),scale[2],length=100),seq(scale[1],scale[2],length=10)/10,
              matrix(rep(seq(-log(maxpvalue),scale[2],length=100),10),100,10),
              yaxt="n",xaxt="n",xlab="", ylab="",zlim=c(-log(maxpvalue),scale[2]), col=scalecol)
        lines(c(value,value),scale,col=1)
        lines(c(-log(0.01),-log(0.01)),scale,col=2)
        text(-log(0.01),scale[1]+0.01*diff(scale),pos=4,"1e-2")
        lines(c(-log(0.001),-log(0.001)),scale,col=2)
        text(-log(0.001),scale[1]+0.01*diff(scale),pos=4,"1e-3")
        lines(c(-log(0.0001),-log(0.0001)),scale,col=2)
        text(-log(0.0001),scale[1]+0.01*diff(scale),pos=4,"1e-4")
        lines(c(-log(0.00001),-log(0.00001)),scale,col=2)
        text(-log(0.00001),scale[1]+0.01*diff(scale),pos=4,"1e-5")
        lines(c(-log(0.000001),-log(0.000001)),scale,col=2)
        text(-log(0.000001),scale[1]+0.01*diff(scale),pos=4,"1e-6")
        lines(c(-log(0.0000001),-log(0.0000001)),scale,col=2)
        text(-log(0.0000001),scale[1]+0.01*diff(scale),pos=4,"1e-7")
        lines(c(-log(0.00000001),-log(0.00000001)),scale,col=2)
        text(-log(0.00000001),scale[1]+0.01*diff(scale),pos=4,"1e-8")
        lines(c(-log(0.000000001),-log(0.000000001)),scale,col=2)
        text(-log(0.000000001),scale[1]+0.01*diff(scale),pos=4,"1e-9")
      }
      # create the Tk-widget
      tkrplot(tt, f, hscale=ext, vscale=ext)
    })
  }
  
  # create window
  tt <- tktoplevel(bg="white")

  # create slider and images
  if (type == "data") {
    s <- lapply(1:4, fmri.slider)
  } else if (type == "spm") {
    s <- c(lapply(1:3, fmri.slider),fmri.threshold(5))
  } else {
    s <- lapply(1:3, fmri.slider)
  }
  img <- list(fmri.image("x",dt[3]/dt[1]*weights[3]),
              fmri.image("y",dt[3]/dt[2]*weights[3]),
              fmri.image("z",1),
              fmri.scale(type,scale,scalecol))

  # place the images and scales
  tkgrid(img[[2]], img[[1]])
  tkgrid(s[[2]], s[[1]])
  tkgrid(img[[3]], img[[4]])
  if (type == "data") {
    tkgrid(s[[3]], s[[4]])
  } else if (type == "spm") {
    tkgrid(s[[3]], s[[4]])
  } else {
    tkgrid(s[[3]])
  }
  
  # return the window object to the master
  tt
}

# select ROI from fmri dataset
cutroi <- function(data,
                    xind=1:data$dim[1],
                    yind=1:data$dim[2],
                    zind=1:data$dim[3],
                    tind=1:data$dim[4]) {

  if (("fmridata" %in% class(data)) & (!any(c("fmrispm","fmripvalue") %in% class(data)))) {
    ttt <- extract.data(data)[xind,yind,zind,tind]
    data$ttt <- writeBin(as.numeric(ttt),raw(),4)
    data$dim <- c(length(xind),length(yind),length(zind),length(tind))
    data$mask <- data$mask[xind,yind,zind]

    roixa <- (data$roixa:data$roixe)[xind[1]];
    roixe <- (data$roixa:data$roixe)[xind[length(xind)]];
    roiya <- (data$roiya:data$roiye)[yind[1]];
    roiye <- (data$roiya:data$roiye)[yind[length(yind)]];
    roiza <- (data$roiza:data$roize)[zind[1]];
    roize <- (data$roiza:data$roize)[zind[length(zind)]];
    roit <- data$roit[tind];

    data$roixa <- roixa
    data$roixe <- roixe
    data$roiya <- roiya
    data$roiye <- roiye
    data$roiza <- roiza
    data$roize <- roize
    data$roit <- roit
  }
  invisible(data)
}

# show a slice of pvalues with anatomical overlay!
# should this really use adimpro???
show.slice <- function(x, anatomic, maxpvalue = 0.05, slice = 1, view = "axial", col.u, col.o, zlim.u =
                    NULL, zlim.o = NULL) {

  pvalue <- x$pvalue
  pvalue[pvalue>0.05] <- 1
  pvalue[pvalue == 0] <- min(pvalue[pvalue>0])
  pvalue <- -log(pvalue)
  mask <- pvalue > 0

  ind2pos.ana <- conv.ip(anatomic, what="i2p")
  pos2ind.ana <- conv.ip(anatomic, what="p2i")
  ind2pos.func <- conv.ip(x, what="i2p")
  pos2ind.func <- conv.ip(x, what="p2i")

  pixdim.ana <- pixdim(anatomic$header,anatomic$format)
  pixdim.func <- pixdim(x$header,x$format)

  ttt.ana <- extract.data(anatomic)
  ddim.ana <- dim(ttt.ana) <- dim(ttt.ana)[1:3]

  if (view == "axial") {
    dfunc <- dim(pvalue)[1:2]
    if ((slice >= 1) & (slice <= dim(pvalue)[3])) {
      imgdata.o <- pvalue[,,slice]
      mask <- mask[,,slice]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[1:2]))/min(abs(pixdim.ana)))
  } else if (view == "coronal") {
    dfunc <- dim(pvalue)[c(1,3)]
    if ((slice >= 1) & (slice <= dim(pvalue)[2])) {
      imgdata.o <- pvalue[,slice,]
      mask <- mask[,slice,]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[c(1,3)]))/min(abs(pixdim.ana)))
  } else if (view == "sagittal") {
    dfunc <- dim(pvalue)[c(2,3)]
    if ((slice >= 1) & (slice <= dim(pvalue)[1])) {
      imgdata.o <- pvalue[slice,,]
      mask <- mask[slice,,]
    } else {
      mask <- imgdata.o <- array(0,dim=dfunc)
    }
    scale <- ceiling(max(abs(pixdim.func[2:3]))/min(abs(pixdim.ana)))
  } else {
    stop("unknown view",view)
  }

  imgdata.n <- array(0,dim=c(scale*dim(imgdata.o)))
  mask.n <- array(FALSE,dim=c(scale*dim(imgdata.o)))
  for (i in 1:dim(imgdata.o)[1]) {
    for (j in 1:dim(imgdata.o)[2]) {
      imgdata.n[(i-1)*scale+c(1:scale),(j-1)*scale+c(1:scale)] <- imgdata.o[i,j]
      mask.n[(i-1)*scale+c(1:scale),(j-1)*scale+c(1:scale)] <- imgdata.o[i,j]
    }
  }
  imgdata.o <- imgdata.n
  mask <- mask.n

  imgdata.u <- array(0, dim=dfunc*scale)
  for (i in 1:(dfunc[1]*scale)) {
    for (j in 1:(dfunc[2]*scale)) {
      if (view == "axial") {
        pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5, x$roiya+(2*j-1)/(2*scale)-0.5, x$roiza + slice - 1) )
      } else if (view == "coronal") {
        pos <- ind2pos.func( c(x$roixa+(2*i-1)/(2*scale)-0.5, x$roiya + slice - 1, x$roiza+(2*j-1)/(2*scale)-0.5) )
      } else if (view == "sagittal") {
        pos <- ind2pos.func( c(x$roixa + slice -1, x$roiya+(2*i-1)/(2*scale)-0.5, x$roiza+(2*j-1)/(2*scale)-0.5) )
      }
      ind.ana <- pos2ind.ana(pos) # this is real(!) index for anatomic image
      ii <- ind.ana[1]
      jj <- ind.ana[2]
      kk <- ind.ana[3]

      iint <- ceiling(ind.ana[1]) # these are the integer indices
      jint <- ceiling(ind.ana[2])
      kint <- ceiling(ind.ana[3])

#       if ((iint-1 >= 1) & (jint -1 >= 1) & (kint -1 >= 1) &
#           (iint <= ddim.ana[1]) & (jint <= ddim.ana[2]) & (kint <= ddim.ana[3])) {
#         imgdata.u[i,j] <-
#           ttt.ana[iint-1,jint-1,kint-1] * (iint - ii) * (jint - jj) * (kint - kk) +
#             ttt.ana[iint-1,jint,kint-1] * (iint - ii) * (jj - jint + 1) * (kint - kk) +
#               ttt.ana[iint,jint-1,kint-1] * (ii - iint + 1) * (jint - jj) * (kint - kk) +
#                 ttt.ana[iint,jint,kint-1] * (ii - iint + 1) * (jj - jint + 1) * (kint - kk) +
#                   ttt.ana[iint-1,jint-1,kint] * (iint - ii) * (jint - jj) * (kk - kint + 1) +
#                     ttt.ana[iint-1,jint,kint] * (iint - ii) * (jj - jint + 1) * (kk - kint + 1) +
#                       ttt.ana[iint,jint-1,kint] * (ii - iint + 1) * (jint - jj) * (kk - kint + 1) +
#                         ttt.ana[iint,jint,kint] * (ii - iint + 1) * (jj - jint + 1) * (kk - kint + 1) 
#       }
      if ((iint >= 1) & (jint >= 1) & (kint >= 1) &
          (iint <= ddim.ana[1]) & (jint <= ddim.ana[2]) & (kint <= ddim.ana[3])) {
        imgdata.u[i,j] <- ttt.ana[iint,jint,kint] * (ii - iint + 1) * (jj - jint + 1) * (kk - kint + 1)
        if (kint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint,kint-1] * (ii - iint + 1) * (jj - jint + 1) * (kint - kk)
        if (jint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint-1,kint] * (ii - iint + 1) * (jint - jj) * (kk - kint + 1)
        if (iint > 1) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint,kint] * (iint - ii) * (jj - jint + 1) * (kk - kint + 1)
        if ((iint > 1) & (jint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint-1,kint] * (iint - ii) * (jint - jj) * (kk - kint + 1)
        if ((iint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint,kint-1] * (iint - ii) * (jj - jint + 1) * (kint - kk)
        if ((jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint,jint-1,kint-1] * (ii - iint + 1) * (jint - jj) * (kint - kk)
        if ((iint > 1) & (jint > 1) & (kint > 1)) imgdata.u[i,j] <- imgdata.u[i,j] + ttt.ana[iint-1,jint-1,kint-1] * (iint - ii) * (jint - jj) * (kint - kk)
      }
    }
  }

  if (is.null(zlim.o)) {
    zlim.o <- range(imgdata.o)
  } else {
    if (length(zlim.o) != 2) stop("zlim.o not length 2")
    if (zlim.o[2] < zlim.o[1]) stop("zlim.o[2] < zlim.o[1]")
    imgdata.o[imgdata.o > zlim.o[2]] <- zlim.o[2]
    imgdata.o[imgdata.o < zlim.o[1]] <- zlim.o[1]
  }
  if (is.null(zlim.u)) {
    zlim.u <- range(imgdata.u)
  } else {
    if (length(zlim.u) != 2) stop("zlim.u not length 2")
    if (zlim.u[2] < zlim.u[1]) stop("zlim.u[2] < zlim.u[1]")
    imgdata.u[imgdata.u > zlim.u[2]] <- zlim.u[2]
    imgdata.u[imgdata.u < zlim.u[1]] <- zlim.u[1]
  }
  
  img <- array(0, dim=c(dim(imgdata.u),3))
  for (i in 1:dim(imgdata.u)[1]) {
    for (j in 1:dim(imgdata.u)[2]) {
      if (mask[i,j]) { # use overlay
        level <- length(col.o) * (imgdata.o[i,j] - zlim.o[1]) / diff(zlim.o)
        level <- ceiling(level) # now in 0:length(col.o)
        if (is.na(level)) level <- 1
        if (level == 0) level <- 1 # now in 1:length(col.o)
        img[i,j,] <- as.integer(col2rgb(col.o[level])) * 256
      } else { # use underlay
        level <- length(col.u) * (imgdata.u[i,j] - zlim.u[1]) / diff(zlim.u)
        level <- ceiling(level) # now in 0:length(col.u)
        if (is.na(level)) level <- 1
        if (level == 0) level <- 1 # now in 1:length(col.u)
        img[i,j,] <- as.integer(col2rgb(col.u[level])) * 256

      }
    }
  }

  invisible(img)
}

pixdim <- function(header,format) {
  if (format == "NIFTI") {
    return(header$pixdim[2:4])
  } else if (format == "ANALYZE") {
    return(header$pixdim[2:4])
  } else if (format == "HEAD/BRIK") {
    return(header$DELTA)
  } else {
    stop("Not implemented for this data format:", format)
  }
}

conv.ip <- function(data, what="i2p") {
  if (!("fmridata" %in% class(data))) stop("Cannot evaluate real-space position for this dataset. Not type fmridata!")

  if (data$format == "NIFTI") {

    if (data$header$qform > 0) {
      origin <- c(data$header$qoffsetx, data$header$qoffsety, data$header$qoffsetz)
      b <- data$header$quaternb
      c <- data$header$quaternc
      d <- data$header$quaternd
      a <- sqrt(pmax(0,1-b*b-c*c-d*d))
      R <- t(matrix(c(a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c,
                    2*b*c+2*a*d, a*a+c*c-b*b-d*d, 2*c*d -2*a*b,
                    2*b*d-2*a*c, 2*c*d+2*a*b, a*a+d*d-c*c-b*b),3,3))
      pixdim <- data$header$pixdim[2:4]
      qfac <- data$header$pixdim[1]

      if (what == "i2p") {
        return(function(ind) R %*% (c(1,1,qfac) * pixdim * (ind-1)) + origin)
      } else {
        return(function(pos) (solve(R) %*% (pos - origin))/(c(1,1,qfac) * pixdim) + 1)
      }
    } else if (data$header$sform > 0) {
      origin <- c(data$header$srowx[4],data$header$srowy[4],data$header$srowz[4])
      SR <- matrix(c(data$header$srowx[1],data$header$srowy[1],data$header$srowz[1],
                     data$header$srowx[2],data$header$srowy[2],data$header$srowz[2],
                     data$header$srowx[3],data$header$srowy[3],data$header$srowz[3]),3,3)
      if (what == "i2p") {
        return(function(ind) SR %*% (ind-1) + origin)
      } else {
        return(function(pos) solve(SR) %*% (pos - origin) + 1)
      }
    } else if (data$header$qform == 0) {
      warning("This method is specified only for compatibility reasons to ANALYZE 7.5. May not deliver useful results")

      pixdim <- data$header$pixdim[2:4]
      if (what == "i2p") {
        return(function(ind) pixdim * (ind-1))
      } else {
        return(function(pos) pos/pixdim + 1)
      }
    } else {
      stop("Neither Method 1, 2, nor 3 for real-space position applicable. See NIFTI specification!")
    }

  } else if (data$format == "ANALYZE") {

    stop("Not yet implemented real-space position evaluation for this data format:", data$format)

  } else if (data$format == "HEAD/BRIK") {

    orientation <- data$header$ORIENT_SPECIFIC
    if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation",orientation,"found! \n")

    rxyz <- (orientation)%/%2+1
    xyz <- rxyz[rxyz]
    pixdim <- data$header$DELTA[xyz]
    origin <- data$header$ORIGIN[xyz]
    if (what == "i2p") {
      return(function(ind) pixdim * (ind[xyz]-1) + origin)
    } else {
      return(function(pos) ((pos-origin)/pixdim + 1)[rxyz])
    }

  } else if (data$format == "DICOM") {

    stop("Not yet implemented real-space position evaluation for this data format:", data$format)

  } else {

    stop("Not implemented real-space position evaluation for this data format (not in fmri package):", data$format)

  }
}

summary.fmridata <- function(object,...) {
  if ("fmripvalue" %in% class(object)) {
    dt <- dim(object$pvalue)
    cat("Data Dimension  :", dt,"\n")
    values <- range(object$pvalue)
    cat("Data Range      :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(object, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(object, "design")), "\n")
    switch(attr(object, "white"),cat("Prewhitening performed with smoothed map\nof autocorrelation parameter in AR(1) model for time series!\n"),
                                 cat("Prewhitening performed with map of autocorrelation parameter in AR(1) model for time series\n"),
                                 cat("No prewhitening performed!\n"))
    if (!is.null(attr(object, "smooth"))) cat(attr(object, "smooth"),"\n")
    cat(attr(object, "mode"), "\n")
    invisible(list(dim=dt,values=values, files=attr(object, "read"),
                   z=attr(object, "design")))
  } else if ("fmrispm" %in% class(object)) {
    dt <- object$dim
    cat("Data Dimension  :", dt,"\n")
    values <- range(object$cbeta)
    cat("Data Range      :", values[1], "...", values[2], "\n")
    cat("File(s)         :", attr(object, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(object, "design")), "\n")
    switch(attr(object, "white"),cat("Prewhitening performed with smoothed map\nof autocorrelation parameter in AR(1) model for time series!\n"),
                                 cat("Prewhitening performed with map of autocorrelation parameter in AR(1) model for time series\n"),
                                 cat("No prewhitening performed!\n"))
    if (!is.null(attr(object, "smooth"))) cat(attr(object, "smooth"))
    invisible(list(dim=dt,values=values, files=attr(object, "read"),
              z=attr(object, "design")))
  } else {
    dt <- object$dim
    cat("Data Dimension:", dt,"\n")
    values <- range(extract.data(object))
    cat("Data Range    :", values[1], "...", values[2], "\n")
    delta <- object$delta
    cat("Voxel Size    :", delta,"\n")
    cat("File(s)", attr(object, "file"),"\n")
    invisible(list(dim=dt,delta=delta,values=values, files=attr(object, "read")))
  }

}

print.fmridata <- function(x,...) {
  if ("fmripvalue" %in% class(x)) {
    cat("Data Dimension:", dim(x$pvalue),"\n")
    values <- range(x$pvalue)
    cat("Data Range    :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(x, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(x, "design")), "\n")
    switch(attr(x, "white"),cat("Prewhitening performed with smoothed map\nof autocorrelation parameter in AR(1) model for time series!\n"),
                                 cat("Prewhitening performed with map of autocorrelation parameter in AR(1) model for time series\n"),
                                 cat("No prewhitening performed!\n"))
    if (!is.null(attr(x, "smooth"))) cat(attr(x, "smooth"),"\n")
    cat(attr(x, "mode"), "\n")
  } else if ("fmrispm" %in% class(x)) {
    cat("Data Dimension:", x$dim,"\n")
    values <- range(x$cbeta)
    cat("Data Range    :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(x, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(x, "design")), "\n")
    switch(attr(x, "white"),cat("Prewhitening performed with smoothed map\nof autocorrelation parameter in AR(1) model for time series!\n"),
                                 cat("Prewhitening performed with map of autocorrelation parameter in AR(1) model for time series\n"),
                                 cat("No prewhitening performed!\n"))
    if (!is.null(attr(x, "smooth"))) cat(attr(x, "smooth"))
#    lmcall <- attr(x, "lm")
#    cat("Linear Model - Number of stimuli
  } else {
    cat("Data Dimension: ", x$dim,"\n")
    cat("Voxel Size    :", x$delta,"\n")
    values <- range(extract.data(x))
    cat("Data Range    :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(x, "file"),"\n")
  }
  invisible(NULL)
}
