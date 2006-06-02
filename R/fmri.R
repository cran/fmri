fmri.smooth <- function(spm,hmax=4,adaptive=TRUE,lkern="Triangle",skern="Triangle") {
  cat("fmri.smooth: entering function\n")
  
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
  if (is.null(spm$scorr)) {
    scorr <- 0
  } else {
    scorr <- spm$scorr
  }

  cat("fmri.smooth: smoothing the Statistical Paramteric Map\n")
  if (adaptive) {
    ttthat <- vaws3D(y=spm$cbeta, sigma2=variance, hmax=hmax,
                     wghts=weights, scorr=scorr, vwghts = spm$vwghts,
                     lkern=lkern,skern=skern)
  } else {
    ttthat <- vaws3D(y=spm$cbeta, sigma2=variance, hmax=hmax,
                     qlambda = 1, wghts=weights, scorr=scorr,
                     vwghts = spm$vwghts,lkern=lkern,skern=skern)
  }
  cat("\n")
  
  cat("fmri.smooth: determine local smoothness\n")
  bw <- get3Dh.gauss(ttthat$vred,weights)
  rxyz <- c(resel(1,bw[,1]), resel(1,bw[,2]), resel(1,bw[,3]))
  dim(rxyz) <- c(dim(bw)[1],3)
  bw0 <- get3Dh.gauss(ttthat$vred0,weights)
  rxyz0 <- c(resel(1,bw0[,1]), resel(1,bw0[,2]), resel(1,bw0[,3]))
  dim(rxyz0) <- c(dim(bw0)[1],3)
  cat("fmri.smooth: exiting function\n")
    
  if (dim(ttthat$theta)[4] == 1) {
    z <- list(cbeta = ttthat$theta[,,,1], var = ttthat$var, rxyz =
              rxyz, rxyz0 = rxyz0, scorr = spm$scorr, weights =
              spm$weights, vwghts = spm$vwghts,
              hmax = ttthat$hmax, dim = spm$dim, hrf = spm$hrf)
  } else {
    z <- list(cbeta = ttthat$theta, var = ttthat$var, rxyz = rxyz, rxyz0 = rxyz0, 
              scorr = spm$scorr, weights = spm$weights, vwghts = spm$vwghts,
              hmax = ttthat$hmax, dim = spm$dim, hrf = spm$hrf)
  }    

  class(z) <- c("fmridata","fmrispm")

  attr(z, "file") <- attr(spm, "file")
  attr(z, "white") <- attr(spm, "white")
  attr(z, "design") <- attr(spm, "design")

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

fmri.pvalue <- function(spm, mode="basic", delta=NULL) {
  cat("fmri.pvalue: entering function\n")

  if (!("fmrispm" %in% class(spm)) ) {
    warning("fmri.pvalue: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (length(dim(spm$cbeta)) < 4) {

    stat <- spm$cbeta/sqrt(spm$var)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm")
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm")
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm")
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm")
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }        
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm")
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm")
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
  mask <- rep(1,length=prod(spm$dim[1:3]))
  mask[stat < thresh] <- 0
  pv[!mask] <- 1
  dim(pv) <- spm$dim[1:3]
  
  cat("fmri.pvalue: exiting function\n")

  z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim, hrf = spm$hrf)
  
  class(z) <- c("fmridata","fmripvalue")

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




plot.fmridata <- function(x, anatomic = NULL , maxpvalue =
                          0.05, spm = TRUE,
                            pos = c(-1,-1,-1), type="slice",
                            device="X11", file="plot.png",...) {
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

    quant <- if (x$smooth) qnorm(1-maxpvalue) else qt(1-maxpvalue,length(hrf))
    
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
    signal <- x$ttt
    
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
    pos <- c(pos,1,scal[1])
  }
  posv <- lapply(pos, tclVar)


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

summary.fmridata <- function(object,...) {
  if ("fmripvalue" %in% class(object)) {
    dt <- dim(object$pvalue)
    cat("Data Dimension  :", dt,"\n")
    values <- range(object$pvalue)
    cat("Data Range      :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(object, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(object, "design")), "\n")
    cat(attr(object, "white"), "\n")
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
    cat(attr(object, "white"), "\n")
    if (!is.null(attr(object, "smooth"))) cat(attr(object, "smooth"))
    invisible(list(dim=dt,values=values, files=attr(object, "read"),
              z=attr(object, "design")))
  } else {
    dt <- object$dim
    cat("Data Dimension:", dt,"\n")
    values <- range(object$ttt)
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
    cat(attr(x, "white"), "\n")
    if (!is.null(attr(x, "smooth"))) cat(attr(x, "smooth"),"\n")
    cat(attr(x, "mode"), "\n")
  } else if ("fmrispm" %in% class(x)) {
    cat("Data Dimension:", x$dim,"\n")
    values <- range(x$cbeta)
    cat("Data Range    :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(x, "file"),"\n\n")
    cat("Design Dimension:", dim(attr(x, "design")), "\n")
    cat(attr(x, "white"), "\n")
    if (!is.null(attr(x, "smooth"))) cat(attr(x, "smooth"))
#    lmcall <- attr(x, "lm")
#    cat("Linear Model - Number of stimuli
  } else {
    cat("Data Dimension: ", x$dim,"\n")
    cat("Voxel Size    :", x$delta,"\n")
    values <- range(x$ttt)
    cat("Data Range    :", values[1], "...", values[2], "\n")
    cat("File(s)", attr(x, "file"),"\n")
  }
  invisible(NULL)
}
