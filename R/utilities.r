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
