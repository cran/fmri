fwhm2bw <- function(hfwhm) hfwhm/sqrt(8*log(2))
bw2fwhm <- function(h) h*sqrt(8*log(2))


fdr <- function(pval, alpha=0.05){
  n <- length(pval)
  ind <- 1:n
  oind <- order(pval)
  nind <- length(ind[pval[oind] <= ind/n*alpha])
  oind[1:nind]
}



get3Dh.gauss <- function(vred,vwghts,step=1.002,interv=1){
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing,
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding
#  discretisation into voxel)
#
  ## NOTE this function requires hvred in fmri/R/sysdata.rda
  n <- length(vred)
  hv <- matrix(0,3,n)
  fixed <- rep(FALSE,n)
  i <- 1
  while (any(!fixed) && i < dim(hvred)[1]) {
    ind<-(1:n)[!fixed][vred[!fixed] >= hvred[i,2]]
    hv[,ind] <- hvred[i,1]
    fixed[ind] <- TRUE
    i <- i + 1
  }
  hv[,!fixed] <- hvred[i,1]

  invisible(t(hv))
}



corrrisk <- function(bw,lag,data){
  mean((data-thcorr3D(bw,lag))^2)
}

thcorr3D <- function(bw,lag=rep(5,3)){
  g <- trunc(fwhm2bw(bw)*4)
  gw1 <- dnorm(-(g[1]):g[1],0,fwhm2bw(bw[1]))
  gw2 <- dnorm(-(g[2]):g[2],0,fwhm2bw(bw[2]))
  gw3 <- dnorm(-(g[3]):g[3],0,fwhm2bw(bw[3]))
  gwght <- outer(gw1,outer(gw2,gw3,"*"),"*")
  gwght <- gwght/sum(gwght)
  dgw <- dim(gwght)
  scorr <- .Fortran(C_thcorr,as.double(gwght),
                    as.integer(dgw[1]),
                    as.integer(dgw[2]),
                    as.integer(dgw[3]),
                    scorr=double(prod(lag)),
                    as.integer(lag[1]),
                    as.integer(lag[2]),
                    as.integer(lag[3]),
                    PACKAGE="fmri")$scorr
  # bandwidth in FWHM in voxel units
  dim(scorr) <- lag
  scorr
}

connect.mask <- function(mask){
dm <- dim(mask)
n1 <- dm[1]
n2 <- dm[2]
n3 <- dm[3]
n <- n1*n2*n3
mask1 <- .Fortran(C_lconnect,
                 as.integer(mask),
                 as.integer(n1),
                 as.integer(n2),
                 as.integer(n3),
                 as.integer((n1+1)/2),
                 as.integer((n2+1)/2),
                 as.integer((n3+1)/2),
                 integer(n),
                 integer(n),
                 integer(n),
                 integer(n),
                 mask=integer(n))$mask
mask1 <- as.logical(mask1)
dim(mask1) <- dm
mask1
}


##
##  create fmri-object from nifti
##
niftiImage2fmri <- function(niftiobj, level=0.75, mask=NULL, setmask=TRUE,
                            indx=NULL,indy=NULL,indz=NULL,avoidnegs=FALSE) {
  ## Convert niftiImage object to "fmridata" S3 object
  ddim <- attributes(niftiobj)$dim
  pixdim <- attributes(niftiobj)$pixdim
  pixunits <- attributes(niftiobj)$pixunits
  dx <- ddim[1]
  dy <- ddim[2]
  dz <- ddim[3]
  dt <- ddim[4]
  if(is.null(indx)) indx <- 1:dx
  if(is.null(indy)) indy <- 1:dy
  if(is.null(indz)) indz <- 1:dz
  dx <- length(indx)
  dy <- length(indy)
  dz <- length(indz)
  ddim[1:3] <- c(dx,dy,dz)
  weights <- abs(pixdim[1:3])/min(abs(pixdim[1:3]))
  ttt <- array(0,ddim)
  ttt[,,,] <- niftiobj[indx,indy,indz,]
  if(avoidnegs&any(ttt<0)){
    ttt <- ttt-min(ttt)
    warning("changed mean of data to avoid negatives")
  }
  mask0 <- array(TRUE, c(dx,dy,dz))
  if (setmask) {
    mask0[ttt[,,,1] < quantile(ttt[,,,1], level, na.rm=TRUE)] <- FALSE
    dim(ttt) <- c(prod(dim(ttt)[1:3]), dim(ttt)[4])
    na <- ttt %*% rep(1, dim(ttt)[2])
    mask0[is.na(na)] <- FALSE
    ttt[is.na(na), ] <- 0
    dim(mask0) <- c(dx, dy, dz)
    mask0 <- connect.mask(mask0)|mask0
  }
  datascale <- max(abs(ttt))/32767
  ttt <- as.integer(ttt/datascale)
  z <- list(ttt = writeBin(ttt, raw(), 2),
            format = "NIFTI",
            delta = rep(1,3),
            origin = NULL,
            orient = NULL,
            dim = ddim,
            datascale = datascale,
            roixa = 1,
            roixe = dx,
            roiya = 1,
            roiye = dy,
            roiza = 1,
            roize = dz,
            roit = 1,
            weights = weights,
            header = list(NULL),
            mask = mask0)
  class(z) <- "fmridata"
  attr(z, "file") <- ""
  if(!is.null(mask)) z <- setmask(z,mask)
  return(z)
}

setmask <- function(fmriobj, mask){
   if(!inherits(fmriobj, "fmridata")) stop("setmask 1st argument needs to be of class fmridata\n")
   if(!is.array(mask)) stop("setmask 2nd argument needs to be an array\n")
   if(any(fmriobj$dim[1:3]!=dim(mask))) stop("setmask incompatible dimensions of arguments\n")
   fmriobj$mask <- mask
}
