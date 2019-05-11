fwhm2bw <- function(hfwhm) hfwhm/sqrt(8*log(2))
bw2fwhm <- function(h) h*sqrt(8*log(2))

getvofh <- function(bw,lkern,wght){
.Fortran(C_getvofh,
         as.double(bw),
         as.integer(lkern),
         as.double(wght),
         vol=double(1))$vol
}
gethani <- function(x,y,lkern,value,wght,eps=1e-2){
.Fortran(C_gethani,
         as.double(x),
         as.double(y),
         as.integer(lkern),
         as.double(value),
         as.double(wght),
         as.double(eps),
         bw=double(1))$bw
}

fdr <- function(pval, alpha=0.05){
  n <- length(pval)
  ind <- 1:n
  oind <- order(pval)
  nind <- length(ind[pval[oind] <= ind/n*alpha])
  oind[1:nind]
}

Varcor.gauss<-function(h,interv = 1){
#
#   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
#
#   in case of colored noise that was produced by smoothing with lkern and bandwidth h
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing,
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding
#  discretisation into voxel)
#
#  h in FWHM
#
  h<-fwhm2bw(h)*interv
  ih<-trunc(4*h)+1
  dx<-2*ih+1
  d<-length(h)
  penl <- dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl <- outer(penl,dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl <- outer(penl,outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  2*sum(penl)^2/sum(diff(penl,interv)^2)/interv^(d)
}

Spatialvar.gauss<-function(h,h0,d,interv=1){
#
#   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in
#
#   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
#
#   Spatialvar.gauss(lkern,h,h0,d)/Spatialvar.gauss(lkern,h,1e-5,d) gives the
#   a factor for lambda to be used with bandwidth h
#
#
#  interv allows for further discretization of the Gaussian Kernel, result depends on
#  interv for small bandwidths. interv=1  is correct for kernel smoothing,
#  interv>>1 should be used to handle intrinsic correlation (smoothing preceeding
#  discretisation into voxel)
#
#  h and h0 in FWHM
#
  h0 <- pmax(h0,1e-5)
  h <- pmax(h,1e-5)
  h<-fwhm2bw(h)*interv
  if(length(h)==1) h<-rep(h,d)
  ih<-trunc(4*h)
  ih<-pmax(1,ih)
  dx<-2*ih+1
  penl<-dnorm(((-ih[1]):ih[1])/h[1])
  if(d==2) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
  if(d==3) penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),outer(dnorm(((-ih[2]):ih[2])/h[2]),dnorm(((-ih[3]):ih[3])/h[3]),"*"),"*")
  dim(penl)<-dx
  h0<-fwhm2bw(h0)*interv
  if(length(h0)==1) h0<-rep(h0,d)
  ih<-trunc(4*h0)
  ih<-pmax(1,ih)
  dx0<-2*ih+1
  x<- ((-ih[1]):ih[1])/h0[1]
  penl0<-dnorm(((-ih[1]):ih[1])/h0[1])
  if(d==2) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
  if(d==3) penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),outer(dnorm(((-ih[2]):ih[2])/h0[2]),dnorm(((-ih[3]):ih[3])/h0[3]),"*"),"*")
  dim(penl0)<-dx0
  penl0<-penl0/sum(penl0)
  dz<-dx+dx0-1
  z<-array(0,dz)
  if(d==1){
    for(i1 in 1:dx0) {
      ind1<-c(0:(i1-1),(dz-dx0+i1):dz+1)
      ind1<-ind1[ind1<=dz][-1]
      z[-ind1]<-z[-ind1]+penl*penl0[i1]
    }
  } else if(d==2){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
    }
  } else if(d==3){
    for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]) for(i3 in 1:dx0[3]){
      ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
      ind1<-ind1[ind1<=dz[1]][-1]
      ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
      ind2<-ind2[ind2<=dz[2]][-1]
      ind3<-c(0:(i3-1),(dz[3]-dx0[3]+i3):dz[3]+1)
      ind3<-ind3[ind3<=dz[3]][-1]
      z[-ind1,-ind2,-ind3]<-z[-ind1,-ind2,-ind3]+penl*penl0[i1,i2,i3]
    }
  }
  sum(z^2)/sum(z)^2*interv^d
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


get.corr.gauss <- function(h,interv=1) {
    #
    #   Calculates the correlation of
    #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
    #   Result does not depend on d for "Gaussian" kernel !!
    #
    #   h in FWHM
    #
    h <- fwhm2bw(h)*interv
    ih <- trunc(4*h+ 2*interv-1)
    dx <- 2*ih+1
    penl <- dnorm(((-ih):ih)/h)
    sum(penl[-(1:interv)]*penl[-((dx-interv+1):dx)])/sum(penl^2)
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
                 as.logical(mask),
                 as.integer(n1),
                 as.integer(n2),
                 as.integer(n3),
                 as.integer((n1+1)/2),
                 as.integer((n2+1)/2),
                 as.integer((n3+1)/2),
                 integer(n),
                 integer(n),
                 integer(n),
                 logical(n),
                 mask=logical(n))$mask
dim(mask1) <- dm
mask1
}

spatial.corr <- function(residuals){
  lags <- c(5,5,3)
  dy <- dim(residuals)
  mask <- array(TRUE,dy[2:4])
  corr <- .Fortran(C_mcorr,as.double(residuals),
                     as.logical(mask),
                     as.integer(dy[2]),
                     as.integer(dy[3]),
                     as.integer(dy[4]),
                     as.integer(dy[1]),
                     scorr=double(prod(lags)),
                     as.integer(lags[1]),
                     as.integer(lags[2]),
                     as.integer(lags[3]))$scorr
  dim(corr) <- lags
  bw <- optim(c(2,2,2),corrrisk,method="L-BFGS-B",lower=c(.25,.25,.25),upper=c(20,20,20),lag=lags,data=corr)$par
  bw[bw<=.25] <- 0
  if( (max(bw) > 2.5 ) || (corr[lags[1],1,1]+corr[1,lags[2],1]+corr[1,1,lags[3]] >0.5) ) warning(paste("Local smoothness characterized by large bandwidth ",bw," check residuals for structure",collapse=","))
  rxyz <- c(resel(1,bw[1]), resel(1,bw[2]), resel(1,bw[3]))
  dim(rxyz) <- c(1,3)
list(scorr=corr,bw=bw)
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
            dim0 = ddim,
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
   if(class(fmriobj)!="fmridata") stop("setmask 1st argument needs to be of class fmridata\n")
   if(!is.array(mask)) stop("setmask 2nd argument needs to be an array\n")
   if(any(fmriobj$dim[1:3]!=dim(mask))) stop("setmask incompatible dimensions of arguments\n")
   fmriobj$mask <- mask
}
