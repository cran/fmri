fmri.sICA <- function(data, mask=NULL, ncomp=20,
               alg.typ=c("parallel","deflation"),
               fun=c("logcosh","exp"), alpha=1, detrend=TRUE,
               degree=2, nuisance= NULL, ssmooth=TRUE, tsmooth=TRUE, bwt=4,
               bws=8, unit=c("FWHM","SD")){
  if(detrend) data <- fmri.detrend(data,degree,nuisance)
  bwvoxel <- if(length(data$header$pixdim)>=4)
        bws/data$header$pixdim[2:4] else bws
  bwtime <- if(length(data$header$pixdim)>=5)
        bwt/data$header$pixdim[5] else bwt
  if(ssmooth) data <- smooth.fmridata(data,bwvoxel,unit)
  if(tsmooth) data <- smooth.fmridata(data,bwtime,unit,what="temporal")
  ttt <- extractData(data)
  if(is.null(mask)) mask <- data$mask
  if(!is.logical(mask)) mask <- as.logical(mask)
  ddim <- data$dim[1:3]
  dim(mask) <- ddim
  dim(ttt) <- c(prod(ddim),dim(ttt)[4])
  ttt0 <- ttt[mask,]
  if(!requireNamespace("fastICA",quietly=TRUE)){
     stop("package fastICA needed for this functionality, please install")
  }
  cat("Computing independent components with fastICA \n")
  set.seed(1)
  fsica <- fastICA::fastICA(ttt0,n.comp=ncomp,
       alg.typ=alg.typ, fun=fun, alpha=alpha, method="C",
       maxit=500, tol=1e-5)
  cimgs <- array(0,c(prod(ddim),ncomp))
  cimgs[mask,] <- fsica$S
  dim(cimgs) <- c(ddim,ncomp)
  z <- list(scomp=cimgs,X=fsica$X,k=fsica$K,W=fsica$W,A=fsica$A,
       mask=mask,voxdim=data$header$pixdim[2:4],TR=data$header$pixdim[5])
  class(z) <- "fmriICA"
  z
  }

ICAfingerprint <- function(icaobj,nbin=256,plot=FALSE){
##
##  Martino-Goebel(2006)
##
   scomp <- icaobj$scomp
   dscomp <- dim(scomp)
   ncomp <- dscomp[4]
   dim(scomp) <- c(prod(dscomp[1:3]),ncomp)
   scomp <- scomp[icaobj$mask,]
   nvox <- dim(scomp)[1]
   tcomp <- icaobj$A
   ntime <- dim(tcomp)[2]
   icafp <- matrix(0,11,ncomp)
   dimnames(icafp) <- list(c("kurt","skew","sentropy","dclust",
                        "ar1","tentropy","power1","power2",
                        "power3","power4","power5"),
                        paste0("Component",1:ncomp))
   for(i in 1:ncomp){
      sdcompi <- sd(scomp[,i])
      scomp[,i] <- scomp[,i]/sdcompi
      icafp[1,i] <- sum(scomp[,i]^4)/nvox - 3
      icafp[2,i] <- sum(scomp[,i]^3)/nvox
      z <- hist(scomp[,i],nbin,plot=FALSE)$counts
      z <- z[z>0]/sum(z)
      icafp[3,i] <- -sum(z*log(z,2))
      scompi <- icaobj$scomp[,,,i]
      scompi <- abs(scompi-mean(scompi))/sdcompi
      z <- findclusters(scompi,2.5)
      icafp[4,i] <- sum(z$size>270/prod(icaobj$voxdim))/nvox
      icafp[5,i] <- abs(sum(tcomp[i,-1]*tcomp[i,-ntime]))/sum(tcomp[i,]^2)*ntime/(ntime-1)
      z <- hist(tcomp[i,],nbin,plot=FALSE)$counts
      z <- z[z>0]/sum(z)
      icafp[6,i] <- -sum(z*log(z,2))
      z <- spectrum(tcomp[i,],plot=FALSE)
      ns <- length(z$freq)
      zfreq <- z$freq/icaobj$TR
      ind1 <- (1:ns)[zfreq<=.008]
      ind2 <- (1:ns)[zfreq>.008&zfreq<=.02]
      ind3 <- (1:ns)[zfreq>.02&zfreq<=.05]
      ind4 <- (1:ns)[zfreq>.05&zfreq<=.1]
      ind5 <- (1:ns)[zfreq>.1&zfreq<=.25]
      sumsp <- sum(z$spec)
      icafp[7,i] <- sum(z$spec[ind1])/sumsp
      icafp[8,i] <- sum(z$spec[ind2])/sumsp
      icafp[9,i] <- sum(z$spec[ind3])/sumsp
      icafp[10,i] <- sum(z$spec[ind4])/sumsp
      icafp[11,i] <- sum(z$spec[ind5])/sumsp
   }
   ## original normalizations from Martino/Goebel are invalid
   ## replaced by normalizations that are monotone and don't produce NaN's
   icafp[1,] <- log(icafp[1,]-min(icafp[1,])+1)
   icafp[1,] <- icafp[1,]/max(icafp[1,])
   icafp[2,] <- log(icafp[2,]-min(icafp[2,])+1)
   icafp[2,] <- icafp[2,]/max(icafp[2,])
   icafp[3,] <- log(icafp[3,]+1)
   icafp[3,] <- icafp[3,]/max(icafp[3,])
   icafp[4,] <- icafp[4,]/max(icafp[4,])
   icafp[5,] <- icafp[5,]/max(icafp[5,])
   icafp[6,] <- log(icafp[6,]+1)
   icafp[6,] <- icafp[6,]/max(icafp[6,])
   icafp[7:11,] <- sweep(icafp[7:11,],2,apply(icafp[7:11,],2,max),"/")
   icafp[7:11,] <- sweep(icafp[7:11,],1,apply(icafp[7:11,],1,max),"/")
   if(plot) stars(t(icafp),key.xpd=NA)
   icaobj$fingerprint <- t(icafp)
   invisible(icaobj)
}

plot.fmriICA <- function(x,comp=1,center=NULL,thresh=1.5,...){
#
#  in Anlehnung an Martini et al PNAS 2007
#
  mask = x$mask
  ddim <- dim(x$scomp)
   nt <- dim(x$A)[2]
   ncomp <- ddim[4]
   # use neurological view convention
   scomp <- x$scomp[,,,comp]
   if(is.null(center)){
   # select orthographic view with maximum information
      z <- abs(scomp)
      z[z<thresh] <- 0
      center <- rep(0,3)
      zs <- apply(z,1,sum)
      center[1] <- (1:ddim[1])[zs==max(zs)][1]
      zs <- apply(z,2,sum)
      center[2] <- (1:ddim[2])[zs==max(zs)][1]
      zs <- apply(z,3,sum)
      center[3] <- (1:ddim[3])[zs==max(zs)][1]
   }
   if(!mask[center[1],center[2],center[3]]) stop("Please specify center within brain mask\n")
   if(is.null(x$fingerprint)) x <- ICAfingerprint(x)
   n1 <- ddim[1]
   n2 <- ddim[2]
   n3 <- ddim[3]
   scomp[!x$mask] <- 0
   indx <- (1:n1)[apply(x$mask,1,any)]
   indy <- (1:n2)[apply(x$mask,2,any)]
   indz <- (1:n3)[apply(x$mask,3,any)]
   scomp <- scomp[indx,indy,indz]
   n1 <- length(indx)
   n2 <- length(indy)
   n3 <- length(indz)
   center[1] <- (1:n1)[indx==center[1]]
   center[2] <- (1:n2)[indy==center[2]]
   center[3] <- (1:n3)[indz==center[3]]
   n23 <- max(n2,n3)
   wh <- 2*n1+n2+n2/8 + n23/1.5
   wv <- n23 + n23/2
   mat <- matrix(c(1,1,7,
                   2,2,7,
                   3,3,8,
                   4,5,8,
                   6,6,8),3,5)
   layout(mat,widths=c(n2,n1,n1,n2/8,n23/1.5)/wh,
              heights=c(n23/2,n23/2,n23/2)/wv)
   oldpar <- par(mar=c(2.5,2.5,2.5,.1),mgp=c(1.5,.75,0))
   scompp <- scompn <- scomp
   scompp[scomp < thresh] <- NA
   scompn[scomp > -thresh] <- NA
   rs <- range(scomp,-thresh,thresh)
   rsp <- c(thresh,rs[2])
   rsn <- c(rs[1],-thresh)
   rs <- c(-thresh,thresh)
   indx <- indx*x$voxdim[1]
   indy <- indy*x$voxdim[2]
   indz <- indz*x$voxdim[3]
   image(-indy[n2:1],indz,scomp[center[1],n2:1,],zlim=rs,col=grey(0:255/255),asp=TRUE,xlab="-y (mm)", ylab="z (mm)")
   title(paste("Component",comp,"sagittal"))
   lines(-indy[c(1,n2)],rep(indz[center[3]],2),col=2)
   lines(rep(-indy[center[2]],2),indz[c(1,n3)],col=2)
   image(-indy[n2:1],indz,scompp[center[1],n2:1,],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
   image(-indy[n2:1],indz,scompn[center[1],n2:1,],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
   image(indx,indz,scomp[,center[2],],zlim=rs,col=grey(0:255/255),asp=TRUE, xlab="x (mm)", ylab="z (mm)")
   title(paste("Component",comp,"coronal"))
   lines(indx[c(1,n1)],rep(indz[center[3]],2),col=2)
   lines(rep(indx[center[1]],2),indz[c(1,n3)],col=2)
   image(indx,indz,scompp[,center[2],],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
   image(indx,indz,scompn[,center[2],],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
   image(indx,indy,scomp[,,center[3]],zlim=rs,col=grey(0:255/255),asp=TRUE, xlab="x (mm)", ylab="y (mm)")
   title(paste("Component",comp,"axial"))
   lines(indx[c(1,n1)],rep(indy[center[2]],2),col=2)
   lines(rep(indx[center[1]],2),indy[c(1,n2)],col=2)
   image(indx,indy,scompp[,,center[3]],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
   image(indx,indy,scompn[,,center[3]],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
   scalep <- seq(thresh,max(scomp,thresh),length=256)
   scalen <- seq(min(scomp,-thresh),-thresh,length=256)
   scalep <- t(matrix(scalep,length(scalep),10))
   scalen <- t(matrix(scalen,length(scalen),10))
   image(1:10,scalep[1,],scalep,col=heat.colors(256),xaxt="n",xlab="",ylab="signal")
   image(1:10,scalen[1,],scalen,col=rainbow(256,start=.4,end=.7)[256:1],xaxt="n",xlab="",ylab="signal")
   stars(rbind(x$fingerprint[comp,],rep(0,11)),ncol=1,scale=FALSE,
   labels=c(paste("comp",comp),""),key.loc=c(2.3,2.2))
   title("IC fingerprint")
   plot((1:nt)*x$TR,x$A[comp,],xlab="time(s)",ylab="signal",type="l",main="Time series")
   cspectr <- spectrum(x$A[comp,],plot=FALSE)
   plot(cspectr$freq/x$TR,cspectr$spec,xlab="frequency(Hz)",ylab="spectral density",type="l",main="Spectral density")
   par(oldpar)
   invisible(NULL)
}

fmri.sgroupICA <- function(icaobjlist,thresh=.75,minsize=2){
   nobj <- length(icaobjlist)
   ddim <- dim(icaobjlist[[1]]$scomp)
   voxdim <- icaobjlist[[1]]$voxdim
   ncomp <- 0
   for(i in 1:nobj){
      ncomp <- ncomp + dim(icaobjlist[[i]]$A)[1]
   }
   scompall <- array(0,c(ddim[1:3],ncomp))
   label <- rep(0,ncomp)
   last <- 0
   for(i in 1:nobj){
      nci <- dim(icaobjlist[[i]]$A)[1]
      label[last+1:nci] <- i
      scompall[,,,last+1:nci] <- icaobjlist[[i]]$scomp
      last <- last+nci
   }
   dim(scompall) <- c(prod(ddim[1:3]),ncomp)
   scompall <- scompall[icaobjlist[[1]]$mask,]
   CCs <- cor(scompall)
   SM <- abs(CCs)
   dim(SM) <- dim(CCs)
   DM <- sqrt(1-SM)
   dim(DM) <- dim(CCs)
   DM[outer(label,label,"==")] <- 1
   DM <- as.dist(DM)
# we will be using hclust with complete linkage so that components
# from same runs / subjects will not be in the same cluster
   hdm <- hclust(DM)
   nsteps <- sum(hdm$height<thresh)
   cluster <- -(1:ncomp)
   height <- rep(0,ncomp)
   for(i in 1:nsteps){
     height[cluster==hdm$merge[i,1]] <- hdm$height[i]
     height[cluster==hdm$merge[i,2]] <- hdm$height[i]
      cluster[cluster==hdm$merge[i,1]] <- i
      cluster[cluster==hdm$merge[i,2]] <- i
   }
   cl <- sort(unique(cluster))
   cl <- cl[cl>0]
   icacomp <- array(0,c(prod(ddim[1:3]),length(cl)))
   size <- numeric(length(cl))
   for(i in 1:length(cl)){
       ind <- (1:ncomp)[cluster==cl[i]]
       csign <- sign(CCs[ind[1],ind])
       sizei <- length(ind)
       icacomp[icaobjlist[[1]]$mask,i] <- scompall[,ind]%*%csign/sizei
       size[i] <- sizei
   }
   ind <- size>=minsize
   cl <- cl[ind]
   size <- size[ind]
   icacomp <- icacomp[,ind]
   dim(icacomp) <- c(ddim[1:3],length(cl))
   z<-list(icacomp=icacomp, cluster=cluster, height=height, hdm=hdm, cl=cl,
           size=size, voxdim=voxdim)
   class(z) <- "fmrigroupICA"
   invisible(z)
}

plot.fmrigroupICA <- function(x,comp=1,center=NULL,thresh=1.5,...){
ddim <- dim(x$icacomp)
oind <- order(x$size,decreasing=TRUE)
size <- x$size[oind[comp]]
scomp <- x$icacomp[,,,oind[comp]]
if(is.null(center)){
# select orthographic view with maximum information
   z <- abs(scomp)
   z[z<thresh] <- 0
   center <- rep(0,3)
   zs <- apply(z,1,sum)
   center[1] <- (1:ddim[1])[zs==max(zs)][1]
   zs <- apply(z,2,sum)
   center[2] <- (1:ddim[2])[zs==max(zs)][1]
   zs <- apply(z,3,sum)
   center[3] <- (1:ddim[3])[zs==max(zs)][1]
}
mask <- scomp!=0
if(!mask[center[1],center[2],center[3]]) stop("Please specify center within brain mask\n")
n1 <- ddim[1]
n2 <- ddim[2]
n3 <- ddim[3]
indx <- (1:n1)[apply(mask,1,any)]
indy <- (1:n2)[apply(mask,2,any)]
indz <- (1:n3)[apply(mask,3,any)]
scomp <- scomp[indx,indy,indz]
n1 <- length(indx)
n2 <- length(indy)
n3 <- length(indz)
center[1] <- (1:n1)[indx==center[1]]
center[2] <- (1:n2)[indy==center[2]]
center[3] <- (1:n3)[indz==center[3]]
n23 <- max(n2,n3)
wh <- 2*n1+n2+n2/8
mat <- matrix(c(1,1,
                2,2,
                3,3,
                4,5),2,4)
layout(mat,widths=c(n2,n1,n1,n2/8,n23)/wh,
           heights=c(1/2,1/2))
oldpar <- par(mar=c(2.5,2.5,2.5,.1),mgp=c(1.5,.75,0))
scompp <- scompn <- scomp
scompp[scomp < thresh] <- NA
scompn[scomp > -thresh] <- NA
rs <- range(scomp,-thresh,thresh)
rsp <- c(thresh,rs[2])
rsn <- c(rs[1],-thresh)
rs <- c(-thresh,thresh)
indx <- indx*x$voxdim[1]
indy <- indy*x$voxdim[2]
indz <- indz*x$voxdim[3]
image(-indy[n2:1],indz,scomp[center[1],n2:1,],zlim=rs,col=grey(0:255/255),asp=TRUE,xlab="-y (mm)",ylab="z (mm)")
title(paste("Component",comp,"size",size,"sagittal"))
lines(-indy[c(1,n2)],rep(indz[center[3]],2),col=2)
lines(rep(-indy[center[2]],2),indz[c(1,n3)],col=2)
image(-indy[n2:1],indz,scompp[center[1],n2:1,],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
image(-indy[n2:1],indz,scompn[center[1],n2:1,],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
image(indx,indz,scomp[,center[2],],zlim=rs,col=grey(0:255/255),asp=TRUE,xlab="x (mm)",ylab="z (mm)")
title(paste("Component",comp,"size",size,"coronal"))
lines(indx[c(1,n1)],rep(indz[center[3]],2),col=2)
lines(rep(indx[center[1]],2),indz[c(1,n3)],col=2)
image(indx,indz,scompp[,center[2],],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
image(indx,indz,scompn[,center[2],],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
image(indx,indy,scomp[,,center[3]],zlim=rs,col=grey(0:255/255),asp=TRUE,xlab="x (mm)",ylab="y                (mm)")
title(paste("Component",comp,"size",size,"axial"))
lines(indx[c(1,n1)],rep(indy[center[2]],2),col=2)
lines(rep(indx[center[1]],2),indy[c(1,n2)],col=2)
image(indx,indy,scompp[,,center[3]],zlim=rsp,add=TRUE,col=heat.colors(256),asp=TRUE)
image(indx,indy,scompn[,,center[3]],zlim=rsn,add=TRUE,col=rainbow(256,start=.4,end=.7)[256:1],asp=TRUE)
scalep <- seq(thresh,max(scomp,thresh),length=256)
scalen <- seq(min(scomp,-thresh),-thresh,length=256)
scalep <- t(matrix(scalep,length(scalep),10))
scalen <- t(matrix(scalen,length(scalen),10))
image(1:10,scalep[1,],scalep,col=heat.colors(256),xaxt="n",xlab="",ylab="signal")
image(1:10,scalen[1,],scalen,col=rainbow(256,start=.4,end=.7)[256:1],xaxt="n",xlab="",ylab="signal")
par(oldpar)
invisible(NULL)
}
