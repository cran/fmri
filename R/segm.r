###################################################################################################
#
#    R - function  segm3D  for Adaptive Weights Segmentation in 3D SPM's
#
#    exact value for Variance reduction from smoothed residuals
#
#    emaphazises on the propagation-separation approach
#
#    Copyright (C) 2010-12 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#

segm3D <- function(y,weighted=TRUE,
                   sigma2=NULL,mask=NULL,hinit=NULL,hmax=NULL,
                   ladjust=1,wghts=NULL,
                   df=100,h0=c(0,0,0),residuals=NULL, resscale=NULL,
                   delta=0,alpha=.05,restricted=FALSE) {
#
#
#  Auxilary functions
      getkrval <- function(df,ladj,fov,k,alpha){
#  define threshold for segmentation
      dfinv <- 1/df
      lfov <- log(fov)
      idffov <- dfinv*log(fov)
      lalpha <- log(alpha)
      lk <- log(k)
      ladj <- ladj-1
      kv <- 1.0105241 + 3.2036997*dfinv + 0.5442412*idffov -
            0.0002686*lalpha -0.0172885*ladj +0.0009396*lk -1.7976621*dfinv*log(df) -
            0.7856579*dfinv*lalpha+0.0765149*idffov*ladj
      kv
      }
#
# first check arguments and initialize
#
   args <- match.call()
   dmask <- dim(mask)
   n1 <- dmask[1]
   n2 <- dmask[2]
   n3 <- dmask[3]
   nvoxel <- sum(mask)
   nt <- dim(residuals)[1]
   position <- array(0,dmask)
   position[mask] <- 1:nvoxel
   if (length(y)!=nvoxel) {
     stop("segm3D: y needs to have length sum(mask)")
   }
# set the code for the kernel (used in lkern) and set lambda
   lkern <- 1
   skern <- 1
# define lambda
   lambda <- ladjust*(exp(2.6-3.17*log(df)+8.4*log(log(df)))+16)
# to have similar preformance compared to skern="Exp"
   if (is.null(hinit)||hinit<1) hinit <- 1
# define hmax
   if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points
# estimate variance in the gaussian case if necessary
# deal with homoskedastic Gaussian case by extending sigma2
   if (length(sigma2)==1) sigma2<-array(sigma2,nvoxel)
   if (length(sigma2)!=nvoxel) stop("sigma2 does not have length 1 or same length as y")
#  in these points sigma2 probably contains NA's
   sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes
# deal with homoskedastic Gaussian case by extending sigma2
  cat("\nfmri.smooth: first variance estimate","\n")
  varest0 <- aws::residualVariance(residuals, mask, resscale=1, compact=TRUE)
   vq <- varest0*sigma2
   if (is.null(wghts)) wghts <- c(1,1,1)
   hinit <- hinit/wghts[1]
   hmax <- hmax/wghts[1]
   wghts <- wghts[2:3]/wghts[1]
   tobj <- list(bi= sigma2)
   theta <- y
   segm <- numeric(nvoxel)
   varest <- varest0
   maxvol <- aws::getvofh(hmax,lkern,wghts)
   fov <- sum(mask)
   kstar <- as.integer(log(maxvol)/log(1.25))
   steps <- kstar+1
   k <- 1
   hakt <- hinit
   hakt0 <- hinit
   lambda0 <- lambda
   if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
   scorr <- numeric(3)
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
#   for(i in 10:kstar) thresh <- max(thresh,tadj*getkrval(df,ladjust,fov,i,alpha))
   thresh <- getkrval(df,ladjust,fov,kstar,alpha)
#  just to ensure monotonicity of thresh with kmax, there exist a few parameter configurations
#  where the approximation formula does not ensure monotonicity
   cat("FOV",fov,"delta",delta,"thresh",thresh,"ladjust",ladjust,"lambda",lambda,"df",df,"\n")
#
#   need these values to compute variances
#
   while (k<=kstar) {
      hakt0 <- aws::gethani(1,10,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- aws::gethani(1,10,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"bandwidth",signif(hakt,3)," ")
      dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:3]
      hakt0 <- hakt
      bi0 <- tobj$bi
      tobj <- .Fortran(C_segm3d,
                       as.double(y),
                       as.double(residuals),
                       as.double(sigma2),
                       as.integer(position),
                       as.integer(weighted),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.integer(nt),
                       as.double(df),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi0),
                       thnew=double(nvoxel),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nt),#swres
                       double(nvoxel),#pvalues
                       segm=as.integer(segm),
                       as.double(delta),
                       as.double(thresh),
                       as.double(fov),
                       as.double(vq),
                       as.double(varest0),
                       varest=as.double(varest),
                       as.integer(restricted))[c("bi","thnew","hakt","segm","varest")]
      gc()
      theta <- tobj$thnew
      segm <- tobj$segm
      varest <- tobj$varest
      if (max(total) >0) {
         cat(signif(total[k],2)*100,"%                 \r",sep="")
      }
      k <- k+1
      lambda0 <- lambda
      gc()
   }
  z <- list(theta=theta,ni=tobj$bi,var=varest,y=y,segm=segm,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),
            scorr=scorr,mask=mask)
  class(z) <- "aws.gaussian"
  invisible(z)
}
