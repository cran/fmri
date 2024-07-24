#########################################################################################################################
#
#    R - function  aws3D  for vector-valued  Adaptive Weights Smoothing (AWS) in 3D
#
#    reduced version for qtau==1, heteroskedastic variances only, exact value for Variance reduction
#
#    emaphazises on the propagation-separation approach
#
#    Copyright (C) 2005 Weierstrass-Institut fuer
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

aws3D <- function(y, qlambda=NULL, lkern="Gaussian", skern="Plateau", weighted=TRUE,
                  sigma2=NULL, mask=NULL, hmax=NULL, ladjust=1, u=NULL, wghts=NULL,
                  h0=c(0,0,0), testprop=FALSE, res=NULL){
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  uses sum of weights and correction factor C(h,g) in statistical penalty
#
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,Plateau=.9945,
		                                            Triangle=.9988,
																								Exponential=.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- ladjust*qchisq(qlambda,1)
  } else {
    lambda <- 1e50
  }

  z <- aws::aws3Dmask(y, mask, lambda, hmax, res, sigma2, lkern, skern,
	                    weighted, u, wghts, h0, testprop)

  ## "compress" the residuals
  scale <- max(abs(range(z$residuals)))/32767
  z$residuals <- writeBin(as.integer(z$residuals/scale), raw(), 2)
	z$resscale <- scale
	z$maskOnly <- TRUE
  #   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
  class(z) <- "aws.gaussian"
  invisible(z)
}

aws3Dfull <- function(y, qlambda=NULL, lkern="Gaussian", skern="Plateau", weighted=TRUE,
                   sigma2=NULL, mask=NULL, hmax=NULL, ladjust=1, u=NULL, wghts=NULL,
                   testprop=FALSE, res=NULL, resscale=NULL) {
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  implements estimation of true variances from residuals (slower) than aws3D
#
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
    # MAE
  mae <- NULL


  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,
		            Plateau=.9945,
								Triangle=.9988,
								Exponential=.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- ladjust*qchisq(qlambda,1)
  } else {
    lambda <- 1e50
  }

	z <- aws::aws3Dmaskfull(y, mask, lambda, hmax, res, sigma2, lkern, skern,
											weighted, u, wghts, testprop)

	## "compress" the residuals
	scale <- max(abs(range(z$residuals)))/32767
	z$residuals <- writeBin(as.integer(z$residuals/scale), raw(), 2)
	z$resscale <- scale
	z$maskOnly <- TRUE
	#   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
	class(z) <- "aws.gaussian"
	invisible(z)
}
