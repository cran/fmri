\name{fmri.lm}
\alias{fmri.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for fMRI data }
\description{
  Estimate the parameters and variances in a linear model. 
}
\usage{
  fmri.lm(ds, z, mask = NULL,
          actype = c("smooth", "noac", "ac", "accalc"),
          contrast = c(1), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{ Data object of class "fmridata" }
  \item{z}{ Designmatrix specifying the expected BOLD response(s) 
            and additional components for trend and other effects. }
  \item{mask}{ Array of dimensionality of the data describing a (brain) mask 
               the computation should be restricted to. The default is the 
               mask given with the data. }
  \item{actype}{ String describing the type of handling autocorrelation
                 of time series. One of "smooth", "nonac", "ac", "accalc".}
  \item{contrast}{ Contrast vector for the covariates.}
  \item{verbose}{ Verbose mode, default is \code{FALSE}. }
}
\details{
  This function performs parameter estimation in the linear model.
  It implements a two step procedure. After primary estimation of the
  parameters in the first step residuals are obtained. If \code{actype} 
  \code{\%in\%} \code{c("ac", "accalc", "smooth")} an AR(1) model is 
  fitted, in each voxel, to the time series of residuals. The estimated 
  AR-coefficients are corrected for bias. If \code{actype=="smooth"}
  the estimated AR-coefficients are spatially smoothed. If \code{actype} 
  \code{\%in\%} \code{c("ac", "smooth")} the linear model is pre-withened
  using the estimated (and possibly smoothed) AR-coefficients. Parameter 
  and variance estimates are then obtained from the pre-withened
  data. The argument \code{keep} describes the amount of data which is
  returned. The estimated effects 
  \deqn{\tilde{\gamma}_i = C^T\tilde{\beta}_i}
  and their estimated variances are returned as well as the 
  residuals and temporal autocorrelation.
  \code{cbeta} then contains the corresponding parameter
  estimates and thus is a vector of corresponding length in each voxel.

  If warning "Local smoothness characterized by large bandwidth" occurs,
  check \code{scorr} elements. If correlation drops with lag towards
  zero, data has been pre-smoothed. Adaptive smoothing the SPM can then 
  only be of limited use. If correlation does not go to zero, check the 
  residuals of the linear model for unexplained structure (spin saturation 
  in first scans? discard them!).
}
\value{
  object with class attributes "fmrispm" and "fmridata"
  \item{beta}{estimated parameters }
  \item{cbeta}{estimated contrast of parameters}
  \item{var}{estimated variance of the contrast of parameters.}
  \item{varm}{covariance matrix of the parameters given by \code{vvector}}
  \item{res}{raw (integer size 2) vector containing residuals of the estimated linear model
             up to scale factor resscale.}
  \item{resscale}{\code{resscale*extract.data(object,"residuals")} are the residuals. }
  \item{dim}{dimension of the data cube and residuals}
  \item{arfactor}{estimated autocorrelation parameter}
  \item{rxyz}{array of smoothness from estimated correlation for each
    voxel in resel space (for analysis without smoothing)}
  \item{scorr}{array of spatial correlations with maximal lags 5, 5, 3 in x,y and z-direction. }
  \item{bw}{vector of bandwidths (in FWHM) corresponding to the spatial correlation within the data.}
  \item{weights}{ratio of voxel dimensions }
  \item{vwghts}{ratio of estimated variances for the stimululi given by
  \code{vvector}}
  \item{mask}{head mask.}
  \item{df}{Degrees of freedom for t-statistics.}
  \item{hrf}{expected BOLD response for contrast }
}
\references{ Worsley, K.J. (2005). Spatial smoothing of autocorrelations
  to control the degrees of freedom in fMRI analysis. NeuroImage,
  26:635-641.

  Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1-15.

  Tabelow, K., Polzehl, J., Voss, H.U., and Spokoiny, V. (2006). 
\emph{Analysing {fMRI}
    experiments with structure adaptive smoothing procedures},
  NeuroImage, 33:55-62. 

}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ The argument \code{vvector} is no longer supported. 
}
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.stimulus}} }
\examples{
  \dontrun{
  # Example 1
  data <- list(ttt=writeBin(rnorm(32*32*32*107), raw(), 4),
               mask=array(1, c(32, 32, 32)), dim=c(32, 32, 32, 107))
  class(data) <- "fmridata"
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  model <- fmri.lm(data, z, verbose=TRUE)
  plot(extract.data(data)[16, 16, 16,])
  lines(extract.data(data)[16, 16, 16, ] - extract.data(model, "residuals")[16, 16, 16, ], col=2)
  }
}
\keyword{ regression }% at least one, from doc/KEYWORDS

