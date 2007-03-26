\name{fmri.lm}
\alias{fmri.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model For FMRI Data }
\description{
  Estimate the parameters and variances in a linear model. 
}
\usage{
  fmri.lm(data, z, actype = "accalc", hmax = 3.52, vtype = "var",
          step = 0.01, contrast = c(1), vvector = c(1),
          keep = "essential")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ object of class "fmridata" }
  \item{z}{ designmatrix specifying the expected BOLD response(s) 
  and additional components for trend and other effects. }
  \item{actype}{ string describing the type of handling autocorrelation
    of time series. "nonac", "ac", "accalc", "smooth"}
  \item{hmax}{ bandwidth for smoothing autocorrelation parameter if
    actype = "smooth" }
  \item{vtype}{ method of estimating residual variance (only
  "var" implemented) }
  \item{step}{ step size for binning autocorrelations (see
    details) }
  \item{contrast}{ contrast vector }
  \item{vvector}{ vector defining the parameters for which the
    covariance matrix is returned as well as the corresponding length of
    the vector \code{cbeta} in each voxel}
  \item{keep}{ string describing the amount of data returned:
    "essential", "diagnostic", "all" }
}
\details{
  This function performs parameter estimation in the linear model.
  It implements a two step procedure. After primary estimation of the
  parameters in the first step residuals
  are obtained. If \code{actype} \code{\%in\%} \code{c("ac", "accalc", "smooth")} 
   an AR(1) model is fitted, in each voxel, to 
  the time series of residuals. The estimated AR-coefficient is corrected for bias. 
  If \code{actype=="smooth"}
  the estimated AR-coefficients are spatially smoothed using bandwidth \code{hmax}.
  If \code{actype} \code{\%in\%} \code{c("ac", "smooth")} the linear model is prewithened
  using the estimated (smoothed) AR-coefficients. Parameter 
  and variance estimates are then obtained from the prewithened
  data. The argument \code{keep} describes the amount of data which is
  returned. If "essential" only the estimated effects 
  \deqn{\tilde{\gamma}_i = C^T\tilde{\beta}_i}
   and their
  estimated variances are returned. "all" gives the full data, including
  residuals, temporal autocorrelation.
  If \code{vvector} is given and has length greater than 1, the
  covariance matrix for the stimuli given therein are returned
  (\code{varm}) and \code{vwghts} contains an estimate for the ratio of
  the variances of the parameter for the stimuli indicated in
  \code{vvector}. \code{cbeta} then contains the corresponding parameter
  estimates and thus is a vector of corresponding length in each voxel.
}
\value{
  object with class attributes "fmrispm" and "fmridata"
  \item{beta}{estimated parameters }
  \item{cbeta}{estimated contrast of parameters}
  \item{var}{estimated variance of the contrast of parameters.}
  \item{varm}{covariance matrix of the parameters given by \code{vvector}}
  \item{res}{residuals of the estimated linear model}
  \item{arfactor}{estimated autocorrelation parameter}
  \item{scorr}{spatial correlation of data }
  \item{weights}{ratio of voxel dimensions }
  \item{vwghts}{ratio of estimated variances for the stimululi given by
  \code{vvector}}
  \item{rxyz}{array of smoothness from estimated correlation for each
    voxel in resel space (for analysis without smoothing)}
  \item{hrf}{expected BOLD response for contrast }
}
\references{ Worsley, K.J. (2005). Spatial smoothing of autocorrelations
  to control the degrees of freedom in fMRI analysis. NeuroImage,
  26:635-641.

  Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1-15.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ \code{vvector} is intended to be used for delay of the HRF using
  its first derivative. Do not mix with the \code{contrast} argument, since
  unexpected side effects may occur. Look out for updates of this package.
}
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.stimulus}} }
\examples{
  # Example 1
  data <- list(ttt=array(rnorm(32*32*32*107),c(32,32,32,107)),
               mask=array(1,c(32,32,32)),dim=c(32,32,32))
  class(data) <- "fmridata"
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  model <- fmri.lm(data,z,keep="all")
  plot(data$ttt[16,16,16,])
  lines(data$ttt[16,16,16,] - model$res[16,16,16,],col=2)
}
\keyword{ regression }% at least one, from doc/KEYWORDS

