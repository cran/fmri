\name{fmri.smooth}
\alias{fmri.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing Statistical Parametric Maps }
\description{
  Perform the adaptive weights smoothing procedure 
}
\usage{
fmri.smooth(spm, hmax = 4, adaptive=TRUE,
            lkern="Triangle", skern="Triangle", na.rm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ object of class \code{fmrispm} }
  \item{hmax}{ maximum bandwidth to smooth }
  \item{adaptive}{ logical. TRUE (default) for adaptive smoothing }
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to
    "Triangle", other choices are "Gaussian", "Quadratic", "Cubic" and
    "Uniform". Note that the location kernel is applied to
    \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the
    Epanechnicov kernel in nonparametric kernel regression. }
  \item{skern}{ \code{skern} specifies the kernel for the statistical
    penalty. Defaults to "Triangle", the alternative is "Exp".
    \code{lkern="Triangle"} allows for much faster computation (saves up
    to 50\%).}
  \item{na.rm}{\code{na.rm} specifies how NA's in the SPM are handled. NA's may occur
   in voxel where the time series information did not allow for estimating parameters and their variances
   or where the time series information where constant over time. A high (1e19) value of the variance 
   and a parameter of 0 are used to characterize NA's. If \code{na.rm=FALSE} these values are simply downweigthed
   by their high variance and estimates in these voxels are produced as nonadaptive averages
   of values from neighboring voxels. If \code{na.rm=TRUE} no estimates are computed (the values 
   for mean and variance are kept) in these voxels.}
}
\details{
  This function performs the smoothing on the Statistical Parametric Map
  spm.

  \code{hmax} is the (maximal) bandwidth used in the last iteration. Choose
  \code{adaptive} as \code{FALSE} for non adaptive
  smoothing. \code{lkern} can be used for specifying the
  localization kernel. For comparison with non adaptive methods use
  "Gaussian" (hmax given in FWHM), for better adaptation use "Triangle"
  (default, hmax given in voxel). \code{skern} can be used for specifying the
  kernel for the statistical penalty.  

  The function handles zero variances by assigning a large value (1e20)
  to these variances. 
}
\value{
  object with class attributes "fmrispm" and "fmridata"
  \item{cbeta}{ smoothed parameter estimate}
  \item{var}{ variance of the parameter}
  \item{hmax}{ maximum bandwidth used}
  \item{rxyz}{ smoothness in resel space. all directions }
  \item{rxyz0}{ smoothness in resel space as would be achieved by
    a Gaussian filter with the same bandwidth. all directions }
  \item{scorr}{ spatial correlation of original data }
  \item{weights}{ ratio of voxel dimensions }
  \item{vwghts}{ ratio of estimated variances for the stimuli given by
  \code{vvector}}
  \item{hrf}{ Expected BOLD response for the specified effect }
}
\references{
  Tabelow, K., Polzehl, J., Voss, H.U., and Spokoiny, V. (2005). 
  \emph{Analysing {fMRI}
    experiments with structure adaptive smoothing procedures},
  NeuroImage, accepted (2006). 

  Polzehl, J. and Spokoiny, V. (2006). \emph{Propagation-Separation
    Approach for Local Likelihood Estimation}, Probab. Theory
  Relat. Fields 135, 335-362.
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de},
  Karsten Tabelow \email{tabelow@wias-berlin.de}  }
\examples{
\dontrun{fmri.smooth(spm, hmax = 4, lkern = "Gaussian")}
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
