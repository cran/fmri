\name{fmri.smooth}
\alias{fmri.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing Statistical Parametric Maps }
\description{
  Perform the adaptive weights smoothing procedure 
}
\usage{
fmri.smooth(spm, hmax = 4, adaptive=TRUE,
            lkern="Gaussian", skern="Plateau", weighted=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ object of class \code{fmrispm} }
  \item{hmax}{ maximum bandwidth to smooth }
  \item{adaptive}{ logical. TRUE (default) for adaptive smoothing }
  \item{lkern}{ \code{lkern} specifies the location kernel. Defaults to
    "Gaussian", other choices are "Triangle" and "Plateau". Note that the location kernel is applied to
    \code{(x-x_j)^2/h^2}, i.e. the use of "Triangle" corresponds to the
    Epanechnicov kernel in nonparametric kernel regression. "Plateau" specifies a kernel that is equal to 1 in the interval (0,.3),
    decays linearly in (.5,1) and is 0 for arguments larger than 1.
    }
  \item{skern}{ \code{skern} specifies the kernel for the statistical
    penalty. Defaults to "Plateau", the alternatives are "Triangle" and "Exp".
    "Plateau" specifies a kernel that is equal to 1 in the interval (0,.3),
    decays linearly in (.3,1) and is 0 for arguments larger than 1.
    \code{lkern="Plateau"} and \code{lkern="Triangle"} allow for much faster computation (saves up
    to 50\% CPU-time). \code{lkern="Plateau"} produces a less random weighting scheme.}
  \item{weighted}{ \code{weighted} (logical) determines if weights contain the inverse of local 
   variances as a factor (Weighted Least Squares). \code{weighted=FALSE} does not employ the
   heteroscedasticity of variances for the weighting scheme and is preferable if variance estimates
   are highly variable, e.g. for short time series.} 
}
\details{
  This function performs the smoothing on the Statistical Parametric Map
  spm.

  \code{hmax} is the (maximal) bandwidth used in the last iteration. Choose
  \code{adaptive} as \code{FALSE} for non adaptive
  smoothing. \code{lkern} can be used for specifying the
  localization kernel. For comparison with non adaptive methods use
  "Gaussian" (hmax given in FWHM), for better adaptation use "Plateau" or "Triangle"
  (default, hmax given in voxel). For \code{lkern="Plateau"} and \code{lkern="Triangle"} thresholds may be inaccurate, due to a violation of
    the Gaussian random field assumption under homogeneity. \code{lkern="Plateau"} 
    is expected to provide best results with adaptive smoothing.
   
  \code{skern} can be used for specifying the
  kernel for the statistical penalty. "Plateau" is expected to provide the best results,
  due to a less random weighting scheme.

  The function handles zero variances by assigning a large value (1e20)
  to these variances. Smoothing is restricted to voxel with \code{spm$mask}.
}
\value{
  object with class attributes "fmrispm" and "fmridata"
  \item{cbeta}{ smoothed parameter estimate}
  \item{var}{ variance of the parameter}
  \item{hmax}{ maximum bandwidth used}
  \item{rxyz}{ smoothness in resel space. all directions }
  \item{rxyz0}{ smoothness in resel space as would be achieved by
    a Gaussian filter with the same bandwidth. all directions }
  \item{scorr}{array of spatial correlations with maximal lags 5, 5, 3 in x,y and z-direction. }
  \item{bw}{vector of bandwidths (in FWHM) corresponding to the spatial correlation within the data.}
  \item{dim}{dimension of the data cube and residuals}
  \item{weights}{ ratio of voxel dimensions }
  \item{vwghts}{ ratio of estimated variances for the stimuli given by
  \code{vvector}}
  \item{hrf}{ Expected BOLD response for the specified effect }
}
\references{
  Tabelow, K., Polzehl, J., Voss, H.U., and Spokoiny, V.. 
\emph{Analysing {fMRI}
    experiments with structure adaptive smoothing procedures},
  NeuroImage, 33:55-62 (2006). 

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
