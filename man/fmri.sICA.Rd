\name{fmri.sICA}
\alias{fmri.sICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spacial ICA for fmri data}
\description{ Uses fastICA to perform spatial ICA on fMRI data.}
\usage{
fmri.sICA(data, mask=NULL, ncomp=20,
  alg.typ=c("parallel","deflation"), fun=c("logcosh","exp"),
  alpha=1, detrend=TRUE, degree=2, nuisance= NULL, ssmooth=TRUE,
  tsmooth=TRUE, bwt=4, bws=8, unit=c("FWHM","SD"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ fMRI dataset of class ''\code{fmridata}'' }
  \item{mask}{ Brain mask, if \code{NULL} then \code{data$mask} is used. }
  \item{ncomp}{Number of ICA components to compute.}
  \item{alg.typ}{Alg. to be used in \code{fastICA}.}
  \item{fun}{Test functions to be used in \code{fastICA}.}
  \item{alpha}{Scale parameter in test functions, see \code{fastICA}.}
  \item{detrend}{Trend removal (polynomial)}
  \item{degree}{degree of polynomial trend}
  \item{nuisance}{ Matrix of additional nuisance parameters to regress against. }
  \item{ssmooth}{ Should spatial smoothing be used for variance reduction}
  \item{tsmooth}{ Should temporal smoothing be be applied}
  \item{bws}{Bandwidth for spatial Gaussian kernel}
  \item{bwt}{Bandwidth for temporal Gaussian kernel}
  \item{unit}{Unit of bandwidth, either standard deviation (SD) of Full Width Half Maximum (FWHM).}
}
\details{
If specified polynomial trends and effects due to nuisance parameters, e.g.,
motion parameters, are removed. If \code{smooth==TRUE} the resulting residual series is
spatially smoothed using a Gaussian kernel with specified bandwidth.
ICA components are the estimated using fastICA based on data within brain mask.
The components of the result are related as \code{XKW=scomp[mask,]} and \code{X=scomp[mask,]*A}.
}
\value{object of class ''\code{fmriICA}''
  list with components
  \item{scomp}{4D array with ICA component images. Last index varies over components. }
  \item{X}{pre-processed data matrix}
  \item{K}{pre-processed data matrix}
  \item{W}{estimated un-mixing matrix}
  \item{A}{estimated mixing matrix}
  \item{mask}{Brain mask}
  \item{pixdim}{voxelsize}
  \item{TR}{Repetition Time (TR)}
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{\code{\link{plot.fmriICA}},\code{\link{ICAfingerprint}}, \code{\link[fastICA]{fastICA}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "multivariate" }
