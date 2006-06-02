\name{fmri.pvalue}
\alias{fmri.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ P-values }
\description{
  Determines p-values.
}
\usage{
fmri.pvalue(spm, mode="basic", delta=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{ \code{fmrispm} object }
  \item{mode}{ type of pvalue definition }
  \item{delta}{ physically meaningful range of latency for HRF }
}
\details{
  If only a contrast is given in \code{spm}, we simply use a t-statistic
  and define p-values according to random field theory for the resulting gaussian
  field (sufficiently large number of df - see ref.). If \code{spm} is a
  vector of length larger then one for each voxel, a chisq field is calculated and evaluated (see
  Worsley and Taylor (2006)). If \code{delta} is given, a cone statistics is
  used.

  The parameter \code{mode} allows for different kinds of p-value
  calculation. "basic" corresponds to a global definition of the
  resel counts based on the amount of smoothness achieved by an equivalent
  Gaussian filter. The propagation condition ensures, that under the
  hypothesis
  \deqn{\hat{\Theta} = 0}
  adaptive smoothing perform like a
  non adaptive filter with the same kernel function which justifies this
  approach. "local"
  corresponds to a more conservative setting, where the p-value is
  derived from the estimated local resel counts that has been achieved by
  adaptive smoothing. In contrast to "basic", "global" takes a global
  median to adjust for the randomness of the weighting scheme generated
  by adaptive smoothing. "global" and "local" are more conservative than
  basic, that is, they generate sligthly larger p-values. }
\value{
  Object with class attributes "fmripvalue" and "fmridata"
  \item{pvalue}{ p-value. use with \code{\link{plot}} for thresholding.}
  \item{weights}{ voxelsize ratio }
  \item{dim}{ data dimension }
  \item{hrf}{ expected BOLD response for contrast (single stimulus only)}
}
\references{
  Tabelow, K., Polzehl, J. and Spokoiny, V. (2005). \emph{Analysing {fMRI}
    experiments with structure adaptive smoothing procedures},
  WIAS-Preprint 1079. 
  
  Worsley, K.J., and Taylor, J.E., \emph{Detecting fMRI activation
    allowing for unknown latency of the hemodynamic response},
  NeuroImage 29:649-654 (2006).
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{Unxepected side effects may occur if spm does not meet the
  requirements, especially if a parameter estimate vector of length greater than 2 through
  argument \code{vvector} in \code{\link{fmri.lm}} has beeen produced for every voxel.}
\seealso{ \code{\link{fmri.smooth}}, \code{\link{plot.fmridata}}  }
\examples{
\dontrun{fmri.pvalue(smoothresult)}
}
\keyword{ htest }% at least one, from doc/KEYWORDS
