\name{ICAfingerprint}
\alias{ICAfingerprint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IC fingerprinting
}
\description{
Implements ICA fingerprinting mainly following De Martino et.al., Neuroimage 2007.
}
\usage{
ICAfingerprint(icaobj, nbin = 256, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{icaobj}{
  object returned by function \code{fmri.sICA}.
}
  \item{nbin}{
  number of bins for entropy estimation
}
  \item{plot}{
  provide results as star plots.
}
}
\details{
For some characteristics normalization of values differs from De Martino et. al..
Frequency bands are obtained from periodogram estimated instead of using Welch's method.
}
\value{
object of class ''\code{fmriICA}''
  list with components
  \item{scomp}{4D array with ICA component images. Last index varies over components. }
  \item{X}{pre-processed data matrix}
  \item{K}{pre-processed data matrix}
  \item{W}{estimated un-mixing matrix}
  \item{A}{estimated mixing matrix}
  \item{mask}{Brain mask}
  \item{pixdim}{voxelsize}
  \item{TR}{Repetition Time (TR)}
  \item{fingerprint}{matrix of IC characteristics. Colums correspond to IC's .}
}
\references{
De Martino et. al., Classification of fMRI independent components using IC-fingerprints
and support vector machine classifieres, Neuroimage 34 (2007) 177-194.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmri.sICA}}, \code{\link{plot.fmriICA}}, \code{\link[fastICA]{fastICA}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ utilities }% __ONLY ONE__ keyword per line
