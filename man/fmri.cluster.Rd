\name{fmri.cluster}
\alias{fmri.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Cluster thresholding.
}
\description{
Detection of activated regions using cluster thresholding.

}
\usage{
fmri.cluster(spm, alpha = 0.05, ncmin = 2, minimum.signal = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{\code{fmrispm} object
}
  \item{alpha}{
   multiple test (over volume and cluster sizes) adjusted significance level
   used for thresholds.
}
  \item{ncmin}{
   minimal cluster size used. An activation is detected if for any
   clustersize in \code{nvmin:20} the size specific threshold is exceeded.
}
  \item{minimum.signal}{
  allows to specify a (positive) minimum value for detected signals. If minimum.signal >0
   the thresholds are to conservative, this case needs further improvements.
}
}
\details{ Approximate thresholds for the existence of a cluster with mean t-statistics
exceeding a \code{1-beta} threshold \code{k_{nc}} for cluster size \code{nc}
are based on a simulation study under the hypothesis and adjusted for both degrees
of freedom, number of voxel in mask and spatial correlation.
\code{beta} is choosen such that under the hypothesis the combined (over cluster sizes
\code{ncmin:20}) test has approximate significance level \code{alpha}.
}
\value{
Object with class attributes "fmripvalue" and "fmridata"
\item{pvalue}{ p-values of voxelwise t-statistics for voxel that were detected
   for any cluster size, a value of \code{1} otherwise.}
\item{weights}{ voxelsize ratio }
\item{dim}{ data dimension }
\item{hrf}{ expected BOLD response for contrast (single stimulus only)}
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmri.lm}}, \code{\link{fmri.pvalue}}, \code{\link{fmri.searchlight}}
}
\examples{
  \dontrun{fmri.cluster(fmrispmobj)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
