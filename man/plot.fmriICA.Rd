\name{plot.fmriICA}
\alias{plot.fmriICA}
\alias{plot.fmrigroupICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostics plots for objects of class ''\code{fmriICA}''
}
\description{
The function generates plots for inspecting independent components
obtained by spatial indipendent component analysis.
}
\usage{
\method{plot}{fmriICA}(x, comp = 1, center = NULL, thresh = 1.5, ...)
\method{plot}{fmrigroupICA}(x, comp = 1, center = NULL, thresh = 1.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   object returned by function \code{fmri.sICA} or preferably
   function \code{ICAfingerprinting} in case of \code{plot.fmriICA}
   and object returned by function \code{fmri.sgroupICA} in case of
    \code{plot.fmrigroupICA}
}
  \item{comp}{
  number of the independent component to inspect.
}
  \item{center}{
  coordinates for central point to determine axial, coronal and sattigal
  slices for display. If NULL the central point of the image cube is selected.
  center needs to be within the brain mask.
}
  \item{thresh}{
  Threshold value
}
\item{...}{currently not used}
}
\details{ The function generates diagnostic plots for the independent component specified
in \code{comp}. It
 provides axial, coronal and sattigal images as determined by \code{center}.
 Values exceeding the threshold are displayed using a color scale.
 An IC fingerprint is given as a star plot.
 Additionally the time series corresponding to the spatial IC and its spectral density are plottet.
}
\value{ nothing returned.
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
\code{\link{fmri.sICA}}, \code{\link{ICAfingerprint}}, \code{\link[fastICA]{fastICA}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
\keyword{ utilities }% __ONLY ONE__ keyword per line
