\name{getSearchlightPattern}
\alias{getSearchlightPattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Extract searchlight pattern from a SPM
}
\description{
 For a provided spm object and a mask of voxel the function extracts
 the values of the parameter estimates within the searchlight region
 and for all voxel in the mask.
}
\usage{
getSearchlightPattern(spm, voxel, radius)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{
   an object of class 'fmrispm'
}
  \item{voxel}{
   a mask (logical) with dimensionality compatible to the spm
}
  \item{radius}{
   radius of the searchlight
}
}
\value{
%%  ~Describe the value returned
an array of dimension  c(nb, nsl, nvox)
with nb the number of estimated parameters in spm$beta,
nsl the number of voxel in the searchlight and
nvox the number of voxel in the mask provided as second argument
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmri.searchlight}}, \code{\link{fmri.lm}}~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ misc }% __ONLY ONE__ keyword per line
