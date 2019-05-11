\name{setmask}
\alias{setmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add or replace mask in an fmridata object
}
\description{
The function replaces the information in the mask component of an fmridata object.
}
\usage{
setmask(fmriobj, mask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmriobj}{
  object of class 'fmridata'
}
  \item{mask}{
  object of class 'array' or 'nifti'
}
}
\details{
   The dimensions of both objects supplied as arguments need to be compatible.
}
\value{
on object of class 'fmridata'.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
 \code{\link{oro2fmri}}, \code{\link{niftiImage2fmri}}, \code{\link{read.NIFTI}},
 \code{\link{read.AFNI}}, \code{\link{read.ANALYZE}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% use one of  RShowDoc("KEYWORDS")
\keyword{ utilities }% __ONLY ONE__ keyword per line
