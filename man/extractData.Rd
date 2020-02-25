\name{extractData}
\alias{extractData}
\alias{condensefMRI}
\alias{expandfMRI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract data or residuals from a fmridata object }
\description{
  The function extracts data stored as raw within an object of class
'fmridata'.
}
\usage{
extractData(z, what = "data", maskOnly = FALSE)
expandfMRI(z)
condensefMRI(z, mask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ an object of class 'fmridata' }
  \item{what}{ either \code{"data"} or \code{"residuals"}. }
  \item{maskOnly}{logical: if TRUE only values within the brain mask will be returned.}
  \item{mask}{logical brain mask}
}
\details{
  The function \code{extractData} extracts data stored as raw within an object of class
'fmridata'. Functions \code{expandfMRI} and \code{condensefMRI} change the way
data and residuals are stored between full 3D data and data within a brain mask.
\code{condensefMRI} can also be used to set a more restrictive brain mask.
}
\value{
  In case of function \code{extractData} an array of dimension \code{data$dim} containing either the
  fmri-data or residuals. The other two functions return an object of class 'fmridata'.
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{\code{\link{fmri.lm}} }
\keyword{ utilities }% at least one, from doc/KEYWORDS
