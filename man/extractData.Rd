\name{extractData}
\alias{extractData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract data or residuals from a fmridata object }
\description{
  The function extracts data stored as raw within an object of class
'fmridata'.
}
\usage{
extractData(z, what = "data")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ an object of class 'fmridata' }
  \item{what}{ either \code{"data"} or \code{"residuals"}. }
}
\details{
  The function extracts data stored as raw within an object of class
'fmridata'.
}
\value{
  an array of dimension \code{data$dim} containing either the
  fmri-data or residuals.
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{\code{\link{fmri.lm}} }
\keyword{ utilities }% at least one, from doc/KEYWORDS
