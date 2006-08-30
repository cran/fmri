\name{read.DICOM}
\alias{read.DICOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O function }
\description{
  Read DICOM file.
}
\usage{
read.DICOM(filename,includedata = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the file }
  \item{includedata}{ logical. should data be read too? defaults to \code{TRUE}. }
}
\details{
  The function reads a DICOM file.
}
\value{
  Object with the following list entries:
  \item{header}{header information as raw data}
  \item{ttt}{ image data if requested}
  \item{format}{data source. string "DICOM"}
  \item{delta}{voxel size in mm}
  \item{series}{series identifier}
  \item{image}{image number within series}
  \item{dim}{dimension of the data if available}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\references{http://medical.nema.org}
\note{ Since the DICOM standard is rather complicated, there may be
  cases where this function cannot read a DICOM file. Known issue: it
  cannot read header with implicit VR. Return value may change in future
  version!
}
\seealso{ \code{\link{read.AFNI}}, \code{\link{read.ANALYZE}} }
\examples{
\dontrun{dicom <- read.DICOM("dicomfile")}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
