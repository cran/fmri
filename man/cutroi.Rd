\name{cutroi}
\alias{cutroi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O function }
\description{
  This functions cuts a region-of-interest (ROI) from input data.
}
\usage{
cutroi(data, xind = 1:data$dim[1], yind = 1:data$dim[2], 
             zind = 1:data$dim[3], tind = 1:data$dim[4])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Object of class fmridata. }
  \item{xind}{ vector of roi-indices for first data index }
  \item{yind}{ vector of roi-indices for second data index }
  \item{zind}{ vector of roi-indices for third data index }
  \item{tind}{ vector of roi-indices for 4th data index }
}
\details{
  Cut a region of interest from frmidata.
}
\value{
  Corresponding cutted fmridata object.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.AFNI}}, \code{\link{read.ANALYZE}}, \code{\link{read.NIFTI}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ utilities }% __ONLY ONE__ keyword per line
