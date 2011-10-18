\name{summary.fmridata}
\alias{summary.fmridata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  'summary' method for class '"fmridata"'.
}
\usage{
\method{summary}{fmridata}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{fmridata}, usually, a result of a call to
          \code{fmri.lm}, \code{fmri.smooth}, \code{fmri.pvalue},
	  \code{read.AFNI}, or \code{read.ANALYZE}. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  The method tries to print information on data, like data dimension,
  voxel size, value range.
}
\value{
  A list with the following elements:
  \item{dim }{data dimension}
  \item{delta }{voxel dimension, if available}
  \item{values }{value range}
  \item{z}{ design matrix }
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{print.fmridata}} }
\examples{
\dontrun{summary(data)}
}
\keyword{ print }% at least one, from doc/KEYWORDS
