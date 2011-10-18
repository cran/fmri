\name{print.fmridata}
\alias{print.fmridata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  'print' method for class '"fmridata"'.
}
\usage{
\method{print}{fmridata}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{fmridata}, usually, a result of a call to
          \code{fmri.lm}, \code{fmri.smooth}, \code{fmri.pvalue},
	  \code{read.AFNI}, or \code{read.ANALYZE}. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  The method tries to print information on data, like data dimension,
  voxel size, value range.
}
\value{
  none
}
\references{
  Polzehl, J. and Tabelow, K. (2007)
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{summary.fmridata}} }
\examples{
\dontrun{print(data)}
}
\keyword{ print }% at least one, from doc/KEYWORDS
