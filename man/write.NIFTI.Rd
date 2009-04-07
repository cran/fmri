\name{write.NIFTI}
\alias{write.NIFTI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write a 4 dimensional datacube in NIFTI file format.
}
\usage{
write.NIFTI(ttt, header=NULL, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ 4 dimensional datacube }
  \item{header}{ header information }
  \item{filename}{ file name }
}
\details{
  Writes the datacube \code{ttt} to a file named \code{file} in NIFTI file
  format. \code{header} is a list that contains the header
  information.
  
  See NIFTI documentation for details.
}
\value{
  Nothing is returned.
}
\references{
  Polzehl, J. and Tabelow, K. (2007)
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.ANALYZE}}, \code{\link{write.AFNI}} }
\examples{
## Example 1
write.NIFTI(array(as.integer(65526*runif(10*10*10*20)),c(10,10,10,20)),
              file="niftifile")
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
