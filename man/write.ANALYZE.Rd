\name{write.ANALYZE}
\alias{write.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Write a 4 dimensional datacube in ANALYZE file format.
}
\usage{
write.ANALYZE(ttt, header=NULL, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttt}{ 4 dimensional datacube }
  \item{header}{ header information }
  \item{filename}{ file name }
}
\details{
  Writes the datacube \code{ttt} to a file named \code{file} in ANALYZE file
  format. \code{header} is a list that contains the header information
  as documented by the Mayo Foundation. We give here a short summary. If
  a value is not provided, it will be tried to fill it with reasonable
  defaults, but do not expect fine results, if the entry has a special
  important meaning (h.i. pixdim).
  \tabular{rlrl}{
 [1] \tab datatype1 -- 10 byte character \tab
 [2] \tab dbname -- 18 byte character\cr
 [3] \tab extents -- integer \tab
 [4] \tab sessionerror -- integer\cr
 [5] \tab regular -- character \tab
 [6] \tab hkey -- character\cr
 [7] \tab dimension -- 8 integers, dimensions ... \tab
 [8] \tab unused -- 7 integers\cr
 [9] \tab datatype -- integer, datatype usually "4" \tab
 [10] \tab bitpix -- integer\cr
 [11] \tab dimun0 -- integer \tab
 [12] \tab pixdim -- 8 floats, voxel dimensions ...\cr
 [13] \tab voxoffset -- float \tab
 [14] \tab funused -- 3 floats\cr
 [15] \tab calmax -- float \tab
 [16] \tab calmin -- float\cr
 [17] \tab compressed -- float \tab
 [18] \tab verified -- float\cr
 [19] \tab glmax -- integer \tab
 [20] \tab glmin -- integer\cr
 [21] \tab describ -- 80 byte character \tab
 [22] \tab auxfile -- 24 byte character\cr
 [23] \tab orient -- character \tab
 [24] \tab originator -- 5 integers \cr
 [25] \tab generated -- 10 byte character \tab
 [26] \tab scannum -- 10 byte character\cr
 [27] \tab patientid -- 10 byte character \tab
 [28] \tab expdate -- 10 byte character\cr
 [29] \tab exptime -- 10 byte character \tab
 [30] \tab histun0 -- 3 byte character\cr
 [31] \tab views -- integer \tab
 [32] \tab voladded -- integer\cr
 [33] \tab startfield -- integer \tab
 [34] \tab fieldskip -- integer\cr
 [35] \tab omax -- integer \tab
 [36] \tab omin -- integer\cr
 [37] \tab smax -- integer \tab
 [38] \tab smin -- integer       }
 See ANALYZE documentation for details.
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
write.ANALYZE(array(as.integer(65526*runif(10*10*10*20)),c(10,10,10,20)),
              file=file.path(tempdir(),"analyzefile"))
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
