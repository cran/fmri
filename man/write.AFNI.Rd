\name{write.AFNI}
\alias{write.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Write BRIK/HEAD files.
}
\usage{
write.AFNI(filename, ttt, label, note = "", origin = c(0, 0, 0), 
               delta = c(4, 4, 4), idcode = "WIAS_noid")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the file }
  \item{ttt}{ datacube }
  \item{label}{ labels (BRICK\_LABS) }
  \item{note}{ notes on data (HISTORY\_NOTE) }
  \item{origin}{ origin of datacube (ORIGIN) }
  \item{delta}{ voxel dimensions (DELTA) }
  \item{idcode}{ idcode of data (IDCODE\_STRING) }
}
\details{
  Write out BRIK/HEAD files as required by AFNI. Most arguments
  correspond to entries in the HEAD file. 
}
\value{
  Nothing is returned.
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.AFNI}},\code{\link{write.ANALYZE}} }
\examples{
 write.AFNI("afnifile", array(as.integer(65526*runif(10*10*10*20)),
     c(10,10,10,20)), c("signal"), note="random data",
     origin=c(0,0,0), delta=c(4,4,5), idcode="unique ID")
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
