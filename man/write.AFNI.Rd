\name{write.AFNI}
\alias{write.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Write BRIK/HEAD files.
}
\usage{
write.AFNI(filename, ttt, label = NULL, note = NULL, origin = NULL, 
               delta = NULL, idcode = NULL, header = NULL, taxis = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the file }
  \item{ttt}{ datacube }
  \item{label}{ labels (BRICK\_LABS), depreciated - see header}
  \item{note}{ notes on data (HISTORY\_NOTE), depreciated - see header }
  \item{origin}{ origin of datacube (ORIGIN), depreciated - see header }
  \item{delta}{ voxel dimensions (DELTA), depreciated - see header }
  \item{idcode}{ idcode of data (IDCODE\_STRING), depreciated - see header }
  \item{header}{ This is a list of header information such as 
                 DATASET\_RANK to be written to the .HEAD file.
                 Arguments \code{label}, ... are depreciated and to be substituted
                 by a corresponding list entry. For backward compatibility
                 the use of the old arguments is still supported and should 
                 give the same results. This will be removed in some future 
                 release! Since AFNI does not read any dataset with a header choose 
                 carefully what is written. There are some basic tests in this function,
                 but this may not be sufficient.}
  \item{taxis}{ logical (defaults to \code{FALSE}. Are the sub-bricks time series? 
                This results in writing TAXIS attributes to the header file.}
}
\details{
  Write out BRIK/HEAD files as required by AFNI. Most arguments
  correspond to entries in the HEAD file, but use is depreciated. Use \code{header} and \code{taxis} instead!
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
\seealso{ \code{\link{read.AFNI}},\code{\link{write.ANALYZE}} }
\examples{
 \dontrun{write.AFNI("afnifile", array(as.integer(65526*runif(10*10*10*20)),
     c(10,10,10,20)), c("signal"), note="random data",
     origin=c(0,0,0), delta=c(4,4,5), idcode="unique ID")}
 write.AFNI("afnifile", array(as.integer(65526*runif(10*10*10*20)),
     c(10,10,10,20)), header=list(HISTORY_NOTE="random data",
     ORIGIN=c(0,0,0), DELTA=c(4,4,5), IDCODE_STRING="unique ID"),taxis=FALSE)
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
