\name{read.AFNI}
\alias{read.AFNI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O function }
\description{
  Read HEAD/BRIK file.
}
\usage{
read.AFNI(filename,vol=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the file (without extension) }
  \item{vol}{vector of volumes of the dataset to be read}
}
\details{
  The function reads a HEAD/BRIK file. If \code{vol} is given (defaults to
  \code{NULL}), only volumes in this vector are read, in order to save 
  memory.
}
\value{
  Object of class "fmridata" with the following list entries:
  \item{ttt}{raw vector (numeric size 4) containing the four dimensional data cube (the first three dimensions are voxel
    dimensions, the fourth dimension denotes the time).}
  \item{header}{header information list}
  \item{format}{data source. string "HEAD/BRIK"}
  \item{delta}{voxel size in mm}
  \item{origin}{position of the datacube origin}
  \item{orient}{data orientation code. see AFNI documentation}
  \item{dim}{dimension of the datacube}
  \item{weights}{weights vector coding the relative voxel sizes in x,
    y, z-direction.}
  \item{mask}{head mask}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\references{R. W. Cox (1996). AFNI: Software for analysis and visualization of
  functional magnetic resonance neuroimages. Computers and Biomed. Res. 
  29:162-173.
}
\seealso{ \code{\link{write.AFNI}}, \code{\link{read.ANALYZE}} }
\examples{
\dontrun{afni <- read.AFNI("afnifile")}
}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
