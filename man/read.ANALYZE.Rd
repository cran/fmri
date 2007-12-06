\name{read.ANALYZE}
\alias{read.ANALYZE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Read fMRI data from ANALYZE file(s).
}
\usage{
read.ANALYZE(prefix = "", numbered = FALSE, postfix = "",
             picstart = 1, numbpic = 1, level = 0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prefix}{ string(s). part of the file name before the
    number or vector of strings for filename (if \code{numbered} is \code{FALSE}) }
  \item{numbered}{ logical. if \code{FALSE} only \code{prefix} is taken
    as file name (default).}
  \item{postfix}{ string. part of the file name after the number}
  \item{picstart}{ number of the first image to be read. }
  \item{numbpic}{ number of images to be read }
  \item{level}{Quantile level defining the mask}
}
\details{
  This function reads fMRI data files in ANALYZE format.
  If \code{numbered} is \code{FALSE}, only the vector of strings in \code{prefix}
  is used for file name (default).

  If \code{numbered} is \code{TRUE}, it takes the first string in \code{prefix} and \code{postfix} and
  a number of the form "007" in between to create the file name.

  The number is assumed to be 3 digits (including leading zeros). First
  number is given in \code{picstart}, while \code{numbpic} defines the
  total number of images to be read. Data in multiple files
  will be combined into a four dimensional datacube.
}
\value{
  Object of class "fmridata" with the following list entries:
  \item{ttt}{raw vector (numeric size 4) containing the four dimensional data cube (the first three dimensions are voxel
    dimensions, the fourth dimension denotes the time).}
  \item{header}{header information of the data}
  \item{format}{data source. string "ANALYZE"}
  \item{delta}{voxel size in mm}
  \item{origin}{position of the datacube origin}
  \item{orient}{data orientation code}
  \item{dim}{dimension of the datacube}
  \item{weights}{weights vector coding the relative voxel sizes in x,
    y, z-direction}
  \item{mask}{head mask}
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\references{Biomedical Imaging Resource (2001). Analyze Program. Mayo Foundation.
}
\note{   Since numbering and naming of ANALYZE files widely vary, this
  function may not meet your personal needs. See Details section above for 
  a description.
}
\seealso{ \code{\link{write.ANALYZE}}, \code{\link{read.AFNI}} }
\examples{ \dontrun{analyze <- read.ANALYZE("analyze",TRUE,"file",31,107)}}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
