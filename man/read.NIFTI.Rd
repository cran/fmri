\name{read.NIFTI}
\alias{read.NIFTI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O Functions }
\description{
  Read fMRI data from NIFTI file(s).
}
\usage{
read.NIFTI(filename, level = 0.75, mask=NULL, setmask=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the NIfTI file }
  \item{level}{Quantile level defining the mask}
  \item{mask}{array or nifti-object containing the mask. If set this replaces
  the mask defined by argument level.}
  \item{setmask}{Logical (default \code{TRUE}), whether to define a suitable mask based on \code{level}}
}
\details{
  This function reads fMRI data files in NIfTI format.

  The filename can be given with or without extension. If extension is
  not included, the function searches for the ".nii" file and then for
  the "hdr/img" pair.
}
\value{
  Object of class "fmridata" with the following list entries:
  \item{ttt}{raw vector (numeric size 4) containing the four dimensional data cube (the first three dimensions are voxel
    dimensions, the fourth dimension denotes the time).}
  \item{header}{header information of the data}
  \item{format}{data source. string "NIFTI"}
  \item{delta}{voxel size in mm}
  \item{origin}{position of the datacube origin}
  \item{orient}{data orientation code}
  \item{dim}{dimension of the datacube}
  \item{weights}{weights vector coding the relative voxel sizes in x,
    y, z-direction}
  \item{mask}{head mask}
}
\references{
  Polzehl, J. and Tabelow, K. (2007)
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{read.ANALYZE}}, \code{\link{read.AFNI}} }
\examples{ \dontrun{analyze <- read.NIFTI("niftifile.nii")}}
\keyword{ IO }% at least one, from doc/KEYWORDS
\keyword{ utilities }% __ONLY ONE__ keyword per line
