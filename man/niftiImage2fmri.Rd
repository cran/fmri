\name{niftiImage2fmri}
\alias{niftiImage2fmri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create fmridata object from niftiImage
}
\description{
Transforms a niftiImage (created by readNifti from package RNiftyReg) into
an object with class fmridata
}
\usage{
niftiImage2fmri(niftiobj, level = 0.75, mask=NULL, setmask = TRUE, indx = NULL,
   indy = NULL, indz = NULL, avoidnegs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niftiobj}{
     an object of class niftiImage
     }
  \item{level}{
    quantile used in mask definition
}
\item{mask}{array or nifti-object containing the mask. If set this replaces
the mask defined by argument level.}
  \item{setmask}{
    if \code{TRUE} create a brain mask
}
  \item{indx}{
    index vector for subcube definition
}
  \item{indy}{
    index vector for subcube definition
}
  \item{indz}{
    index vector for subcube definition
}
  \item{avoidnegs}{
    if \code{TRUE} change the mean to avoid negative image intensities
}
}
\details{
   This function can be used in connection with readNifti from package RNiftyReg
   to read large fMRI series from nifti files. The resulting fmridata-object
   stores the image data as 2 byte integer in raw format, in contrast for the
   4 byte real used with other functions.
}
\value{
   an object of class fmridata
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }

\seealso{
\code{\link{read.AFNI}}, \code{\link{read.DICOM}},
\code{\link{read.ANALYZE}}, \code{\link{read.NIFTI}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }% use one of  RShowDoc("KEYWORDS")
\keyword{ utilities }% __ONLY ONE__ keyword per line
