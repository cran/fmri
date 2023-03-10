\name{slicetiming}
\alias{slicetiming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
slicetiming for fmridata-objects
}
\description{
Perform slicetiming for fMRI data, ideally before preprocessing (registration).
Recording times for slices are assumed to be equispaced between scans with
argument \code{sliceorder} providing the order of slice acquisitions.
Interpolation between slices is performed using a sinc filter.
}
\usage{
slicetiming(fmridataobj, sliceorder = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmridataobj}{
 object of class fmridata
}
  \item{sliceorder}{
 order of lice acquisitions
}
}

\value{
 an object of class fmridata
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

\seealso{
  \code{\link{fmri.stimulus}}, \code{\link{fmri.design}},\code{\link{fmri.lm}},~~~
}
\examples{
  \dontrun{
# Example 1
  data <- list(ttt=writeBin(rnorm(32*32*32*107), raw(), 4),
               mask=array(TRUE, c(32, 32, 32)), dim=c(32, 32, 32, 107))
  class(data) <- "fmridata"
  data <- slicetiming(data,sliceorder=1:32)
  ## provides data corrected for sequential slice acquisition in linear order
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip  }% use one of  RShowDoc("KEYWORDS")
\keyword{utilities }% __ONLY ONE__ keyword per line
