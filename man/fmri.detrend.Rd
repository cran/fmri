\name{fmri.detrend}
\alias{fmri.detrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Detrend fMRI time series }
\description{
  Detrend fMRI dataset with a polynomial of given degree
}
\usage{
fmri.detrend(data, degree = 1, nuisance=NULL, accoef = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ fMRI dataset of class ''\code{fmridata}'' }
  \item{degree}{ Degree of the polynomial used to detrend the data. defaults to 1 (linear trends). }
  \item{nuisance}{ Matrix of additional nuisance parameters to regress against. }
  \item{accoef}{ Coefficient of AR(1) model used for prewhitening. default 0. }
}
\details{
  The function can be used to detrend the time series of an fMRI dataset \code{data} (of class ''\code{fmridata}'' using polynomials. If the argument \code{degree} is larger than 0 (default: 1) the polynomial trends up to the given degree are removed from the data. If the argument \code{accoef} is larger than 0 (default: 0) prewhitening using an AR(1) model is performed.
}
\value{
  Detrended data object of class ''\code{fmridata}''.
}
\references{
  Polzehl, J. and Tabelow, K. (2007).
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17.
}
\author{ Joerg Polzehl \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{fmri.lm}} }
\examples{
  # Example 1
  data <- list(ttt=writeBin(rnorm(32*32*32*107),raw(),4),
               mask=array(1,c(32,32,32)),dim=c(32,32,32,107))
  class(data) <- "fmridata"
  data <- fmri.detrend(data,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
