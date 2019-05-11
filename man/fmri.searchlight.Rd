\name{fmri.searchlight}
\alias{fmri.searchlight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Searchlight signal detection
}
\description{
  Detection of activated regions using searchlights.
}
\usage{
fmri.searchlight(spm, alpha = 0.05, radius, minimum.signal = 0,
   kind = c("abs", "squared"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spm}{
  \code{fmrispm} object
}
  \item{alpha}{
  multiple test (over volume) adjusted significance level.
}
  \item{radius}{
  radius of searchlight. Value needs to be larger or eaqual than 1.
}
  \item{minimum.signal}{
  allows to specify a (positive) minimum value for detected signals. If minimum.signal >0
   the thresholds are to conservative, this case needs further improvements.
}
  \item{kind}{
  Kind of statistics used for aggregation over search light region.
  \code{"abs"} specifies averaging of absolute voxelwise t-statistics
  while \code{"squared"} corresponds to averaging of squares of these statistics.
}
}
\details{
   The function computes mean statistics (depending on \code{kind}) over a
   searchlight region of radius \code{radius}.
   Approximate voxelwise p-values are determined with respect an empirical
   (simulated) distribution of the searchlight statistics under the null hypothesis
   a central t-distributed spm. Thresholding used \code{FDR} with rate \code{alpha}.
}
\value{
Object with class attributes "fmripvalue" and "fmridata"
\item{pvalue}{ voxelwise p-value if exceeding FDR-critical value, 1 otherwise.}
\item{weights}{ voxelsize ratio }
\item{dim}{ data dimension }
\item{hrf}{ expected BOLD response for contrast (single stimulus only)}
}
\references{
Kriegeskorte, N.; Goebel, R. & Bandettini, P. (2006) \emph{Information-based functional brain mapping}, PNAS 103:3863-3868.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{
\code{\link{fmri.lm}}, \code{\link{fmri.pvalue}}, \code{\link{fmri.cluster}}
}
\examples{
  \dontrun{fmri.searchlight(fmrispmobj)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
