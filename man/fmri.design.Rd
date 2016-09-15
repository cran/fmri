\name{fmri.design}
\alias{fmri.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Return a design matrix for a linear model with given stimuli and
  possible polynomial drift terms.
}
\usage{
  fmri.design(stimulus, order = 2, cef = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stimulus}{ matrix containing expected BOLD repsonse(s) for the linear model as columns}
  \item{order}{ order of the polynomial drift terms }
  \item{cef}{ confounding effects }
  \item{verbose}{ Report more if \code{TRUE} }
}
\details{
  The stimuli given in \code{stimulus} are used as first columns in
  the design matrix.

  The order of the polynomial drift terms is given
  by \code{order}, which defaults to 2.

  Confounding effects can be included in a matrix \code{cef}.

  The polynomials are defined orthogonal to the stimuli given in
  \code{stimulus}. 
}
\value{
  design matrix of the linear model
}
\references{
  Polzehl, J. and Tabelow, K.(2007).
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmri.stimulus}}, \code{\link{fmri.lm}} }
\examples{
  # Example 1
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf, 2)
  par(mfrow=c(2, 2))
  for (i in 1:4) plot(z[, i], type="l")
}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line
