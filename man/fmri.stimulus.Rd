\name{fmri.stimulus}
\alias{fmri.stimulus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Create the expected BOLD response for a given task indicator function.
}
\usage{
  fmri.stimulus(scans = 1, onsets = c(1), length = 1, rt = 3,
                mean = TRUE,
                a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, cc = 0.35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scans}{ number of scans }
  \item{onsets}{ vector of onset times (in scans) }
  \item{length}{ length of ON stimulus (in scans) }
  \item{rt}{ time between scans in seconds (TR) }
  \item{mean}{ logical. if TRUE the mean is substracted from the resulting vector }
  \item{a1}{ parameter of the hemodynamic response function (see details) }
  \item{a2}{ parameter of the hemodynamic response function (see details) }
  \item{b1}{ parameter of the hemodynamic response function (see details) } 
  \item{b2}{ parameter of the hemodynamic response function (see details) }
  \item{cc}{ parameter of the hemodynamic response function (see details) }
}
\details{
  The functions calculates the expected BOLD response for the task
  indicator function given by the argument as a convolution with the
  hemodynamic response function. The latter is modelled by the difference
  between two gamma functions as given in the reference (with the defaults
  for a1, a2, b1, b2, cc given therein):
  
  \deqn{\left(\frac{t}{d_1}\right)^{a_1} \exp \left(-\frac{t-d_1}{b_1}\right) 
- c \left(\frac{t}{d_2}\right)^{a_2} \exp
\left(-\frac{t-d_2}{b_2}\right) }{(x/d1)^a1 * exp(-(x - d1)/b1) - c * (x/d2)^a2 * exp(-(x - d2)/b2)}

  The parameters of this function can be changed through the arguments
  \code{a1}, \code{a2}, \code{b1}, \code{b2}, \code{cc}.
  
  The dimension  of the function value is set to \code{c(scans,1)}.

  If \code{mean} is TRUE (default) the resulting vector is
  corrected to have zero mean.
}
\value{
  Vector with dimension \code{c(scans,1)}.
}
\references{ Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1-15. }
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.lm}} }
\examples{
  # Example 1
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf,2)
  par(mfrow=c(2,2))
  for (i in 1:4) plot(z[,i],type="l")

}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line