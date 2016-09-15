\name{fmri.stimulus}
\alias{fmri.stimulus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Model for FMRI Data }
\description{
  Create the expected BOLD response for a given task indicator function.
}
\usage{
  fmri.stimulus(scans = 1, onsets = c(1), durations = c(1), TR = 2,
                times = FALSE, type = c("canonical", "gamma", "boxcar", "user"),
                par = NULL, scale = 10, hrf = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scans}{ number of scans }
  \item{onsets}{ vector of onset times (in scans) }
  \item{durations}{ vector of duration of ON stimulus in scans or
    seconds (if \code{!is.null(times)})  }
  \item{TR}{ time between scans in seconds (TR) }
  \item{times}{ onset times in seconds. If present \code{onsets}
    arguments is ignored. }
  \item{type}{One of \code{"canonical"}, \code{"gamma"}, \code{"boxcar"}, 
     \code{"user"}}
  \item{par}{ Possible parameters to the HRF. }
  \item{scale}{ Temporal undersampling factor }
  \item{hrf}{ If \code{type} is \code{"user"} this should be a function evaluating the hemodynamic response function }
  \item{verbose}{ Report more if \code{TRUE} }
}
\details{
  The functions calculates the expected BOLD response for the task
  indicator function given by the argument as a convolution with the
  hemodynamic response function. 
  
  For \code{type} is \code{"canonical"} the latter is modelled by the difference
  between two gamma functions as given in the reference (with the defaults
  for a1, a2, b1, b2, cc given therein):
  
  \deqn{\left(\frac{t}{d_1}\right)^{a_1} \exp \left(-\frac{t-d_1}{b_1}\right) 
- c \left(\frac{t}{d_2}\right)^{a_2} \exp
\left(-\frac{t-d_2}{b_2}\right) }{(x/d1)^a1 * exp(-(x - d1)/b1) - c * (x/d2)^a2 * exp(-(x - d2)/b2)}

  The parameters \code{a1}, \code{a2}, \code{b1}, \code{b2}, \code{cc} of this function 
  can be changed through the argument \code{par} in this order.

  Other choices are a simple gamma function 
  
  \deqn{\frac{1}{k\tau_h (k-1)!} \left( \frac{t}{\tau_h} \right)^k 
        \exp \left( - \frac{t}{\tau_h} \right)}
        
  or the \code{"boxcar"} stimulus, or a user defined function \code{hrf}.

  The dimension  of the function value is set to \code{c(scans, 1)}.

  If \code{!is.null(times)} durations are specified in seconds. 

}
\value{
  Vector with dimension \code{c(scans, 1)}.
}
\references{ Worsley, K.J., Liao, C., Aston, J., Petre, V., Duncan,
  G.H., Morales, F., Evans, A.C. (2002). A general statistical analysis
  for fMRI data. NeuroImage, 15:1-15. 

  Polzehl, J. and Tabelow, K. (2007)
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\seealso{ \code{\link{fmri.design}}, \code{\link{fmri.lm}} }
\examples{
  # Example 1
  hrf <- fmri.stimulus(107, c(18, 48, 78), 15, 2)
  z <- fmri.design(hrf, 2)
  par(mfrow=c(2, 2))
  for (i in 1:4) plot(z[, i], type="l")

}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ design }% __ONLY ONE__ keyword per line
