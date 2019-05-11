\name{sincfilter}
\alias{sincfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for sinc-interpolation
}
\description{
Performs sinc interpolation for a equidistant time series \code{x}
to times \code{t}.
}
\usage{
sincfilter(t, x, wr=8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
  vector of new time points
}
  \item{x}{
  observed time series at times \code{1:length(x)}.
}
 \item{wr}{
  determines truncation of series expansion
 }
}
\value{
a vector of interpolated values of the time series at time points given in
\code{t}.
}
\author{
         Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{slicetiming}}
}
\examples{
  x <- 1:107
  y <- rnorm(x)
  z <- sincfilter(seq(1,107,.01),y)
  plot(x, y, ylim=range(y,z))
  lines(seq(1,107,.01),z,col=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip  }% use one of  RShowDoc("KEYWORDS")
\keyword{utilities }% __ONLY ONE__ keyword per line
