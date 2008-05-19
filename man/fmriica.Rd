\name{fmriica}
\alias{fmriica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Independent Component Analysis }
\description{
  The function performs Independent Component Analysis.
}
\usage{
fmriica(data, m = 3, method = "spatial", xind = NULL, yind = NULL, zind = NULL, tind = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{  Observation matrix (dimension Nxd) }
  \item{m}{ Number of independent components. }
  \item{method}{ Either \code{"spatial"} or \code{"temporal"}. Specifies the type of ICA to perform.  }
  \item{xind}{ index of x-coordinates to use }
  \item{yind}{ index of y-coordinates to use }
  \item{zind}{ index of z-coordinates to use }
  \item{tind}{ index of time points to use }
  \item{\dots}{ further arguments to fastICA }
}
\details{
  This is still experimental code based on the package fastICA. The package fastICA seems limited in
the data-size it can handle. \code{xind}, \code{yind}, \code{zind} and \code{tind} may be use to restrict
the analysis to a cube in space and certain time points.
}
\value{The function returns a list with components
    \item{ihat }{ Matrix containing the first m ICA directions as columns.}
    \item{sdev }{ Standard deviations of the principal components of 
                 the thresholded ICA directions}
    \item{xhat }{first m components of the rotated data}
    \item{v }{If \code{keepv==TRUE} the set of directions \code{v^{(k)}}}
    \item{normv }{If \code{keepv==TRUE} the norm of each \code{v^{(k)}}.}
}
\author{ J\"org Polzehl polzehl@wias-berlin.de  }
\seealso{ \code{\link{ngca}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }
