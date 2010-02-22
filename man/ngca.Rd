\name{ngca}
\alias{ngca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Non-Gaussian Component Analysis }
\description{
  The function performs Non-Gaussian Component Analysis as described in Blanchard et.al. (2005).
}
\usage{
ngca(data, L=c(1000,1000,1000), T=10, m=3, eps=1.5, npca=min(dim(x)[2],
     dim(x)[1]), filter.time="None", filter.space=FALSE, method="temporal",
     dg.trend = 2, h.space=3, h.time=3, keepv=TRUE, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Observation matrix (dimension Nxd) }
  \item{L}{ Number basis functions in each of four classes. }
  \item{T}{ Number of Fast ICA iterations  }
  \item{m}{ Number of non-Gaussian components. }
  \item{eps}{ Threshold (defaults to 1.5)}
  \item{npca}{ Reduce space to \code{npca} principal components. This can be 
               used to avoid standardizing by numerically singular covariance matrices.
               In fMRI this allows to reduce the dimensionality assuming that
               the interesting non-Gaussian directions are also characterized by
               larger variances.}
  \item{filter.time}{ Choice of temporal filtering before analysis: \code{"None"}, \code{"Low"}, \code{"Both"}, \code{"High"} (default \code{"None"})}
  \item{filter.space}{Choice of spatial filtering before analysis: logical, default \code{FALSE}}
  \item{method}{ Either \code{"spatial"} or \code{"temporal"}. Specifies the type of NGCA to perform. }
  \item{dg.trend}{ not yet documented }
  \item{h.space}{ bandwidth for spatial filtering. default 3}
  \item{h.time}{ bandwidth for temporal filtering. default 3}
  \item{keepv}{if \code{TRUE} intermediate results from fast ICA step are kept.}
  \item{delta}{ not yet documented }
}
\details{
  The function performs Non-Gaussian Component Analysis as described in Blanchard et.al. (2006).
  The procedure uses four classes of basis functions, i.e. Gauss-Power3, Hyperbolic
  Tangent and the real and complex part of the Fourier class. See 
  Blanchard et.al. (2005) for details.
}
\value{ The function returns a list with components
    \item{ihat }{ Matrix containing the first m NGCA directions as columns.}
    \item{sdev }{ Standard deviations of the principal components of 
                 the thresholded ICA directions}
    \item{xhat }{first m components of the rotated data}
    \item{v }{If \code{keepv==TRUE} the set of directions \code{v^{(k)}}}
    \item{normv }{If \code{keepv==TRUE} the norm of each \code{v^{(k)}}.}
  ...
}
\references{ Blanchard, G., Kawanabe, M., Sugiyama, M., Spokoiny, V. and M\"uller 
K.-R. (2005). In Search of Non-Gaussian Components of a High-Dimensional Distribution.
Journal of Machine Learning Research. pp. 1-48. }
\author{ J\"org Polzehl polzehl@wias-berlin.de }
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
