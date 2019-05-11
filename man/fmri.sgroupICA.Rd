\name{fmri.sgroupICA}
\alias{fmri.sgroupICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial group ICA for fmri
}
\description{
   Combine ICA results from multiple runs or multiple subjects
   in group fMRI studies
}
\usage{
fmri.sgroupICA(icaobjlist, thresh = 0.75, minsize=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{icaobjlist}{
    List of results obtained by function \code{\link{fmri.sICA}}
    for a series of fmri data sets (multiple runs or multiple subjects).
}
  \item{thresh}{
    threshold for cluster aggregation. Needs to be in (0,1).
}
\item{minsize}{
  Minimal size of cluster to consider in IC aggregation. Needs to be
  larger than 1.
}
}
\details{
The fMRI time series need to be preprocessed and registered
before thr ICA decomposition is performed.

The function employes a hierarchical clustering algorithm (complete
linkage) on the combined set of spatial independent components
obtained from the individual time series. A distance matrix is obtained from
correlations of the independent component images. Aggregation of two components from
the same fmri series is prevented in the algorithm.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
An object of class ''\code{fmrigroupICA}'' with components
  \item{icacomp}{Mean IC's over cluster members for cluster of size larger
  or equal \code{minsize}}
  \item{size}{Size of selected clusters}
  \item{cl}{Number of selected clusters}
  \item{cluster}{Cluster membership corresponding to \code{thresh}.}
  \item{height}{Distance value at which the cluster was createt.
  Elements correspond to elements of cluster.}
  \item{hdm}{Object returned by function \code{\link[stats]{hclust}}.}
}
\references{
F. Esposito et al (2005) Independent component analysis of fMRI group studies by
self-organizing clustering, Neuroimage, pp. 193-205.
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmri.sICA}}, \code{\link{plot.fmrigroupICA}}, \code{\link[stats]{hclust}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
