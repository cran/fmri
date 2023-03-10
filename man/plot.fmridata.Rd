\name{plot.fmridata}
\alias{plot.fmridata}
\alias{plot.fmrisegment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ I/O functions }
\description{
  Visualize fMRI data and (intermediate) results.
}
\usage{
\method{plot}{fmridata}(x, anatomic = NULL, maxpvalue = 0.05,
              spm = TRUE, pos = c(-1, -1, -1), type = "slice",
              slice =  1, view = "axial" ,zlim.u =
              NULL, zlim.o = NULL,col.o = heat.colors(256), col.u =
              grey(0:255/255), cutOff = c(0, 1), ...)
\method{plot}{fmrisegment}(x, anatomic = NULL,
              slice =  1, view = c( "axial", "coronal", "sagittal") ,zlim.u =
              NULL, zlim.o = NULL,col.o = c( rainbow( 64, start = 2/6, end = 4/6),
              rainbow( 64, start = 0, end = 1/6)),
              col.u = grey(0:127/127), verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class "fmrisegment", "fmrispm" or "fmridata" }
  \item{anatomic}{ overlay of same dimension as the functional data, or fmridata object
                   (if of x is fmripvalue object) }
  \item{maxpvalue}{ maximum p-value for thresholding }
  \item{spm}{ logical. if class is "fmrispm" decide whether to plot
    the t-statistics for the estimated effect (\code{spm=TRUE}) or the estimated effect
    itself (\code{spm=FALSE}).}
  \item{pos}{ voxel to be marked on output }
  \item{type}{ string. "slice" for slicewise view and "3d" for
    3d view.}
  \item{slice}{ number of slice in x, if anatomic is of "fmridata" class}
  \item{view}{ "axial", "coronal", or "sagittal", if anatomic is of "fmridata" class }
  \item{zlim.u}{ full range for anatomical underlay used for color scale, if anatomic is of "fmridata" class }
  \item{zlim.o}{ full range for functional overlay used for color scale, if anatomic is of "fmridata" class }
  \item{col.u}{ color scale for anatomical underlay, if anatomic is of "fmridata" class, default grey(0:255/255) }
  \item{col.o}{ color scale for functional overlay, if anatomic is of "fmridata" class, default heat.colors(256) }
  \item{cutOff}{ not yet documented }
  \item{verbose}{ tell something on the progress?}
  \item{...}{ additional arguments for plot}
}
\details{
  Provides a slicewise view of "fmridata" objects with anatomic overlay (if
  appropriate, that is for class "fmripvalue"). For objects of class
  "fmrispm" it plots the t-statistics for the estimated effects if \code{spm} is
  \code{TRUE}, or the estimated effect otherwise. For objects of class
  "fmridata" only a plot of the data slices itself is produced. If \code{device} is
  specified as  "png", "jpeg", "ppm" output is done to a file. A grey/color scale
  is provided in the remaining space.

  For objects of class "fmrisegment" the smoothed signal size is shown in the
  activation segments (two-sided test!).

  If \code{type} is "3d" a 3 dimensional interactive view opens. Sliders
  to move in the data cube are given ("x", "y", "z", and "t" if class is
  "fmridata" only). Time series are shown if available. For objects
  of class "fmrispm" a slider is created to remove  information for voxels with
  smaller signals than a cut-off value from the plot.
  Use pvalues for statistical evaluation. If \code{spm} is
  \code{FALSE} the estimated BOLD response together with a confidence
  interval corresponding to \code{maxpvalue} is drawn. For objects of class
  "fmripvalue" the pvalues with overlay are shown.
}
\value{
  If 'type' is "3d" the Tk-object is returned. (Remove the display with \code{tkdestroy(object)})
}
\references{
  Polzehl, J. and Tabelow, K. (2007)
  \emph{fmri: A Package for Analyzing fmri Data},
  R News, 7:13-17 .
}
\author{ Karsten Tabelow \email{tabelow@wias-berlin.de} }
\note{ 3 dimensional plotting requires the \code{tkrplot}
  package.
}
\seealso{ \code{\link{fmri.pvalue}} }
\examples{
\dontrun{plot(pvalue)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ iplot }% at least one, from doc/KEYWORDS
