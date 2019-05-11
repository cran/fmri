\name{plot.fmripvalue}
\alias{plot.fmripvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Visualize fMRI p-value maps
}
\description{ Visualize objects created by function \code{fmri.pvalue}
}
\usage{
\method{plot}{fmripvalue}(x, template = NULL, mask = NULL,
      view = c("axial", "coronal", "sagittal", "orthographic"),
      slices = NULL, ncol = 1, nrow = 1, center = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   object of class 'fmripvalue'
}
  \item{template}{
   Anatomical image of same origin and direction
   as pvalue map in x$pvalue.
}
  \item{mask}{
   optional brain mask
}
  \item{view}{
   Either 'orthographic' or one of 'axial', 'coronal' or 'sagittal'
}
  \item{slices}{
   If \code{view != "orthographic"} vector of slice numbers to use.
   If not provided the \code{ncol*nrow} slices with strongest signals are selected
}
  \item{ncol}{
  If \code{view != "orthographic"} number of slices per row
}
  \item{nrow}{
  If \code{view != "orthographic"} number of rows in display.
}
  \item{center}{
  If \code{view == "orthographic"} center of orthographic view.
  If not provided the center is chosen to provide maximal information.
}
  \item{\dots}{
  additional parameters (not evaluated)
}
}
\value{
list with components
\item{comp1 }{slices, numbers refer to spm}
\item{comp2 }{center, numbers refer to spm}
}
\author{
Joerg Polzehl \email{polzehl@wias-berlin.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmri.pvalue}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
