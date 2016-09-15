\name{fmri.designG}
\alias{fmri.designG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Design matrix for fMRI group analysis }
\description{
This function returns a design matrix for multi-subject fMRI data to fit a Linear Mixed-effects Model (one-stage procedure) with given stimuli, polynomial drift terms and a set of known population parameters.
}
\usage{
fmri.designG(hrf, subj = 1, runs = 1, group = NULL, XG = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrf}{ vector or matrix containing expected BOLD response(s) 
  for one session, typically a \code{\link{fmri.stimulus}} object. }
  \item{subj}{ number of subjects in the study. }
  \item{runs}{ number of repeated measures within subjects. }
  \item{group}{ optional vector to define groups. 
  It is expected one value per subject. A grouping factor can also be part of \code{XG}. }
  \item{XG}{ optionally, a group-level design matrix of class \code{"data.frame"}, which contains population parameters such as ages or gender corresponding to the subjects. It is expected one value per subject. }
}
\details{
Based on the dimensionality of the \code{hrf} object, which provides the total number of scans (time-points) within each session, the entered number of subjects and repeated measures the auxiliary variables: "subj", "run", "scan" and "session" are generated as first part of the returned design matrix.
  
If no \code{group} argument is specified, only one population will be assumed; otherwise group labels are replicated within sessions of the same subject. 
  
First a design matrix for a single run is created by calling: \code{x <- fmri.design(hrf, order = 2)}. Hence the polynomial drift terms are defined orthogonal to the stimuli (see \code{\link{fmri.design}}). This matrix is replicated blockwise to all sessions assuming the same experimental design for all runs. The first drift term, a column of ones, is called "drift0" and models an intercept.
  
If given, further subject characteristics are filled in the design matrix.
}
\value{
A design matrix as a data frame, which contains the following variables:
  \item{subj}{ consecutive subject number: 1 to \code{subj} specified as factor }
  \item{run}{ consecutive run number within the subjects: 1 to \code{runs} specified as factor }
  \item{scan}{ consecutive scan number: 1 to T within each session }
  \item{session}{ consecutive experiment number: 1 to \code{(subj*runs)} specified as factor }
  \item{group}{ grouping variable specified as factor, one group by default }
  \item{hrf, hrf2, ...}{ replicated expected BOLD-response(s) }
  \item{drift0, drift1, drift2}{ replicated polynomial drift terms 
    created with \code{fmri.design(hrf, order = 2)} orthogonal to the stimuli
    given in \code{hrf} }
  \item{...}{ further expanded between-subject factors and covariates }
}
\references{
Polzehl, J. and Tabelow, K.(2007). \emph{fmri: A Package for Analyzing fmri Data}, R News, 7:13-17.  
}
\author{ Sibylle Dames }
\seealso{
  \code{\link{fmri.stimulus}}, \code{\link{fmri.design}}, \code{\link{fmri.lmePar}}
}
\examples{
subj <- 6
runs <- 1
scans <- 121
times <- c(12, 48, 84, 120, 156, 192, 228, 264) 
duration <- 24          
tr <- 2.5

hrf <- fmri.stimulus(scans, times, duration, tr, times = TRUE)
x.group <- fmri.designG(hrf, subj = subj, runs = runs)
# View(x.group)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ design }
