\name{fmri.metaPar}
\alias{fmri.metaPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Mixed-effects Meta-Analysis model for fMRI data }
\description{
Group maps are estimated from BOLD effect estimates and their variances previously determined for each subject. The function \code{\link{rma.uni}} from R package \pkg{metafor} is used to fit mixed-effects meta-analytic models at group level. Voxel-wise regression analysis is accelerated by optional parallel processing using R package \pkg{parallel}.
}
\usage{
fmri.metaPar(Cbold, Vbold, XG = NULL, model = NULL, method = "REML",
             weighted = TRUE, knha = FALSE, mask = NULL, cluster = 2,
             wghts = c(1, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Cbold}{ a 4D-Array with the aggregated individual BOLD contrast estimates in standard space, e.g. all \code{cbeta} maps obtained from single-session analysis with \code{\link{fmri.lm}} may put together. Dimensions 1 to 3 define the voxel space, dimension 4 indicates a subject. If not the whole brain but a region is analyzed, vectors with region-indices can be preserved by adding as attributes (e.g. \code{attr(Cbold, "xind") <- xind}).
}
  \item{Vbold}{ a 4D-Array with the aggregated variance estimates for the contrast parameters in \code{Cbold}, e.g. all \code{var} maps obtained from single-session analysis with \code{\link{fmri.lm}} may put together. Dimensions 1 to 3 define the voxel space, dimension 4 indicates a subject.
}
  \item{XG}{ optionally, a group-level design matrix of class \code{"data.frame"} to include one or more moderators in the model. By default, an intercept is added to the model.
}
  \item{model}{ optionally, a one-sided formula of the form: \code{model <- ~ mod1 + mod2 + mod3} describing a model with moderator variables. Adding \code{"-1"} removes the intercept term.
}
  \item{method}{ a character string specifying whether a fixed- (method = "FE") or a random/mixed-effects model (method = "REML", default) should be fitted. Further estimators for random/mixed-effects models are available, see documentation of \code{\link{rma.uni}} function for more details.
}
  \item{weighted}{ logical indicating whether weighted (\code{weighted = TRUE}, default) or unweighted estimation should be used to fit the model.
}
  \item{knha}{ logical specifying whether the method by Knapp and Hartung (2003) should be used for adjusting standard errors of the estimated coefficients (default is FALSE). The Knapp and Hartung adjustment is only meant to be used in the context of random- or mixed-effects models.
}
  \item{mask}{ if available, a logical 3D-Array of dimensionality of the data (without 4th subject component) describing a brain mask. The computation is restricted to the selected voxels.
}
  \item{cluster}{ number of threads for parallel processing, which is limited to available multi-core CPUs. If you do not know your CPUs, try: \code{detectCores()} from \pkg{parallel} package. Presets are 2 threads. \code{cluster = 1} does not use \pkg{parallel} package.
}
  \item{wghts}{ a vector of length 3 specifying ratio of voxel dimensions. Isotropic voxels (e.g. MNI-space) are set as default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{fmri.metaPar()} fits the configured linear mixed-effects meta-analytic (MEMA) model separately at each voxel and extracts the first regression coefficient (usually the overall group mean), corresponding squared standard errors and degrees of freedom as well as the residuals from resulting \code{rma.uni()} objects, to obtain a statistical parametric map (SPM) for the group. Voxel-by-voxel analysis is performed by either the function \code{\link{apply}} or \code{\link{parApply}} from \pkg{parallel} package, which walks through the \code{Cbold} array.

This two-stage approach reduces the computational burden of fitting a full linear mixed-effects (LME) model, \code{\link{fmri.lmePar}} would do. It assumes first level design is same across subjects and normally distributed not necessarily homogeneous within-subject errors. Warping to standard space has been done before first-stage analyses are carried out. Either no masking or a uniform brain mask should be applied at individual subject analysis level, to avoid loss of information at group level along the edges.

At the second stage, observed individual BOLD effects from each study are combined in a meta-analytic model. There is the opportunity of weighting the fMRI studies by the precision of their respective effect estimate to take account of first level residual heterogeneity (\code{weighted = TRUE}). This is how to deal with intra-subject variability. The REML estimate of cross-subject variability (tau-squared) assumes that each of these observations is drawn independently from the same Gaussian distribution. Since correlation structures cannot be modeled, multi-subject fMRI studies with repeated measures cannot be analyzed in this way.

Spatial correlation among voxels, e.g. through the activation of nearby voxels, is ignored at this stage, but corrects for it, when random field theory define a threshold for significant activation at inference stage.

It is recommended to check your model syntax and residuals choosing some distinct voxels before running the model in loop (see Example). Error handling default is to stop if one of the threads produces an error. When this occurs, the output will be lost from any voxel, where the model has fitted successfully.
}
\value{
An object of class \code{"fmrispm"} and \code{"fmridata"}, basically a list with components:
  \item{beta}{estimated regression coefficients}
  \item{se}{estimated standard errors of the coefficients}
  \item{cbeta}{estimated BOLD contrast parameters for the group. Always the first regression coefficient is taken.}
  \item{var}{estimated variance of the BOLD contrast parameters}
  \item{mask}{brain mask}
  \item{residuals}{raw (integer size 2) vector containing residuals of the estimated linear mixed-effects meta-analytic model up to scale factor \code{resscale}}
  \item{resscale}{\code{resscale*extractData(object,"residuals")} are the residuals.}
  \item{tau2}{estimated amount of (residual) heterogeneity. Always 0 when \code{method = "FE"}.}
  \item{rxyz}{array of smoothness from estimated correlation for each voxel in resel space (for analysis without smoothing).}
  \item{scorr}{array of spatial correlations with maximal lags 5, 5, 3 in x, y and z-direction}
  \item{bw}{vector of bandwidths (in FWHM) corresponding to the spatial correlation within the data}
  \item{weights}{ratio of voxel dimensions}
  \item{dim}{dimension of the data cube and residuals}
  \item{df}{degrees of freedom for t-statistics, df = (n-p-1)}
  \item{sessions}{number of observations entering the meta-analytic model, n}
  \item{coef}{number of coefficients in the meta-analytic model (including the intercept, p+1)}
  \item{method}{estimator used to fit the meta-analytic model. In case of "FE", it is weighted or unweighted least squares.}
  \item{weighted}{estimation with inverse-variance weights}
  \item{knha}{Knapp and Hartung adjustment}
  \item{model}{meta-analytic regression model}
  \item{cluster}{number of threads running in parallel}
  \item{attr(*,"design")}{group-level design matrix}
  \item{attr(*,"approach")}{two-stage estimation method}
}
\references{
Chen G., Saad Z.S., Nath A.R., Beauchamp M.S., Cox R.W. (2012). FMRI group analysis combining effect estimates and their variances. NeuroImage, 60: 747-765.

Knapp G. and Hartung J. (2003). Improved tests for a random effects meta-regression with a single covariate. Statistics in Medicine, 22: 2693-2710.

Viechtbauer W. (2005). Bias and efficiency of meta-analytic variance estimators in the random-effects model. Journal of Educational and Behavioral Statistics, 30: 261-293.

Viechtbauer W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3): 1-48

Viechtbauer W. (2015). \emph{metafor: Meta-Analysis Package for R} R package version 1.9-7.
}
\author{ Sibylle Dames }
\note{
Meta analyses tend to be less powerful for neuroimaging studies, because they only have as many degrees of freedom as number of subjects. If the number of subjects is very small, then it may be impossible to estimate the between-subject variance (tau-squared) with any precision. In this case the fixed effect model may be the only viable option. However, there is also the possibility of using a one-stage model, that includes the full time series data from all subjects and simultaneously estimates subject and group levels parameters (see \code{\link{fmri.lmePar}}). Although this approach is much more computer intensive, it has the advantage of higher degrees of freedom (> 100) at the end.

\itemize{Current Limitations\cr
The function cannot handle:
\item{experimental designs with a within-subject (repeated measures) factor}
\item{paired samples with varying tasks, unless the contrast of the two conditions is used as input}
}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{rma.uni}}, \code{\link{fmri.lm}}, \code{\link{fmri.lmePar}}
}
\examples{
\dontrun{## Generate some fMRI data sets: noise + stimulus
dx <- dy <- dz <- 32
dt <- 107
hrf <- fmri.stimulus(dt, c(18, 48, 78), 15, 2)
stim <- matrix(hrf, nrow= dx*dy*dz, ncol=dt, byrow=TRUE)
mask <- array(FALSE, c(dx, dy, dz))
mask[12:22,12:22,12:22] <- TRUE

ds1 <- list(ttt=writeBin(1.0*rnorm(dx*dy*dz*dt) + as.vector(5*stim),
           raw(), 4), mask = mask, dim = c(dx, dy, dz, dt))
ds2 <- list(ttt=writeBin(1.7*rnorm(dx*dy*dz*dt) + as.vector(3*stim),
           raw(), 4), mask = mask, dim = c(dx, dy, dz, dt))
ds3 <- list(ttt=writeBin(0.8*rnorm(dx*dy*dz*dt) + as.vector(1*stim),
           raw(), 4), mask = mask, dim = c(dx, dy, dz, dt))
ds4 <- list(ttt=writeBin(1.2*rnorm(dx*dy*dz*dt) + as.vector(2*stim),
           raw(), 4), mask = mask, dim = c(dx, dy, dz, dt))
class(ds1) <- class(ds2) <- class(ds3) <- class(ds4) <- "fmridata"

## Stage 1: single-session regression analysis
x <- fmri.design(hrf, order=2)
spm.sub01 <- fmri.lm(ds1, x, mask, actype = "smooth", verbose = TRUE)
spm.sub02 <- fmri.lm(ds2, x, mask, actype = "smooth", verbose = TRUE)
spm.sub03 <- fmri.lm(ds3, x, mask, actype = "smooth", verbose = TRUE)
spm.sub04 <- fmri.lm(ds4, x, mask, actype = "smooth", verbose = TRUE)

## Store observed individual BOLD effects and their variance estimates
subj <- 4
Cbold <- array(0, dim = c(dx, dy, dz, subj))
Cbold[,,,1] <- spm.sub01$cbeta
Cbold[,,,2] <- spm.sub02$cbeta
Cbold[,,,3] <- spm.sub03$cbeta
Cbold[,,,4] <- spm.sub04$cbeta

Vbold <- array(0, dim = c(dx, dy, dz, subj))
Vbold[,,,1] <- spm.sub01$var
Vbold[,,,2] <- spm.sub02$var
Vbold[,,,3] <- spm.sub03$var
Vbold[,,,4] <- spm.sub04$var

## Stage 2: Random-effects meta-regression analysis
## a) Check your model
library(metafor)
M1.1 <- rma.uni(Cbold[16,16,16, ],
                Vbold[16,16,16, ],
                method = "REML",
                weighted = TRUE,
                knha = TRUE,
                verbose = TRUE,
                control = list(stepadj=0.5, maxiter=2000, threshold=0.001))

# Control list contains convergence parameters later used
# at whole data cube. Values were adjusted to fMRI data.

summary(M1.1)
forest(M1.1)
qqnorm(M1.1)

## b) Estimate a group map
## without parallelizing
spm.group1a <- fmri.metaPar(Cbold, Vbold, knha = TRUE,
                            mask = mask, cluster = 1)
## same with 4 parallel threads
spm.group1b <- fmri.metaPar(Cbold, Vbold, knha = TRUE,
                            mask = mask, cluster = 4)}
}
\keyword{ regression }% at least one, from doc/KEYWORDS
