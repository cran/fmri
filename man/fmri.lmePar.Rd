\name{fmri.lmePar}
\alias{fmri.lmePar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Linear Mixed-effects Model for fMRI data }
\description{
Group maps are directly estimated from the BOLD time series data of all subjects using \code{\link{lme}} from R package \pkg{nlme} to fit a Linear Mixed-effects Model with temporally correlated and heteroscedastic within-subject errors. Voxel-wise regression analysis is accelerated by optional parallel processing using R package \pkg{parallel}.
}
\usage{
fmri.lmePar(bold, z, fixed = NULL, random = NULL, mask = NULL, 
            ac = 0.3, vtype = "individual", cluster = 2, 
            wghts = c(1, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bold}{ a large 4D-Array with the aggregated fMRI data of all subjects that were previously registered to a common brain atlas. Be careful with the assembly of this array, the order of the data sets has to be compatible with the design matrix: \code{"z"}. If not the whole brain but a region is analyzed, vectors with region-indices can be preserved by adding as attributes (e.g. \code{attr(bold, "xind") <- xind}).
  }  
  \item{z}{ a design matrix for a multi-subject and/or multi-session fMRI-study of class \code{"data.frame"} specifying the expected BOLD response(s) and additional components for trend and other effects. Typically a \code{\link{fmri.designG}} object. This data frame contains all variables named in the model. There are some indispensable variables: \var{"group"}, \var{"subj"}, \var{"session"} and \var{"run"}, which define the different strata. That information will be used for setting up the residual variance structure.
  }  
  \item{fixed}{ optionally, a one-sided linear formula describing the fixed-effects part of the model. Default settings are: 
  \code{fixed <- ~ 0 + hrf + session + drift1:session + drift2:session} in case of one detected group, and the same but \code{"hrf"} replaced with \code{"hrf:group"} if two group levels in \code{z} are found. Since an intercept would be a linear combination of the session factor-variable modeling session-specific intercepts, it is excluded. 
  }  
  \item{random}{ optionally, a one-sided formula of the form \code{~ x1 + ... + xn | g1/.../gm}, with \code{~ x1 + ... + xn} specifying the model for the random effects and \code{g1/.../gm} the grouping structure.\cr 
  
  Default is always the basic model without covariates, i.e. \cr 
  \code{random <- ~ 0 + hrf|subj}  if no repeated measures in \code{z} are found (\code{nlevels(z$run)==1}),\cr 
  \code{random <- ~ 0 + hrf|subj/session}  if repeated measures and\cr 
  \code{random <- ~ 0 + hrf|session}  if repeated measures but one subject only. \cr 
  In case of two independent groups:\cr 
  \code{random <- list(subj = pdDiag(~ 0 + hrf:group))} is used.
  }  
  \item{mask}{ if available, a logical 3D-Array of dimensionality of the data (without time component) describing a brain mask. The computation is restricted to the selected voxels.
  }  
  \item{ac}{ if available, a numeric 3D-Array of dimensionality of the data (without time component) with spatially smoothed autocorrelation parameters should be used in the AR(1) models fitted in each voxel, e.g. locally estimated and smoothed AR(1)-coefficients from \code{\link{fmri.lm}} applied to the first subject. Alternatively, a global approach with uniform value can be used. In this case enter a number between 0 and 1. Default is \code{0.3} applied to all voxels.
  }  
  \item{vtype}{ a character string choosing the residual variance model. If \code{"equal"}, homoscedastic variance across subjects is assumed setting weights argument in function \code{lme()} to zero, whereas \code{"individual"} allows different within-subject variances. Default method is \code{"individual"} that means subject-specific error variances using formula: \code{weights <- varIdent(form =~ 1|subj)}.
  }  
  \item{cluster}{ number of threads for parallel processing, which is limited to available multi-core CPUs. If you do not know your CPUs, try: \code{detectCores()} from \pkg{parallel} package. Presets are 2 threads. \code{cluster = 1} does not use \pkg{parallel} package.
  }  
  \item{wghts}{ a vector of length 3 specifying ratio of voxel dimensions. Isotropic voxels (e.g. MNI-space) are set as default.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{fmri.lmePar()} fits the configured Linear Mixed-effects Model separately at each voxel and extracts estimated BOLD contrasts, corresponding squared standard errors and degrees of freedom as well as the residuals from resulting \code{lme()} objects to produce a statistical parametric map (SPM) for the group(s). Voxel-by-voxel analysis is performed by either the function \code{\link{apply}} or \code{\link{parApply}} from \pkg{parallel} package, which walks through the \code{bold} array.

If one group is analyzed, from each fitted model the first fixed-effects coefficient and corresponding parameters are stored in results object. This should be the first specified predictor in the fixed-effects part of the model (verify the attribute of \code{"df"} in returned object). However, in two-sample case this principle does not work. The order changes, estimated session-specific intercepts now comes first and the number of these coefficients is not fixed. Therefore in current version it has explicitly been looked for the coefficient names: \var{"hrf:group1"} and \var{"hrf:group2"}. Available functions within the \pkg{nlme} package to extract estimated values from \code{lme()} objects do not operate at contrast matrices.

Spatial correlation among voxels, e.g. through the activation of nearby voxels, is ignored at this stage, but corrects for it, when random field theory define a threshold for significant activation at inference stage.

It is recommended to check your model syntax and residuals choosing some distinct voxels before running the model in loop (see Example, step 1); especially for more advanced designs! Error handling default is to stop if one of the threads produces an error. When this occurs, the output will be lost from any voxel, where the model has fitted successfully.
}
\value{ 
An object of class \code{"fmrispm"} and \code{"fmridata"}, basically a list 
with components:
  \item{cbeta, cbeta2}{estimated BOLD contrast parameters separated for the groups 1 and 2}
  \item{var, var2}{estimated variance of the contrast parameters separated for the groups 1 and 2}
  \item{mask}{ brain mask }
  \item{res, res2}{raw (integer size 2) vector containing residuals of the estimated Linear Mixed-effects Model up to scale factor \code{resscale} separated for the groups 1 and 2}
  \item{resscale, resscale2}{\code{resscale*extract.data(object,"residuals")} are the residuals of group 1 and group 2 respectively. }
  \item{arfactor}{autocorrelation parameters used in AR(1)-model}
  \item{rxyz, rxyz2}{array of smoothness from estimated correlation for each voxel in resel space separated for the groups 1 and 2 (for analysis without smoothing) }
  \item{scorr, scorr2}{array of spatial correlations with maximal lags 5, 5, 3 in x, y and z-direction separated for the groups 1 and 2 }
  \item{bw, bw2}{vector of bandwidths (in FWHM) corresponding to the spatial correlation within the data separated for the groups 1 and 2}
  \item{weights}{ ratio of voxel dimensions }
  \item{dim, dim2}{dimension of the data cube and residuals separated for the groups 1 and 2}
  \item{df, df2}{degrees of freedom for t-statistics reported in \code{lme()} objects for the extracted regression coefficients separated for the groups 1 and 2. The name of the coefficient belonging to this df-value appears as attribute.}
  \item{subjects}{number of subjects in the study}
  \item{subj.runs}{number of repeated measures within subjects}
  \item{sessions}{number of total sessions that were analyzed}
  \item{groups}{number of groups in the study}
  \item{fixedModel}{fixed-effects model}
  \item{randomModel}{random-effects model}
  \item{VarModel}{assumption about the subject error variances}
  \item{cluster}{number of threads run in parallel}
  \item{attr(*,"design")}{design matrix for the multi-subject fMRI-study}
  \item{attr(*,"approach")}{one-stage estimation method}
}
\references{
Pinheiro J. and Bates D. (2000). Mixed-Effects Models in S and S-Plus. Springer.

Pinheiro J., Bates D., DebRoy S., Sarkar D. and the R Core team (2014). \emph{nlme: Linear and Nonlinear Mixed Effects Models} R package version 3.1-117.
}
\author{ Sibylle Dames }
\note{
Maybe the computing power is insufficient to carry out a whole brain analysis. You have two opportunities: either select and analyze a certain brain area or switch to a two-stage model.

\itemize{Current Limitations\cr
The function cannot handle experimental designs with:
\item{more than two independent groups}
\item{more than one stimulus (task)}
\item{paired samples with varying tasks}
\item{user defined contrasts} }

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lme}}, \code{\link{fmri.designG}}, 
  \code{\link{fmri.design}}, \code{\link{fmri.stimulus}}, 
  \code{\link{fmri.metaPar}}
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
            raw(), 4), mask=mask, dim=c(dx, dy, dz, dt))
ds2 <- list(ttt=writeBin(1.7*rnorm(dx*dy*dz*dt) + as.vector(3*stim), 
            raw(), 4), mask=mask, dim=c(dx, dy, dz, dt))
ds3 <- list(ttt=writeBin(0.8*rnorm(dx*dy*dz*dt) + as.vector(1*stim), 
            raw(), 4), mask=mask, dim=c(dx, dy, dz, dt))
ds4 <- list(ttt=writeBin(1.2*rnorm(dx*dy*dz*dt) + as.vector(2*stim), 
            raw(), 4), mask=mask, dim=c(dx, dy, dz, dt))
class(ds1) <- class(ds2) <- class(ds3) <- class(ds4) <- "fmridata"

## Construct a design matrix for a multi-subject study
subj <- 4
runs <- 1
z <-fmri.designG(hrf, subj = subj, runs = runs)

## Assembly of the aggregated BOLD-Array
Bold <- array(0, dim = c(dx,dy,dz,subj*runs*dt))
Bold[1:dx,1:dy,1:dz,1:(dt*1)] <- extract.data(ds1)
Bold[1:dx,1:dy,1:dz,(dt*1+1):(dt*2)] <- extract.data(ds2)
Bold[1:dx,1:dy,1:dz,(dt*2+1):(dt*3)] <- extract.data(ds3)
Bold[1:dx,1:dy,1:dz,(dt*3+1):(dt*4)] <- extract.data(ds4)

## Step 1: Check the model
y <- Bold[16, 16, 16, ] # choose one voxel
M1.1 <-  lme(fixed = y ~ 0 + hrf + session + drift1:session + drift2:session,
            random = ~ 0 + hrf|subj,           
            correlation = corAR1(value = 0.3, form = ~ 1|subj/session, fixed=TRUE),
            weights = varIdent(form =~ 1|subj), 
            method ="REML",
            control = lmeControl(rel.tol=1e-6, returnObject = TRUE),
            data = z)
summary(M1.1)

# Residual plots
plot(M1.1, resid(.,type = "response") ~ scan|subj)
qqnorm(M1.1, ~resid(.,type = "normalized")|subj, abline = c(0,1))

# Testing the assumption of homoscedasticity 
M1.2 <- update(M1.1, weights = NULL, data = z)
anova(M1.2, M1.1) 

# Model fit: observed and fitted values 
fitted.values <- fitted(M1.1)
plot(y[1:dt], type="l", main = "Subject 1", xlab = "scan", 
     ylab = "BOLD-signal", ylim = c(-5,5))
lines(fitted.values[names(fitted.values)==1],lty=1,lwd=2)

plot(y[(dt+1):(2*dt)], type="l", main = "Subject 2", xlab = "scan", 
     ylab = "BOLD-signal", ylim = c(-5,5))
lines(fitted.values[names(fitted.values)==2],lty=1,lwd=2)

plot(y[(2*dt+1):(3*dt)], type="l", main = "Subject 3", xlab = "scan", 
     ylab = "BOLD-signal", ylim = c(-5,5))
lines(fitted.values[names(fitted.values)==3],lty=1,lwd=2)

plot(y[(3*dt+1):(4*dt)], type="l", main = "Subject 4", xlab = "scan", 
     ylab = "BOLD-signal", ylim = c(-5,5))
lines(fitted.values[names(fitted.values)==4],lty=1,lwd=2)

## Step 2: Estimate a group map
## without parallelizing
spm.group1a <- fmri.lmePar(Bold, z, mask = mask, cluster = 1)
# same with 4 parallel threads
spm.group1b <- fmri.lmePar(Bold, z, mask = mask, cluster = 4)
## Example for two independent groups
group <- c(1,1,4,4)
z2 <- fmri.designG(hrf, subj = subj, runs = runs, group = group)
spm.group2 <- fmri.lmePar(Bold, z2, mask = mask, cluster = 4)}
}
\keyword{ regression }% at least one, from doc/KEYWORDS
