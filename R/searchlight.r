searchlight <- function(radius){
   rad <- as.integer(radius)
   nr <- 2*rad+1
   indices <- rbind(rep(-rad:rad,nr*nr),
              rep(-rad:rad,rep(nr,nr)),
              rep(-rad:rad,rep(nr*nr,nr)))
   indices[,apply(indices^2,2,sum)<=radius^2]
}

searchlightdistr.old <- function(df,nregion,kind=0,nsim=500000){
old.seed <- .Random.seed
set.seed(10:100)
rn <- matrix(rt(nsim*nregion,df),nregion,nsim)
.Random.seed <- old.seed
rsl <- apply(if(kind==0) abs(rn) else rn^2, 2, mean)
list(p=seq(1/(2*nsim),(2*nsim-1)/(2*nsim),length=nsim),kvalue=sort(rsl),df=df,kind=kind,nsim=nsim)
}

searchlightdistr <- function(spm,radius,kind=0,sdim=c(200,200,200)){
old.seed <- .Random.seed
set.seed(10:100)
df <- spm$df
bw <- spm$bw
nsim <- prod(sdim)
rn <- array(rt(nsim,df),sdim)
if(any(bw>0)){
# emulate correlation structure in spm
    sdrn <- sd(rn)
    rn <- aws::kernsm(rn, h=bw, unit="FWHM")@yhat
    rn <- rn*sdrn/sd(rn)
}
rn <- if(kind==0) abs(rn) else rn^2
sregion <- searchlight(radius)
nregion <- dim(sregion)[2]
mask <- array(TRUE,sdim)
rsl <- .Fortran(C_slight,
                as.double(rn),
                as.integer(mask),
                as.integer(sdim[1]),
                as.integer(sdim[2]),
                as.integer(sdim[3]),
                as.integer(sregion),
                as.integer(nregion),
                stat=double(nsim))$stat
list(p=seq(1/(2*nsim),(2*nsim-1)/(2*nsim),length=nsim),kvalue=sort(rsl),df=df,kind=kind,nsim=nsim)
}

fmri.searchlight <- function(spm, alpha=.05, radius, minimum.signal=0, kind=c("abs","squared")){
  args <- sys.call()
  args <- c(spm$call,args)

cat("fmri.searchlight: entering function\n")

if (!inherits(spm,"fmrispm")) {
  warning("fmri.searchlight: data not of class <fmrispm>. Try to proceed but strange things may happen")
}
stat <- spm$cbeta
stat[stat>0] <- pmax(0,stat[stat>0]-minimum.signal)
stat[stat<0] <- pmin(0,stat[stat<0]+minimum.signal)
stat <- stat/sqrt(spm$var)
mask <- spm$mask
stat[!mask] <- 0
dim(stat) <- dimspm <- spm$dim[1:3]
sregion <- searchlight(radius)
nregion <- dim(sregion)[2]
cat("fmri.searchlight: using",kind,"size of searchlight:",nregion,"\n")
kind <- if(kind[1]=="abs") 0 else 1
stat <- if(kind==0) abs(stat) else stat^2
##
##  get statistics over searchlights
##
stat <- .Fortran(C_slight,
                as.double(stat),
                as.integer(mask),
                as.integer(dimspm[1]),
                as.integer(dimspm[2]),
                as.integer(dimspm[3]),
                as.integer(sregion),
                as.integer(nregion),
                stat=double(prod(dimspm)))$stat
##
##  approximate distribution of searchlight statistics
##  assumes t-distributed spm
##
cat("fmri.searchlight: get empirical distribution of sl-statistic \n")
# distr <- searchlightdistr.old(spm$df,nregion,kind)
distr <- searchlightdistr(spm,radius,kind)
nvoxel <- sum(mask)
cat("fmri.searchlight: get approximate pvalues \n")
pv <- array(1,dimspm)
pv[mask] <- .Fortran(C_getslpv,
               as.double(stat[mask]),
               as.integer(nvoxel),
               as.double(distr$p),
               as.double(distr$kvalue),
               as.integer(distr$nsim),
               pvalue = double(nvoxel))$pvalue
ind <- fdr(pv[mask], alpha)
thresh <- min(stat[mask][ind])
mask[stat < thresh] <- FALSE
    mask <- mask & spm$mask
    pv[!mask] <- NA
    z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim,
            hrf = spm$hrf, mask = spm$mask)
        class(z) <- c("fmripvalue")
        z$roixa <- spm$roixa
        z$roixe <- spm$roixe
        z$roiya <- spm$roiya
        z$roiye <- spm$roiye
        z$roiza <- spm$roiza
        z$roize <- spm$roize
        z$roit <- spm$roit
        z$header <- spm$header
        z$format <- spm$format
        z$call <- args
        attr(z, "file") <- attr(spm, "file")
        attr(z, "white") <- attr(spm, "white")
        attr(z, "design") <- attr(spm, "design")
        attr(z, "smooth") <- "Not smoothed"
        attr(z, "mode") <- paste("Threshold: searchlight radius=",radius, "\n")
        z
    }

getSearchlightPattern <- function(spm, voxel, radius){
   if(!inherits(spm,"fmrispm")) stop("getSearchlightPattern: first argument needs to be an SPM")
   if(!is.logical(voxel)) stop("getSearchlightPattern: second argument needs to be a voxel mask")
   if(any(dim(voxel)!=spm$dim[1:3])) stop("getSearchlightPattern: Incompatible dimensions")
   if(!is.numeric(radius) ||radius < 1) stop("getSearchlightPattern: illegal radius")
   beta <- spm$beta
   ddim <- dim(beta)
   sregion <- searchlight(radius)
   nsl <- dim(sregion)[2]
   nb <- ddim[4]
   if(any(is.na(voxel))) voxel[is.na(voxel)] <- FALSE
   nvox <- sum(voxel)
   z <- .Fortran(C_extrpatt,
                 as.double(spm$beta),
                 as.integer(voxel),
                 as.integer(ddim[1]),
                 as.integer(ddim[2]),
                 as.integer(ddim[3]),
                 as.integer(ddim[4]),
                 as.integer(sregion),
                 as.integer(nsl),
         pattern=double(nb*nsl*nvox),
                 as.integer(nvox))$pattern
   invisible(array(z,c(nb,nsl,nvox)))
}
