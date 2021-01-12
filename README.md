
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fmri

<!-- badges: start -->

<!-- badges: end -->

The goal of fmri is to perform an fMRI analysis as described in Tabelow
et al. (2006) <https://doi.org/10.1016/j.neuroimage.2006.06.029>,
Polzehl et al. (2010)
<https://doi.org/10.1016/j.neuroimage.2010.04.241>, Tabelow and Polzehl
(2011) <https://doi.org/10.18637/jss.v044.i11>.

## Installation

You can install the released version of fmri from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fmri")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("muschellij2/fmri")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fmri)
## basic example code
```

``` r
gkernsm <- function(y,h=1) {
  grid <- function(d) {
    d0 <- d%/%2+1
    gd <- seq(0,1,length=d0)
    if (2*d0==d+1) gd <- c(gd,-gd[d0:2]) else gd <- c(gd,-gd[(d0-1):2])
    gd
  }
  dy <- dim(y)
  if (is.null(dy)) dy<-length(y)
  ldy <- length(dy)
  if (length(h)!=ldy) h <- rep(h[1],ldy)
  kern <- switch(ldy,dnorm(grid(dy),0,2*h/dy),
                 outer(dnorm(grid(dy[1]),0,2*h[1]/dy[1]),
                       dnorm(grid(dy[2]),0,2*h[2]/dy[2]),"*"),
                 outer(outer(dnorm(grid(dy[1]),0,2*h[1]/dy[1]),
                             dnorm(grid(dy[2]),0,2*h[2]/dy[2]),"*"),
                       dnorm(grid(dy[3]),0,2*h[3]/dy[3]),"*"))
  kern <- kern/sum(kern)
  kernsq <- sum(kern^2)
  list(gkernsm = convolve(y,kern,conj=TRUE),kernsq=kernsq)
}

create.mask <- function(){
  mask <- array(0,dim=c(65,65,26))
  mask[5:10,5:10,] <- 1
  mask[7:8,7:8,] <- 0
  mask[8:10,8:10,] <- 0
  mask[14:17,14:17,] <- 1
  mask[16:17,16:17,] <- 0
  mask[21:23,21:23,] <- 1
  mask[22:23,23,] <- 0
  mask[23,22,] <- 0
  mask[27:28,27:28,] <- 1
  mask[28,28,] <- 0
  mask[5:7,29:33,] <- 1
  mask[7,32:33,] <- 0
  mask[14:15,30:33,] <- 1
  mask[15,30,] <- 0
  mask[21,31:33,] <- 1
  mask[22,33,] <- 1
  mask[27,32:33,] <- 1
  mask[29:33,5:7,] <- 1
  mask[32:33,7,] <- 0
  mask[30:33,14:15,] <- 1
  mask[30,15,] <- 0
  mask[31:33,21,] <- 1
  mask[33,22,] <- 1
  mask[32:33,27,] <- 1
  mask[34:65,1:33,] <- mask[32:1,1:33,]
  mask[1:33,34:65,] <- mask[1:33,32:1,]
  mask[34:65,34:65,] <- mask[32:1,32:1,]
  mask
}

create.sig <- function(signal=1.5,efactor=1.2){
  sig <- array(0,dim=c(65,65,26))
  sig[29:37,38:65,] <- signal
  sig[38:65,38:65,] <- signal * efactor
  sig[38:65,29:37,] <- signal * efactor^2
  sig[38:65,1:28,] <- signal * efactor^3
  sig[29:37,1:28,] <- signal * efactor^4
  sig[1:28,1:28,] <- signal * efactor^5
  sig[1:28,29:37,] <- signal * efactor^6
  sig[1:28,38:65,] <- signal * efactor^7
  sig * create.mask()
}
# some values describing the data
signal <- 1.5
noise <- 20
arfactor <- .3

# maximaum bandwidth for adaptive smoothing
hmax <- 3.06

# datacube dimension 
i <- 65
j <- 65
k <- 26
scans <- 107

# define needed arrays
ttt <- array(0,dim=c(i,j,k,scans))
sig <- array(0,dim=c(i,j,k))

# create the mask for activation
mask <- create.mask()

# assign amplitudes of signals to activated areas 
sig <- create.sig(signal)
```

``` r
# expected BOLD response for some stimulus
hrf <- signal * fmri.stimulus(scans, c(18, 48, 78), 15, 2)

# create time series
dim(sig) <- c(i*j*k,1)
dim(hrf) <- c(1,scans)
sig4 <- sig %*% hrf
dim(sig) <- c(i,j,k)
dim(sig4) <- c(i,j,k,scans)

RNGversion("3.5.0")
#> Warning in RNGkind("Mersenne-Twister", "Inversion", "Rounding"): non-uniform
#> 'Rounding' sampler used
# create noise with spatial and temporal correlation
set.seed(1)
noisy4 <- rnorm(i*j*k*scans,0,noise)
dim(noisy4) <- c(i,j,k,scans)
for (t in 2:scans) noisy4[,,,t] <- noisy4[,,,t] + arfactor*noisy4[,,,t-1]
for (t in 1:scans) noisy4[,,,t] <- gkernsm(noisy4[,,,t],c(0.8,0.8,0.4))$gkernsm

# finally we got the data
ttt <- sig4 + noisy4
data1 <- list(ttt=writeBin(as.numeric(ttt),raw(),4),
              dim=c(i,j,k,scans),weights=c(1,1,2),mask=array(1,c(i,j,k)),
              delta = rep(1, 4))
class(data1) <- "fmridata"
```

``` r
# create design matrix and estimate parameters from linear model
hrf <- fmri.stimulus(scans, c(18, 48, 78), 15, 2)
z <- fmri.design(hrf)
spm <- fmri.lm(data1,z)
```

``` r
# adaptively smooth the spm
resultaws <- fmri.smooth(spm,hmax=hmax,lkern="Gaussian")
detectaws <- fmri.pvalue(resultaws)
pmask <- apply(detectaws$pvalue<0.05,c(1,2),sum)

# smooth non adaptively
resultnonaws <- fmri.smooth(spm,hmax=hmax,adaptation="none",lkern="Gaussian")
detectnonaws <- fmri.pvalue(resultnonaws)
npmask <- apply(detectnonaws$pvalue<0.05,c(1,2),sum)

# at last show some nice images
par(mfrow=c(1,3),mar=c(2.5,2.5,3,0.5))
image(1:i,1:j,sig[,,1],col=grey(0:255/255),xlab="",ylab="",main="True activation")
image(1:i,1:j,pmask,xlab="",ylab="",main="Detected using AWS")
image(1:i,1:j,npmask,xlab="",ylab="",main="Detected using Gaussian filter")
```
