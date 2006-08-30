library(fmri)
source("../R/randomfield.R")

xd <- 32
yd <- 32
zd <- 32
scans <- 107

data <- list()
data$ttt <- array(rnorm(xd*yd*zd*scans),c(xd,yd,zd,scans))

cat("... created data cube")

data$dim <- c(xd,yd,zd,scans)
data$mask <- array(1,c(xd,yd,zd))
class(data) <- "fmridata"
data$weigths <- c(1,1,1)
data$format <- "ARTIFICIAL"

cat("... with appropriate list entries\n")

hrf <- fmri.stimulus(scans,c(18,48,78),15,2)
x <- fmri.design(hrf)

cat("... estimate parameters of the linear model\n")

spm <- fmri.lm(data,x,actype="nocalc")

spmsmoothnonaws <- fmri.smooth(spm,hmax=3.52,adaptive=FALSE)
spmsmoothaws <- fmri.smooth(spm,hmax=3.52)
  
statnonaws <- spmsmoothnonaws$cbeta / sqrt(spmsmoothnonaws$var)
stataws <- spmsmoothaws$cbeta / sqrt(spmsmoothaws$var)

rxyzna0 <- c(median(spmsmoothnonaws$rxyz0[,1]),median(spmsmoothnonaws$rxyz0[,2]),median(spmsmoothnonaws$rxyz0[,3]))
fthreshna05 <- threshold(0.05,spm$dim[1],spm$dim[2],spm$dim[3],rxyzna0[1],rxyzna0[2],rxyzna0[3],type="norm") # eigentlich t!!!!
fthreshna01 <- threshold(0.01,spm$dim[1],spm$dim[2],spm$dim[3],rxyzna0[1],rxyzna0[2],rxyzna0[3],type="norm")
fthreshna001 <- threshold(0.001,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm")
rxyza0 <- c(median(spmsmoothaws$rxyz0[,1]),median(spmsmoothaws$rxyz0[,2]),median(spmsmoothaws$rxyz0[,3]))
fthresha05 <- threshold(0.05,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm") # eigentlich t!!!!
fthresha01 <- threshold(0.01,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm")
fthresha001 <- threshold(0.001,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm")


c(fthreshna05,fthreshna01,fthreshna001)
c(fthresha05,fthresha01,fthresha001)

pvna <- fmri.pvalue(spmsmoothnonaws)
pva <- fmri.pvalue(spmsmoothaws)
