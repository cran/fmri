library(fmri)
source("../R/randomfield.R")

xd <- 64
yd <- 64
zd <- 64
scans <- 107
hrf <- fmri.stimulus(scans,c(18,48,78),15,2)
x <- fmri.design(hrf)
signal <- 4
noise <- 18
factor <- 1.2

noise4d <- array(rnorm(xd*yd*zd*scans,0,noise),c(xd,yd,zd,scans))
signal3d <- array(0,c(xd,yd,zd))

signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),8] <- signal * factor^0
signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),18] <- signal * factor^1
signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),28] <- signal * factor^2
signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),38] <- signal * factor^3
signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),48] <- signal * factor^4
signal3d[c(8,18,28,38,48,58),c(8,18,28,38,48,58),58] <- signal * factor^5

dim(signal3d) <- c(xd*yd*zd,1)
dim(hrf) <- c(1,scans)
signal4d <- signal3d %*% hrf
dim(signal3d) <- c(xd,yd,zd)
dim(signal4d) <- c(xd,yd,zd,scans)

data <- list()
data$ttt <- signal4d + noise4d
data$dim <- c(xd,yd,zd,scans)
data$mask <- array(1,c(xd,yd,zd))
class(data) <- "fmridata"
data$weigths <- c(1,1,1)
data$format <- "ARTIFICIAL"

spm <- fmri.lm(data,x,actype="nocalc")
spmsnonaws <- fmri.smooth(spm,hmax=3.52,adaptive=FALSE)
spmsaws <- fmri.smooth(spm,hmax=3.52)

stat <- spm$cbeta / sqrt(spm$var)
statnonaws <- spmsnonaws$cbeta / sqrt(spmsnonaws$var)
stataws <- spmsaws$cbeta / sqrt(spmsaws$var)

pvaluevoxel <- fmri.pvalue(spm)
pvaluenonaws <- fmri.pvalue(spmsnonaws)
pvalueaws <- fmri.pvalue(spmsaws)

par(mfrow=c(1,3))
image(stat[,,8]);
image(statnonaws[,,8]);
image(stataws[,,8]);
answer <- readline("Return")

image(stat[,,18]);
image(statnonaws[,,18]);
image(stataws[,,18]);
answer <- readline("Return")

image(stat[,,28]);
image(statnonaws[,,28]);
image(stataws[,,28]);
answer <- readline("Return")

image(stat[,,38]);
image(statnonaws[,,38]);
image(stataws[,,38]);
answer <- readline("Return")

image(stat[,,48]);
image(statnonaws[,,48]);
image(stataws[,,48]);
answer <- readline("Return")

image(stat[,,58]);
image(statnonaws[,,58]);
image(stataws[,,58]);
answer <- readline("Return")
