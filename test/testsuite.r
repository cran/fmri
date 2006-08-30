library(fmri)
source("../R/randomfield.R")

cat("Starting testsuite for fmri package ...\n")

xd <- 32
yd <- 32
zd <- 32
scans <- 107

thresh05 <- qt(1-0.05,scans-4)
bonferroni05 <- qt(1-0.05/(xd*yd*yd),scans-4)
thresh01 <- qt(1-0.01,scans-4)
bonferroni01 <- qt(1-0.01/(xd*yd*yd),scans-4)
thresh001 <- qt(1-0.001,scans-4)
bonferroni001 <- qt(1-0.001/(xd*yd*yd),scans-4)


data <- list()

voxeltotal <- 0
falseposvoxel05 <- 0
falseposbonf05 <- 0
falseposfmri05 <- 0
falseposrft05 <- 0
falseposvoxel01 <- 0
falseposbonf01 <- 0
falseposfmri01 <- 0
falseposrft01 <- 0
falseposvoxel001 <- 0
falseposbonf001 <- 0
falseposfmri001 <- 0
falseposrft001 <- 0
count <- 0

while (TRUE) {
  data$ttt <- array(rnorm(xd*yd*zd*scans),c(xd,yd,zd,scans))

#  cat("... created data cube")

  data$dim <- c(xd,yd,zd,scans)
  data$mask <- array(1,c(xd,yd,zd))
  class(data) <- "fmridata"
  data$weigths <- c(1,1,1)
  data$format <- "ARTIFICIAL"

#  cat("... with appropriate list entries\n")
  
  hrf <- fmri.stimulus(scans,c(18,48,78),15,2)
  x <- fmri.design(hrf)

#  cat("... no estimate parameters of the linear model\n")

  #spm <- fmri.lm(data,x)
  spm <- fmri.lm(data,x,actype="nocalc")

#  cat("---> estimated spatial correlation",spm$scorr,"\n")

#  cat("... create t-field\n")

  stat <- spm$cbeta / sqrt(spm$var)

  voxeltotal <- voxeltotal + xd*yd*zd
  count <- count + 1
  falseposvoxel05 <- falseposvoxel05 + sum(stat>thresh05)
  falseposvoxel01 <- falseposvoxel01 + sum(stat>thresh01)
  falseposvoxel001 <- falseposvoxel001 + sum(stat>thresh001)
  if (any(stat>bonferroni05)) falseposbonf05 <- falseposbonf05 + 1
  if (any(stat>bonferroni01)) falseposbonf01 <- falseposbonf01 + 1
  if (any(stat>bonferroni001)) falseposbonf001 <- falseposbonf001 + 1

  rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
  fthresh05 <- threshold(0.05,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm") # eigentlich t!!!!
  fthresh01 <- threshold(0.01,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm")
  fthresh001 <- threshold(0.001,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm")
  if (any(stat>fthresh05)) falseposrft05 <- falseposrft05 + 1
  if (any(stat>fthresh01)) falseposrft01 <- falseposrft01 + 1
  if (any(stat>fthresh001)) falseposrft001 <- falseposrft001 + 1

  pv <- fmri.pvalue(spm)
  if (any(pv$pvalue < 0.05)) falseposfmri05 <- falseposfmri05 + 1
  if (any(pv$pvalue < 0.01)) falseposfmri01 <- falseposfmri01 + 1
  if (any(pv$pvalue < 0.001)) falseposfmri001 <- falseposfmri001 + 1
  
  cat("RESULTS for",count,"run\n")
  cat("  thresholds:",thresh05,thresh01,thresh001,bonferroni05,bonferroni01,bonferroni001,"\n")
  cat("  thresholds:",fthresh05,fthresh01,fthresh001,"\n")
  cat("  t-stat    :",range(stat),"\n")
  cat("---> false positives voxelwise  at p=0.05 :",falseposvoxel05/voxeltotal,"\n")
  cat("---> false positives voxelwise  at p=0.01 :",falseposvoxel01/voxeltotal,"\n")
  cat("---> false positives voxelwise  at p=0.001:",falseposvoxel001/voxeltotal,"\n")
  cat("---> false positives bonferroni at p=0.05 :",falseposbonf05/count,"\n")
  cat("---> false positives bonferroni at p=0.01 :",falseposbonf01/count,"\n")
  cat("---> false positives bonferroni at p=0.001:",falseposbonf001/count,"\n")
  cat("---> false positives rft at p=0.05        :",falseposrft05/count,"\n")
  cat("---> false positives rft at p=0.01        :",falseposrft01/count,"\n")
  cat("---> false positives rft at p=0.001       :",falseposrft001/count,"\n")
  cat("---> false positives fmri at p=0.05       :",falseposfmri05/count,"\n")
  cat("---> false positives fmri at p=0.01       :",falseposfmri01/count,"\n")
  cat("---> false positives fmri at p=0.001      :",falseposfmri001/count,"\n")
  
}
  

##############################################

#arfactor <- 0.4
#for (t in 2:scans) data$ttt[,,,t] <- data$ttt[,,,t] +  arfactor*data$ttt[,,,t-1]
