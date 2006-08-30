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
falseposvoxelna05 <- 0
falseposbonfna05 <- 0
falseposfmrina05 <- 0
falseposvoxelna01 <- 0
falseposbonfna01 <- 0
falseposfmrina01 <- 0
falseposvoxelna001 <- 0
falseposbonfna001 <- 0
falseposfmrina001 <- 0

falseposvoxela05 <- 0
falseposbonfa05 <- 0
falseposfmria05 <- 0
falseposvoxela01 <- 0
falseposbonfa01 <- 0
falseposfmria01 <- 0
falseposvoxela001 <- 0
falseposbonfa001 <- 0
falseposfmria001 <- 0
count <- 0

falsevoxelna05 <- numeric()
falsevoxelna01 <- numeric()
falsevoxelna001 <- numeric()
falsevoxela05 <- numeric()
falsevoxela01 <- numeric()
falsevoxela001 <- numeric()

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

  spmsmoothnonaws <- fmri.smooth(spm,hmax=3.52,adaptive=FALSE)
  spmsmoothaws <- fmri.smooth(spm,hmax=3.52)
  
#  cat("---> estimated spatial correlation",spm$scorr,"\n")

#  cat("... create t-field\n")

  statnonaws <- spmsmoothnonaws$cbeta / sqrt(spmsmoothnonaws$var)
  stataws <- spmsmoothaws$cbeta / sqrt(spmsmoothaws$var)

  voxeltotal <- voxeltotal + xd*yd*zd
  count <- count + 1
  falseposvoxelna05 <- falseposvoxelna05 + sum(statnonaws>thresh05)
  falseposvoxelna01 <- falseposvoxelna01 + sum(statnonaws>thresh01)
  falseposvoxelna001 <- falseposvoxelna001 + sum(statnonaws>thresh001)

  falseposvoxela05 <- falseposvoxela05 + sum(stataws>thresh05)
  falseposvoxela01 <- falseposvoxela01 + sum(stataws>thresh01)
  falseposvoxela001 <- falseposvoxela001 + sum(stataws>thresh001)

  rxyzna0 <- c(median(spmsmoothnonaws$rxyz0[,1]),median(spmsmoothnonaws$rxyz0[,2]),median(spmsmoothnonaws$rxyz0[,3]))
  fthreshna05 <- threshold(0.05,spm$dim[1],spm$dim[2],spm$dim[3],rxyzna0[1],rxyzna0[2],rxyzna0[3],type="norm") # eigentlich t!!!!
  fthreshna01 <- threshold(0.01,spm$dim[1],spm$dim[2],spm$dim[3],rxyzna0[1],rxyzna0[2],rxyzna0[3],type="norm")
  fthreshna001 <- threshold(0.001,spm$dim[1],spm$dim[2],spm$dim[3],rxyzna0[1],rxyzna0[2],rxyzna0[3],type="norm")
  rxyza0 <- c(median(spmsmoothaws$rxyz0[,1]),median(spmsmoothaws$rxyz0[,2]),median(spmsmoothaws$rxyz0[,3]))
  fthresha05 <- threshold(0.05,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm") # eigentlich t!!!!
  fthresha01 <- threshold(0.01,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm")
  fthresha001 <- threshold(0.001,spm$dim[1],spm$dim[2],spm$dim[3],rxyza0[1],rxyza0[2],rxyza0[3],type="norm")

  if (any(statnonaws>bonferroni05)) falseposbonfna05 <- falseposbonfna05 + 1
  if (any(statnonaws>bonferroni01)) falseposbonfna01 <- falseposbonfna01 + 1
  if (any(statnonaws>bonferroni001)) falseposbonfna001 <- falseposbonfna001 + 1

  if (any(stataws>bonferroni05)) falseposbonfa05 <- falseposbonfa05 + 1
  if (any(stataws>bonferroni01)) falseposbonfa01 <- falseposbonfa01 + 1
  if (any(stataws>bonferroni001)) falseposbonfa001 <- falseposbonfa001 + 1

  pvna <- fmri.pvalue(spmsmoothnonaws)
  if (any(pvna$pvalue < 0.05)) falseposfmrina05 <- falseposfmrina05 + 1
  if (any(pvna$pvalue < 0.01)) falseposfmrina01 <- falseposfmrina01 + 1
  if (any(pvna$pvalue < 0.001)) falseposfmrina001 <- falseposfmrina001 + 1
  
  falsevoxelna05[count] <- sum(pvna$pvalue < 0.05)
  falsevoxelna01[count] <- sum(pvna$pvalue < 0.01)
  falsevoxelna001[count] <- sum(pvna$pvalue < 0.001)

  pva <- fmri.pvalue(spmsmoothaws)
  if (any(pva$pvalue < 0.05)) falseposfmria05 <- falseposfmria05 + 1
  if (any(pva$pvalue < 0.01)) falseposfmria01 <- falseposfmria01 + 1
  if (any(pva$pvalue < 0.001)) falseposfmria001 <- falseposfmria001 + 1
  
  falsevoxela05[count] <- sum(pva$pvalue < 0.05)
  falsevoxela01[count] <- sum(pva$pvalue < 0.01)
  falsevoxela001[count] <- sum(pva$pvalue < 0.001)

  cat("RESULTS for",count,"run\n")
  cat("  thresholds:",thresh05,thresh01,thresh001,"\n")
  cat("---> false positives voxelwise nonaws  at p=0.05 :",falseposvoxelna05/voxeltotal,"\n")
  cat("---> false positives voxelwise nonaws  at p=0.01 :",falseposvoxelna01/voxeltotal,"\n")
  cat("---> false positives voxelwise nonaws  at p=0.001:",falseposvoxelna001/voxeltotal,"\n")
  cat("---> false positives voxelwise aws  at p=0.05 :",falseposvoxela05/voxeltotal,"\n")
  cat("---> false positives voxelwise aws  at p=0.01 :",falseposvoxela01/voxeltotal,"\n")
  cat("---> false positives voxelwise aws  at p=0.001:",falseposvoxela001/voxeltotal,"\n")
  cat("  thresholds:",bonferroni05,bonferroni01,bonferroni001,"\n")
  cat("---> false positives bonferroni nonaws at p=0.05 :",falseposbonfna05/count,"\n")
  cat("---> false positives bonferroni nonaws at p=0.01 :",falseposbonfna01/count,"\n")
  cat("---> false positives bonferroni nonaws at p=0.001:",falseposbonfna001/count,"\n")
  cat("---> false positives bonferroni aws at p=0.05 :",falseposbonfa05/count,"\n")
  cat("---> false positives bonferroni aws at p=0.01 :",falseposbonfa01/count,"\n")
  cat("---> false positives bonferroni aws at p=0.001:",falseposbonfa001/count,"\n")
  cat("  thresholds:",fthreshna05,fthreshna01,fthreshna001,"\n")
  cat("---> false positives fmri nonaws at p=0.05       :",falseposfmrina05/count,"\n")
  cat("---> false positives fmri nonaws at p=0.01       :",falseposfmrina01/count,"\n")
  cat("---> false positives fmri nonaws at p=0.001      :",falseposfmrina001/count,"\n")
  cat("  thresholds:",fthresha05,fthresha01,fthresha001,"\n")
  cat("---> false positives fmri aws at p=0.05       :",falseposfmria05/count,"\n")
  cat("---> false positives fmri aws at p=0.01       :",falseposfmria01/count,"\n")
  cat("---> false positives fmri aws at p=0.001      :",falseposfmria001/count,"\n")

  cat(falsevoxela05[count],falsevoxela01[count],falsevoxela001[count],"\n")
  cat(falsevoxelna05[count],falsevoxelna01[count],falsevoxelna001[count],"\n")

  if (length(falsevoxela05)>2) plot(density(falsevoxela05))
}
  

##############################################

#arfactor <- 0.4
#for (t in 2:scans) data$ttt[,,,t] <- data$ttt[,,,t] +  arfactor*data$ttt[,,,t-1]
