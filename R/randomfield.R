resel <- function(voxeld, hmax, hv=1) {
  reselc <- hv * voxeld / hmax # hv=0.919 for larger bandwidths than typically fmri
  reselc[reselc>1] <- 1
  reselc
}


pvalue <- function(z,i,j,k,rx,ry,rz,type="norm",df=4,cone = 0) {
  rho0 <- 0
  rho1 <- 0
  rho2 <- 0
  rho3 <- 0
  rho4 <- 0
  
  if (type == "norm") {
    rho0 <- 1 - pnorm(z)
    rho1 <- (4 * log(2))^.5 / sqrt(2 * pi) * dnorm(z)
    rho2 <- (4 * log(2)) / (2 * pi) * dnorm(z) * z
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^1.5 * dnorm(z) * (z*z -1)
    rho4 <- (4 * log(2))^2 / (2 * pi)^2 * dnorm(z) * (z*z*z - 3*z)
  }
  if (type == "t") {
    rho0 <- 1 - pt(z,df)
    rho1 <- (4 * log(2))^.5 / (2 * pi) * (1+z*z/df)^(-0.5*(df-1))
    rho2 <- (4 * log(2)) / (2 * pi)^1.5 * gamma((df+1)/2)/gamma(df/2)/(df/2)^0.5 * (1+z*z/df)^(-0.5*(df-1)) * z
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^2 * (1+z*z/df)^(-0.5*(df-1)) * ((df-1)/df*z*z-1)
    rho4 <- (4 * log(2))^2 / (2 * pi)^2.5 * gamma((df+1)/2)/gamma(df/2)/(df/2)^0.5 * (1+z*z/df)^(-0.5*(df-1)) * ((df-2)/df*z*z*z - 3*z)
  }
  if (type == "chisq") {
    rho0 <- 1 - pchisq(z,df)
    rho1 <- (4 * log(2))^.5 / (2 * pi)^.5 * 2 * sqrt(z) * dchisq(z, df)
    rho2 <- (4 * log(2)) / (2 * pi) * 2 * dchisq(z, df) * (z-df+1)
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^1.5 * 2 * dchisq(z, df) / sqrt(z) * (z*z - (2*df-1)*z + (df-1) * (df-2))    
    rho4 <- (4 * log(2))^2 / (2 * pi)^2 * 2 * dchisq(z, df) / z * (z*z*z - 3*df*z*z + 3*(df-1)*z - (df-1)*(df-2)*(df-3))
  }
  
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz
  
  rho0 * r0 + rho1 * r1 + rho2 * r2 + rho3 * r3 +
    (rho1 * r0 + rho2 * r1 + rho3 * r2 + rho4 * r3) / sqrt(4*log(2)) * cone
}


threshold <- function(p,i,j,k,rx,ry,rz,type="norm",df=4,step=.001,cone=0) {
  n <- length(rx)
  thr <- numeric(n)
  fixed <- logical(n)
  if (type == "chisq") {
    x <- 10
  } else {
    x <- 3
  }
  pxyz <- rx*ry*rz
  ind <- (1:length(pxyz))[pxyz==min(pxyz)][1]
  pv <- 1
  while (pv>p) {
     pv <- pvalue(x,i,j,k,rx[ind],ry[ind],rz[ind],type,df,cone=cone)
     x <- x+5*step
  }
  x <- x-5*step
  # this runs faster to the lowest level and therefore allows for smaller value of step
  while (any(!fixed)) {  
    pv <- pvalue(x,i,j,k,rx,ry,rz,type,df,cone=cone)
    ind <- (1:n)[!fixed][pv[!fixed]<p]
    thr[ind] <- x
    fixed[ind] <- TRUE
    x <- x+step
  }
  thr
}
