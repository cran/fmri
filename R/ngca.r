ngca <- function(data,L=c(1000,1000,1000),T=10,m=3,eps=1.5,npca=min(dim(x)[2],dim(x)[1]),filter.time="None",filter.space=FALSE,method="temporal",h.space=3,h.time=3,keepv=FALSE){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(all(class(data)=="fmridata")) {
   x <- extract.data(data)
   mask <- data$mask
   fmriobj <- TRUE
   if(filter.time %in% c("High","Both")){
      d <- dim(x)[4]
      x <- x[,,,-1]-x[,,,-d]
   }
   if(filter.time %in% c("Low","Both")){
      dx <- dim(x)
      cat("Start smoothing in time (Bandwidth=",h.time,")\n")
      x <- .Fortran("smtime",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.time),
                    xnew=double(prod(dx)),
                    double(as.integer(2*h.time+1)),
                    as.integer(2*h.time+1),
                    DUP=FALSE,
                    PACKAGE="fmri")$xnew
     cat("Smoothing in time finished)\n")
     dim(x) <- dx
   }
   if(filter.space) {
      dx <- dim(x)
      cat("Start spatial smoothing (Bandwidth=",h.space,")\n")
      wghts <- data$weights
      x <- .Fortran("smspace",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.space),
                    xnew=double(prod(dx)),
                    as.double(wghts),
                    double(prod(as.integer(2*h.space/wghts+1))),
                    as.integer(as.integer(2*h.space/wghts[1]+1)),
                    as.integer(as.integer(2*h.space/wghts[2]+1)),
                    as.integer(as.integer(2*h.space/wghts[3]+1)),
                    DUP=FALSE,
                    PACKAGE="fmri")$xnew
     cat("Spatial moothing finished)\n")
     dim(x) <- dx
   }
} else if (class(data)%in%c("matrix","array")){
   x <- data
   mask <- TRUE
   fmriobj <- FALSE
} else {
   warning("data has incompatible class argument")
   return(data)
}
#  
#  x - data matrix  (Nxd)
#
cat("Start Non-Gaussian Component Analysis\n")
set.seed(1)
xdim <- dim(x)
lxdim <- length(xdim)
d <- xdim[lxdim]
n <- nn <- prod(xdim[1:(lxdim-1)])
dim(x) <- c(n,d)
mask <- as.vector(mask)
if(length(mask)==1) mask <- rep(mask,n)
x <- x[mask,]
n <- sum(mask)
if(is.null(npca)||npca >= min(d,n)) npca <- min(d,n)
x <- switch(method,"spatial"=sweep(x,1,apply(x,1,mean)),
            "temporal"=sweep(x,2,apply(x,2,mean)))
if(method=="spatial"){
x <- t(x)
d <- dim(x)[2]
n <- dim(x)[1]
}
#svdx <- svd(x,nu=0,nv=npca)
svdx <- svd(x,nu=npca,nv=npca)
#cat("Dimension reduced to:",npca,"\n")
#svdxd <- svdxdinv <- abs(svdx$d[1:npca])
#svdxdinv[svdxdinv>1e-10] <- 1/svdxdinv[svdxdinv>1e-10]
#y <- t(x%*%t(svdx$v)%*%diag(svdxdinv)%*%svdx$v)
if(d>n) {
# reduce to space of first npca components
y <- svdx$v[1:npca,]%*%t(svdx$u)*sqrt(n-1)
} else {
y <- svdx$v%*%t(svdx$u)*sqrt(n-1)
}
#
#
Lsum <- L[1]+L[2]+2*L[3]
s <- c(if(L[1]>0) seq(.5,5,length=L[1]), 
       if(L[2]>0) seq(5/L[2],5,length=L[2]), 
       if(L[3]>0) seq(1/L[3],4,length=L[3]),
       if(L[3]>0) seq(0,4,length=L[3]))
Lsum <- L[1]+L[2]+2*L[3]
ifun <- c(rep(1,L[1]),rep(2,L[2]),rep(3,L[3]),rep(4,L[3]))
#
#   now fast ICA
#
omega <- matrix(rnorm(Lsum*npca),npca,Lsum)
omega <- sweep(omega,2,sqrt(apply(omega^2,2,sum)),"/")
fz <- .Fortran("fastica",
              as.double(y),
              as.double(omega),
              as.integer(npca),
              as.integer(n),
              as.integer(Lsum),
              as.integer(ifun),
              as.integer(T),
              double(npca),
              v=double(npca*Lsum),
              normv=double(Lsum),
              as.double(s),
              DUP=FALSE,
              PACKAGE="fmri")[c("v","normv")]
v <- fz$v
normv <- fz$normv
dim(v) <- c(npca,Lsum)
v <- t(v[,normv>eps])
jhat <- prcomp(v)
if(d>n){
ihat <- svdx$v%*%diag(svdx$d)%*%t(svdx$v[1:npca,])%*%jhat$rotation[,1:m]/sqrt(n-1)
} else {
ihat <- svdx$v%*%diag(svdx$d)%*%t(svdx$v)%*%jhat$rotation[,1:m]/sqrt(n-1)
}
xhat <- x%*%ihat
if(fmriobj){
if(method=="spatial"){
z <- matrix(0,nn,m)
z[mask,]<-ihat
ihat <- array(z,c(xdim[1:3],m))
} else {
z <- matrix(0,nn,m)
z[mask,]<-xhat
xhat <- array(z,c(xdim[1:3],m))
}
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
if(keepv) {
z$v <- v
z$normv <- normv
}
class(z) <- "fmringca"
} else {
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat,jhat=jhat,svdx=svdx)
if(keepv) {
z$v <- v
z$normv <- normv
}
class(z) <- "ngca"
}
z
}

fmriica <- function(data,m=3,method="spatial",xind=NULL,yind=NULL,zind=NULL,tind=NULL,...){
#
#  ICA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(!require(fastICA)) stop("Please install package fastICA from CRAN")
if(all(class(data)=="fmridata")) {
   x <- extract.data(data)
   mask <- data$mask
}  else {
   warning("data has incompatible class argument")
   return(data)
}
if(!is.null(xind)) mask[-xind,,] <- FALSE
if(!is.null(yind)) mask[,-yind,] <- FALSE
if(!is.null(zind)) mask[,,-zind] <- FALSE
#
#  x - data matrix  (Nxd)
#
set.seed(1)
xdim <- dim(x)
lxdim <- length(xdim)
d <- dd <- xdim[lxdim]
n <- nn <- prod(xdim[1:(lxdim-1)])
if(is.null(tind)) tind <- (1:d)
dim(x) <- c(n,d)
mask <- as.vector(mask)
if(length(mask)==1) mask <- rep(mask,n)
x <- x[mask,tind]
n <- sum(mask)
if(method=="spatial"){
x <- t(x)
d <- dim(x)[2]
n <- dim(x)[1]
}
#
#   now fast ICA
#
ttt <- fastICA(x,m,method="C",...)
if(method=="spatial"){
z <- matrix(0,nn,m)
z[mask,]<-ttt$A
ihat <- array(z,c(xdim[1:3],m))
xhat <- matrix(0,m,dd)
xhat[,tind]<-ttt$S
} else {
z <- matrix(0,nn,m)
z[mask,]<-ttt$S
xhat <- array(z,c(xdim[1:3],m))
ihat <- matrix(0,m,dd)
ihat[,tind] <- ttt$A
}
z <- list(ihat=ihat,xhat=xhat)
class(z) <- "fmringca"
z
}

