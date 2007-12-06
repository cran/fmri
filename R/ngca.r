ngca <- function(data,L=1000,T=10,m=3,eps=1.5,npca=min(dim(x)[2],dim(x)[1]),method="spatial",sweepmean=NULL,keepv=FALSE){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(all(class(data)=="fmridata")) {
   x <- extract.data(data)
   mask <- data$mask
   fmriobj <- TRUE
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
if(is.null(sweepmean)) sweepmean <- switch(method,"spatial"="temporal","temporal"="spatial")
if(is.null(npca)||npca >= min(d,n)) npca <- min(d,n)
x <- switch(sweepmean,"none"=x,"global"=x - mean(x),"spatial"=sweep(x,2,apply(x,1,mean)),
            "temporal"=sweep(x,2,apply(x,2,mean)))
if(method=="spatial"){
x <- t(x)
d <- dim(x)[2]
n <- dim(x)[1]
}
svdx <- svd(x,nu=0,nv=npca)
#xvar <- var(x)
#z <- svd(xvar)
cat("Dimension reduced to:",npca,"\n")
svdxd <- svdxdinv <- abs(svdx$d[1:npca])
svdxdinv[svdxdinv>0] <- 1/svdxdinv[svdxdinv>0]
y <- t(x%*%svdx$v%*%diag(svdxdinv))
#
#  thats the standardized version of x
#
s <- matrix(0,L,4)
s[,1] <- seq(.5,5,length=L) 
s[,2] <- seq(5/L,5,length=L) 
s[,3] <- seq(4/L,4,length=L) 
s[,4] <- seq(0,4,length=L) 
#
#   now fast ICA
#
omega <- matrix(rnorm(4*L*npca),npca,L*4)
omega <- sweep(omega,2,sqrt(apply(omega^2,2,sum)),"/")
fz <- .Fortran("fastica",
              as.double(y),
              as.double(omega),
              as.integer(npca),
              as.integer(n),
              as.integer(L),
              as.integer(T),
              double(npca),
              v=double(npca*L*4),
              normv=double(L*4),
              as.double(s),
              DUP=FALSE,
              PACKAGE="fmri")[c("v","normv")]
v <- fz$v
normv <- fz$normv
dim(v) <- c(npca,4*L)
v <- t(v[,normv>eps])
jhat <- prcomp(v)
ihat <- svdx$v%*%diag(svdxd)%*%jhat$rotation[,1:m]
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
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
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

