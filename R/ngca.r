ngca <- function(x,L=1000,T=10,m=3,eps=1.5,npca=dim(x)[2],keepv=FALSE){
#
#  NGCA algorithm 
#
#  
#  x - data matrix  (Nxd)
#
set.seed(1)
xdim <- dim(x)
d <- xdim[2]
n <- xdim[1]
xmean <- apply(x,2,mean)
xvar <- var(x)
y <- sweep(x,2,xmean)
z <- svd(xvar)
cat("Dimension reduced to:",npca,"\n")
y <- t(y%*%z$u%*%diag(z$d^(-.5))[,1:npca])
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
ihat <- z$u[,1:npca]%*%jhat$rotation[,1:m]
xhat <- x%*%ihat
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
if(keepv) {
z$v <- v
z$normv <- normv
}
class(z) <- "ngca"
z
}
