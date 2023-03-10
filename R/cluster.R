findclusters <- function(x,thresh){
   dx <- dim(x)
   tx <- as.integer(x>thresh)
   tx[is.na(tx)] <- 0
   z <- .Fortran(C_ccluster,
                 size=as.integer(tx),
                 as.integer(dx[1]),
                 as.integer(dx[2]),
                 as.integer(dx[3]),
                 clusterid=integer(prod(dx)))[c("size","clusterid")]
    dim(z$size) <- dx
    dim(z$clusterid) <- dx
    z
}

getpvalue <- function(x,xclust,n,cc,nc,clusters){
# x - spm array
# xclust result from findclusters
# nc specified actual cluster size
# nca/nce specified min/max cluster size
# cc spatial correlation
   pvalue <- array(1,dim(x))
   clusterids <- unique(xclust$clusterid[xclust$size>=nc])
   for(i in clusterids){
      if(xclust$size[xclust$clusterid==i][1]==nc){
         xvalue <- min(x[xclust$clusterid==i])
      } else {
         vals <- x[xclust$clusterid==i]
         j <- length(vals)-nc+1
         xvalue <- sort(vals,partial=j)[j]
# assuming, that the largest nc values still form a connected cluster
      }
      pvalue[xclust$clusterid==i] <- pvclust(xvalue,n,cc,clusters)[clusters==nc]
   }
   pvalue
}

kvclust <- function(alpha,n,cc,nca,nce=max(nca)){
    # compute critical value for level alpha and 
    # sample size n, cluster sizes nca:nce and spatial correlation cc
    # based on tail approximation Abramovich/Stegun 26.2.14
    nca <- min(nca):nce
    th11 <- -0.5679493  -0.6276305*cc^.3  + 3.0608387*n^(1/9)
    th12 <- -8.47163638 + 0.02560351*n^(1/3) + 84.28136758*cc^.9
    th13 <- 3.18491717 -19.41570878*cc^.9 + 0.00479463*n^(1/3)  
    th14 <- -0.4531994 + 4.7408157*cc^(2/3) -0.1656168*n^(1/15)
    lnc <- log(nca)
    lnc2 <- lnc^2
    lnc3 <- lnc^3
    th1 <- th11 + th12*(-.7961+.3573*lnc) + th13*(1.8-2.108*lnc+0.539*lnc2) + 
      th14*(-4.316+ 8.618*lnc-4.993*lnc2+0.880*lnc3 )
    th2 <- 0.7071393  + 0.5875162*nca^.85  -5.1295360*cc^.75
    nca <- nca[1]
    if(nce>nca){
      ncd <- log(1e-4+nce-nca)
      nca2 <- nca^2
      nca3 <- nca^3
      eta1 <- pmax(0,-0.014020 + 0.008642*n^(1/3) + 0.157975*ncd + 
                     0.114501*nca -0.015963*nca2 +0.000520*nca3)
      eta2 <- pmin(0, -0.2496057 -0.0485483*ncd + 0.0663842*nca -0.0056151*nca2 +0.0001511*nca3)
       th1 <- th1+eta1+eta2*log(alpha)
    }
    x <- 2
    x0 <- 0
    i <- 1
    while(max(abs(x-x0))>1e-5&i<10){
      x0 <- x
      x <- sqrt(pmax(1,(-log(alpha)+th1-log(x))/th2))
#      cat("i",i,"kv",signif(x,3),"\n")
      i <- i+1
    }
    x
  }

pvclust <- function(tvalue,n,cc,nca,nce=max(nca)){
  # compute p-value for test statistic tvalue and 
  # sample size n, cluster size nc and spatial correlation cc
  # based on tail approximation Abramovich/Stegun 26.2.1
  # nce: maximum cluster size 
  # nc: actual cluster size under consideration
  nca <- min(nca):nce
  th11 <- -0.5679493  -0.6276305*cc^.3  + 3.0608387*n^(1/9)
  th12 <- -8.47163638 + 0.02560351*n^(1/3) + 84.28136758*cc^.9
  th13 <- 3.18491717 -19.41570878*cc^.9 + 0.00479463*n^(1/3)  
  th14 <- -0.4531994 + 4.7408157*cc^(2/3) -0.1656168*n^(1/15)
  lnc <- log(nca)
  lnc2 <- lnc^2
  lnc3 <- lnc^3
  th1 <- th11 + th12*(-.7961+.3573*lnc) + th13*(1.8-2.108*lnc+0.539*lnc2) + 
    th14*(-4.316+ 8.618*lnc-4.993*lnc2+0.880*lnc3 )
  th2 <- 0.7071393  + 0.5875162*nca^.85  -5.1295360*cc^.75
  nca <- nca[1]
  if(nce>nca){
    ncd <- log(1e-4+nce-nca)
    nca2 <- nca^2
    nca3 <- nca^3
    eta1 <- pmax(0,-0.014020 + 0.008642*n^(1/3) + 0.157975*ncd + 
                   0.114501*nca -0.015963*nca2 +0.000520*nca3)
    eta2 <- pmin(0, -0.2496057 -0.0485483*ncd + 0.0663842*nca -0.0056151*nca2 +0.0001511*nca3)
  } else {
    eta1 <- eta2 <- 0
  }
  pmin(.5,exp((th1+eta1-th2*tvalue^2-log(tvalue))/(1-eta2)))
}




fmri.cluster <- function(spm, alpha=.05, ncmin=2, ncmax=ncmin, minimum.signal=0, verbose = FALSE){
        args <- sys.call()
        args <- c(spm$call,args)
      if (verbose) cat("fmri.cluster: entering function\n")

      if (!inherits(spm,"fmrispm"))  {
        warning("fmri.cluster: data not of class <fmrispm>. Try to proceed but strange things may happen")
      }

      if (!is.null(attr(spm, "smooth"))) return("fmri.cluster not suitable for smoothed spm's") 
        corr <- mean(spm$scorr)
        if(is.null(corr)) corr <- 0
        stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
        dim(stat) <- spm$dim[1:3]
        pv <- array(1,dim(stat))
      #  clustersizes to use
        clusters <- ncmin:ncmax
        irho <- as.integer(corr/0.05)+1
        if(irho>6) {
           stop("to much spatial correlation")
        }
        # correct for size of multiplicity
        n <- sum(spm$mask)
        kv <- kvclust(alpha,n,corr,clusters)
      # this gives a vector of kritical values corresponding to cluster sizes
        detected <- array(0,spm$dim[1:3])
        for(ic in 1:length(clusters)){
           ttt <- findclusters(stat,kv[ic])
           detected[ttt$size>=clusters[ic]] <- 1
           pv <- pmin(pv,getpvalue(stat,ttt,n,corr,clusters[ic],clusters))
           if (verbose) cat("inspecting cluster size",clusters[ic],"detected voxel",sum(detected),"\n")
        }
        detected <- detected*spm$mask
        if (verbose) cat("fmri.pvalue: thresholding\n")
        mask <- rep(FALSE,length=prod(spm$dim[1:3]))
        mask[as.logical(detected)] <- TRUE
        pv[!mask] <- NA
        dim(pv) <- spm$dim[1:3]

        if (verbose) cat("fmri.pvalue: exiting function\n")

        z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim,
                  hrf = spm$hrf, alpha=alpha, mask = spm$mask)

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
          attr(z, "mode") <- paste("Threshold mode: cluster",ncmin:20,"\n")

        z
      }

simclusterthr <- function(n,distr=c("norm","t"),df=100,bw=0,kern="Gaussian",
                          nsim=1000,seed=1,qgrid = pnorm(seq(1.4,4,.05)),filename="tempres"){
        nq <- length(qgrid)
        nsize <- 20
        results <- array(0,c(nq,nsize))
        set.seed(seed)
        dx <- c(n,n,n)
        quant <- if(distr[1]=="t") qt(qgrid,df) else qnorm(qgrid)
        for(i in 1:nsim){
          x <- if(distr[1]=="t") rt(n^3,df) else rnorm(n^3)
          dim(x) <- dx
          if(bw>0) {
            x <- kernsm(x,bw,kern=kern,unit="FWHM")@yhat
            x <- x/sqrt(var(x))
          }
          for(j in 1:nq){
            mc <- min(20,max(findclusters(x,quant[j])$size))
            results[j,1:mc] <- results[j,1:mc]+1
          }
          cat("Time",format(Sys.time()),"after ",i,"simulations:\n")
          if(i%/%100*100==i) print(results[seq(1,nq,5),]/i)
          if(i%/%1000*1000==i){
            tmp <- list(results=results/i,qgrid=qgrid,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
            save(tmp,file=filename)
          }
        }
        list(results=results/nsim,qgrid=qgrid,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
      }
      
simclusterthr2 <- function(n,kv,distr=c("norm","t"),df=100,bw=0,kern="Gaussian",
                           nsim=1000,seed=1,filename="tempres"){
  #
  #  second order simulations for adjustment of kritical values
  #
  nq <- length(kv)
  nsize <- 20
  results <- numeric(nq)
  set.seed(seed)
  dx <- c(n,n,n)
  for(i in 1:nsim){
    x <- if(distr[1]=="t") rt(n^3,df) else rnorm(n^3)
    dim(x) <- dx
    if(bw>0) {
      x <- kernsm(x,bw,kern=kern,unit="FWHM")@yhat
      x <- x/sqrt(var(x))
    }
    mc <- logical(nq)
    for(j in 1:nq){
      mc[j] <- max(findclusters(x,kv[j])$size)>j
    }
    for(j in 1:nq){
      results[j] <- results[j] + any(mc[j:nq])
    }
    cat("Time",format(Sys.time()),"after ",i,"simulations:\n")
    if(i%/%100*100==i) print(results/i)
    if(i%/%1000*1000==i){
      tmp <- list(results=results/i,kv=kv,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
      save(tmp,file=filename)
    }
  }
  list(results=results/nsim,kv=kv,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
}
