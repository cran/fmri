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

getkv0 <- function(param,mpredf=mpredfactor,irho=1,alpha=.05,ncmin=2){
   nc <- ncmin:20
   kv <- (param[1]+param[2]*log((alpha+1e-5)/(1-alpha+2e-5)))*mpredf[nc,irho]
   kv
      }

      getalphaclust <- function(alpha,clustertable,ncmin=2){
        ## get reference alpha for clusterthreshold with clustersizes ncmin:20
         calpha <- seq(.001,.1,.001)
         if(ncmin <2) ncmin <- 2
         if(ncmin >20) ncmin <- 20
         cta <- clustertable[,ncmin-1]
         ca0 <- max(calpha[calpha<alpha])
         cta0 <- max(cta[cta<alpha])
         ca1 <- min(calpha[calpha>alpha])
         cta1 <- min(cta[cta>alpha])
         ca0+ (cta1-alpha)*(ca1-ca0)/(cta1-cta0)
      }

      fmri.cluster <- function(spm, alpha=.05, ncmin=2, minimum.signal=0){
        args <- sys.call()
        args <- c(spm$call,args)
      cat("fmri.cluster: entering function\n")

      if (!("fmrispm" %in% class(spm)) ) {
        warning("fmri.cluster: data not of class <fmrispm>. Try to proceed but strange things may happen")
      }

      if (!is.null(attr(spm, "smooth"))) {
        if (!is.null(attr(spm, "residuals"))) {
          type <- "t"
          df <- spm$df
          if(is.null(df)) df <- abs(diff(dim(attr(spm, "design"))))
        } else {
          type <- "norm"
          df <- 1000 # this is actually not needed, placeholder
        }
      } else {
        type <- "t"
        df <- spm$df
      }
        corr <- mean(spm$scorr)
        if(is.null(corr)) corr <- 0
        stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
        dim(stat) <- spm$dim[1:3]
        pv <- 1-switch(type,"norm"=pnorm(stat),"t"=pt(stat,df))
        dim(pv) <- spm$dim[1:3]

      #  clustersizes to use
        clusters <- ncmin:20
        alphaclust <- getalphaclust(alpha,clustertable,ncmin)
        irho <- as.integer(corr/0.05)+1
        if(irho>13) {
           stop("to much spatial correlation")
        }
        # correct for size of multiplicity
        n2 <- sum(spm$mask)
        n1 <- 64^3
        alpha <- 1-(1-alpha)^(n2/n1)
        # this reflects the use of n2 instead of 64^3 voxel (simulation)
        # get critical values
        kv <- getkv0(parcoeff,mpredf=mpredfactor,irho=irho,alpha=alpha,ncmin=2)
      # this gives a vector of kritical values corresponding to cluster sizes
      # now adjust for distribution
        if(type=="t") kv <- qt(pnorm(kv),df)
        detected <- array(0,spm$dim[1:3])
        for(ic in 1:length(clusters)){
           ttt <- findclusters(stat,kv[ic])
           detected[ttt$size>=clusters[ic]] <- 1
           cat("inspecting cluster size",clusters[ic],"detected voxel",sum(detected),"\n")
        }
        detected <- detected*spm$mask
        cat("fmri.pvalue: thresholding\n")
        mask <- rep(FALSE,length=prod(spm$dim[1:3]))
        mask[as.logical(detected)] <- TRUE
        pv[!mask] <- NA
        dim(pv) <- spm$dim[1:3]
          pv[spm$var > 9e19] <- 1

        cat("fmri.pvalue: exiting function\n")

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
        z$dim0 <- spm$dim0
        z$call <- args

        attr(z, "file") <- attr(spm, "file")
        attr(z, "white") <- attr(spm, "white")
        attr(z, "design") <- attr(spm, "design")
          attr(z, "smooth") <- "Not smoothed"
          attr(z, "mode") <- paste("Threshold mode: cluster",ncmin:20,"\n")

        z
      }
