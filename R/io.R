read.ANALYZE.header <- function(filename) {
  if (is.na(file.info(filename)$size) | (file.info(filename)$size != 348))
    stop("Hmmm! This does not seem to be an ANALYZE header! Wrong size or does not exist!");

  con <- file(filename,"rb")

  header <- list()
  
  # read 4 bytes and get the endianess by comparing filesize with 348
  endian <- if ((sizeofhdr <- readBin(con,"int",1,4,endian="little")) == 348) "little" else "big"
  header$sizeofhdr <- 348
  header$endian <- endian
  
  header$datatype1 <- readChar(con,10)
  header$dbname <- readChar(con,18)
  header$extents <- readBin(con,"int",1,4,endian=endian)
  header$sessionerror <- readBin(con,"int",1,2,endian=endian)
  header$regular <- readChar(con,1)
  header$hkey <- readChar(con,1)
  header$dimension <- readBin(con,"int",8,2,endian=endian)
  header$unused <- readBin(con,"int",7,2,endian=endian)
  header$datatype <- readBin(con,"int",1,2,endian=endian)
  header$bitpix <- readBin(con,"int",1,2,endian=endian)
  header$dimun0 <- readBin(con,"int",1,2,endian=endian)
  header$pixdim <- readBin(con,"double",8,4,endian=endian)
  header$voxoffset <- readBin(con,"double",1,4,endian=endian)
  header$funused <- readBin(con,"double",3,4,endian=endian)
  header$calmax <- readBin(con,"double",1,4,endian=endian)
  header$calmin <- readBin(con,"double",1,4,endian=endian)
  header$compressed <- readBin(con,"double",1,4,endian=endian)
  header$verified <- readBin(con,"double",1,4,endian=endian)
  header$glmax <- readBin(con,"int",1,4,endian=endian)
  header$glmin <- readBin(con,"int",1,4,endian=endian)
  header$describ <- readChar(con,80)
  header$auxfile<- readChar(con,24)
  header$orient <- readChar(con,1) # is this really a character?
  header$originator <- readBin(con,"int",5,2,endian=endian) # documented as 10 byte character!!
  header$generated <- readChar(con,10)
  header$scannum <- readChar(con,10)
  header$patientid <- readChar(con,10)
  header$expdate <- readChar(con,10)
  header$exptime <- readChar(con,10)
  header$histun0 <- readChar(con,3)
  header$views <- readBin(con,"int",1,4,endian=endian)
  header$voladded<- readBin(con,"int",1,4,endian=endian)
  header$startfield <- readBin(con,"int",1,4,endian=endian)
  header$fieldskip <- readBin(con,"int",1,4,endian=endian)
  header$omax <- readBin(con,"int",1,4,endian=endian)
  header$omin <- readBin(con,"int",1,4,endian=endian)
  header$smax <- readBin(con,"int",1,4,endian=endian)
  header$smin <- readBin(con,"int",1,4,endian=endian)

  close(con)

  header
}

read.ANALYZE.volume <- function(filename) {
  file.name <- substring(filename, 1, nchar(filename) - 4)
  file.hdr <- paste(file.name, ".hdr", sep = "")
  file.img <- paste(file.name, ".img", sep = "")

  header <- read.ANALYZE.header(file.hdr)
  dx <- header$dimension[2]
  dy <- header$dimension[3]
  dz <- header$dimension[4]
  dt <- header$dimension[5]
  endian <- header$endian
  if (header$datatype == 1) { # logical
    what <- "raw"
    signed <- TRUE
    size <- 1
  } else if (header$datatype == 4) { # signed short
    what <- "int"
    signed <- TRUE
    size <- if (header$bitpix) header$bitpix/8 else 2
  } else if (header$datatype == 8) { # signed integer
    what <- "int"
    signed <- TRUE
    size <- if (header$bitpix) header$bitpix/8 else 4
  } else if (header$datatype == 16) { # float
    what <- "double"
    signed <- TRUE
    size <- if (header$bitpix) header$bitpix/8 else 4
  } else if (header$datatype == 32) { # complex
    what <- "complex"
    signed <- TRUE
    size <- if (header$bitpix) header$bitpix/8 else 8
  } else if (header$datatype == 64) { # double
    what <- "double"
    signed <- TRUE
    size <- if (header$bitpix) header$bitpix/8 else 8
  } else { # all other
    what <- "raw"
    signed <- TRUE
    size <- 1
  }
  
  con <- file(filename,"rb")
  if (header$datatype == 2) {
    ttt <- readChar(con,file.info(filename)$size)
  } else {
    ttt <- readBin(con, what, n=dx*dy*dz*dt*size, size=size, signed=signed, endian=endian)
  }
  close(con)

  dim(ttt) <- c(dx,dy,dz,dt)
  
  invisible(list(ttt=ttt,header=header))
}

write.ANALYZE.header <- function(header,filename) {
  con <- file(paste(filename, ".hdr", sep=""), "wb")

  writeBin(as.integer(348), con, 4)
  writeChar(header$datatype1, con, 10, NULL)
  writeChar(header$dbname, con, 18, NULL)
  writeBin(as.integer(header$extents), con, 4)
  writeBin(as.integer(header$sessionerror), con, 2)
  writeChar(header$regular, con, 1, NULL)
  writeChar(header$hkey, con, 1, NULL)
  writeBin(as.integer(header$dimension), con, 2)
  writeBin(as.integer(header$unused), con, 2)
  writeBin(as.integer(header$datatype), con, 2)
  writeBin(as.integer(header$bitpix), con, 2)
  writeBin(as.integer(header$dimun0), con, 2)
  writeBin(header$pixdim, con, 4)
  writeBin(header$voxoffset, con, 4)
  writeBin(header$funused, con, 4)
  writeBin(header$calmax, con, 4)
  writeBin(header$calmin, con, 4)
  writeBin(header$compressed, con, 4)
  writeBin(header$verified, con, 4)
  writeBin(as.integer(header$glmax), con, 4)
  writeBin(as.integer(header$glmin), con, 4)
  writeChar(header$describ, con, 80, NULL)
  writeChar(header$auxfile, con, 24, NULL)
  writeChar(header$orient, con, 1, NULL) # is this really a character?
  writeBin(as.integer(header$originator), con, 2) # documented as 10 byte character!!
  writeChar(header$generated, con, 10, NULL)
  writeChar(header$scannum, con, 10, NULL)
  writeChar(header$patientid, con, 10, NULL)
  writeChar(header$expdate, con, 10, NULL)
  writeChar(header$exptime, con, 10, NULL)
  writeChar(header$histun0, con, 3, NULL)
  writeBin(as.integer(header$views), con, 4)
  writeBin(as.integer(header$voladded), con, 4)
  writeBin(as.integer(header$startfield), con, 4)
  writeBin(as.integer(header$fieldskip), con, 4)
  writeBin(as.integer(header$omax), con, 4)
  writeBin(as.integer(header$omin), con, 4)
  writeBin(as.integer(header$smax), con, 4)
  writeBin(as.integer(header$smin), con, 4)
  
  close(con)
}

write.ANALYZE.volume <- function(ttt,filename) {
  con <- file(paste(filename, ".img", sep=""), "wb")
  dim(ttt) <- NULL
  writeBin(as.integer(ttt), con, 2)
  close(con)
}

read.ANALYZE <- function(prefix = "", numbered = FALSE, postfix = "", picstart = 1, numbpic = 1) {
  counter <- c(paste("00", 1:9, sep=""), paste("0", 10:99, sep=""),paste(100:999,sep=""));

  if (numbered) {
    file.img <- paste(prefix, counter[picstart], postfix, ".img", sep="")
  } else {
    file.img <- paste(prefix, ".img", sep="")
  }
  
  if (length(system(paste("ls",file.img),TRUE,TRUE)) != 0) {
    analyze <- read.ANALYZE.volume(file.img);
    ttt <- analyze$ttt
    dt <- dim(ttt)
    cat(".")
    header <- analyze$header;
    
    if ((numbpic > 1) && !numbered) { 
      for (i in (picstart+1):(picstart+numbpic-1)) {
        filename <- paste(prefix, counter[i], postfix, ".img", sep="")
        a <- read.ANALYZE.volume(filename)$ttt
        if (sum() != 0)
          cat("Error: wrong spatial dimension in picture",i)
        ttt <- c(ttt,a);
        dt[4] <- dt[4] + dim(a)[4]
        cat(".")
      }
    }
    
    cat("\n")
    dim(ttt) <- dt
    
    if (min(abs(header$pixdim[2:4])) != 0) {
      weights <-
        abs(header$pixdim[2:4]/min(abs(header$pixdim[2:4])))
    } else {
      weights <- NULL
    }
    
    z <- list(ttt=ttt,format="ANALYZE",delta=header$pixdim[2:4],origin=NULL,orient=NULL,dim=header$dimension[2:5],weights=weights,header=header)
  } else {
    warning(paste("Error: file",filename,"does not exist!"))
    z <- list(ttt=NULL,format="ANALYZE",delta=NULL,origin=NULL,orient=NULL,dim=NULL,weights=NULL,header=NULL)
  }

  mask <- ttt[,,,1]
  mask[mask < quantile(mask,0.75)] <- 0
  mask[mask != 0] <- 1
  z$mask <- mask

  class(z) <- "fmridata"

  if (numbered) {
    attr(z,"file") <- paste(prefix, counter[picstart], postfix, ".hdr/img", "...",
                            prefix, counter[picstart+numbpic-1], postfix, ".hdr/img",
                            sep="")
  } else {
    attr(z,"file") <- paste(prefix, ".hdr/img", sep="")
  }
  
  invisible(z)
}



write.ANALYZE <- function(ttt, header=NULL, file) {
  if (is.null(header)) header <- list()

  if (!("datatype1" %in% names(header))) header$datatype1 <- paste(rep(" ",10),collapse="")
  if (!("dbname" %in% names(header))) header$dbname <- paste(rep(" ",18),collapse="")
  if (!("extents" %in% names(header))) header$extents <- c(0)
  if (!("sessionerror" %in% names(header))) header$sessionerror <- c(0)
  if (!("regular" %in% names(header))) header$regular <- "r"
  if (!("hkey" %in% names(header))) header$hkey <- " "
  if (!("dimension" %in% names(header))) {header$dimension <- rep(0,8); header$dimension <- c(length(dim(ttt)),dim(ttt))}
  if (length(header$dimension) < 8) header$dimension[(length(header$dimension)+1):8] <- 0
  if (!("unused" %in% names(header))) header$unused <- rep(0,7)
  if (length(header$unused) < 7) header$unused[(length(header$unused)+1):7] <- 0
  if (!("datatype" %in% names(header))) header$datatype <- c(4)
  if (!("bitpix" %in% names(header))) header$bitpix <- c(0)
  if (!("dimun0" %in% names(header))) header$dimun0 <- c(0)
  if (!("pixdim" %in% names(header))) header$pixdim <- c(0,4,4,4,rep(0,4))
  if (length(header$pixdim) < 8) header$pixdim[(length(header$pixdim)+1):8] <- 0
  if (!("voxoffset" %in% names(header))) header$voxoffset <- c(0)
  if (!("funused" %in% names(header))) header$funused <- rep(0,3)
  if (length(header$funused) < 3) header$funused[(length(header$funused)+1):3] <- 0
  if (!("calmax" %in% names(header))) header$calmax <- c(0)
  if (!("calmin" %in% names(header))) header$calmin <- c(0)
  if (!("compressed" %in% names(header))) header$compressed <- c(0)
  if (!("verified" %in% names(header))) header$verified <- c(0)
  if (!("glmax" %in% names(header))) header$glmax <- c(0)
  if (!("glmin" %in% names(header))) header$glmin <- c(0)
  if (!("describ" %in% names(header))) header$describ <- paste(rep(" ",80),collapse="")
  if (!("auxfile" %in% names(header))) header$auxfile <- paste(rep(" ",24),collapse="")
  if (!("orient" %in% names(header))) header$orient <- " "
  if (!("originator" %in% names(header))) header$originator <- rep(0,5)
  if (length(header$originator) < 5) header$originator[(length(header$originator)+1):8] <- 0
  if (!("generated" %in% names(header))) header$generated <- paste(rep(" ",10),collapse="")
  if (!("scannum" %in% names(header))) header$scannum <- paste(rep(" ",10),collapse="")
  if (!("patientid" %in% names(header))) header$patientid <- paste(rep(" ",10),collapse="")
  if (!("expdate" %in% names(header))) header$expdate <- paste(rep(" ",10),collapse="")
  if (!("exptime" %in% names(header))) header$exptime <- paste(rep(" ",10),collapse="")
  if (!("histun0" %in% names(header))) header$histun0 <- paste(rep(" ",3),collapse="")
  if (!("views" %in% names(header))) header$views <- c(0)
  if (!("voladded" %in% names(header))) header$voladded <- c(0)
  if (!("startfield" %in% names(header))) header$startfield <- c(0)
  if (!("fieldskip" %in% names(header))) header$fieldskip <- c(0)
  if (!("omax" %in% names(header))) header$omax <- c(0)
  if (!("omin" %in% names(header))) header$omin <- c(0)
  if (!("smax" %in% names(header))) header$smax <- c(0)
  if (!("smin" %in% names(header))) header$smin <- c(0)

  write.ANALYZE.header(header,file)

  write.ANALYZE.volume(ttt, file)
}



read.AFNI <- function(file) {
  conhead <- file(paste(file,".HEAD",sep=""),"r")
  header <- readLines(conhead)
  close(conhead)

  types <- NULL
  args <- NULL
  counts <- NULL
  values <- NULL
  
  for (i in 1:length(header)) {
    if (regexpr("^type *= *", header[i]) != -1) {
      tmptype <- strsplit(header[i]," *= *")[[1]][2]
      types <- c(types,tmptype)
      args <- c(args,strsplit(header[i+1]," *= *")[[1]][2])
      tmpcounts <- as.numeric(strsplit(header[i+2]," *= *")[[1]][2])
      counts <- c(counts,tmpcounts)
      i <- i+3
      tmpvalue <- ""
      while ((regexpr("^$", header[i]) == -1) && (i <= length(header))) {
        tmpvalue <- paste(tmpvalue,header[i])
        i <- i+1
      }
      tmpvalue <- sub("^ +","",tmpvalue)
      if ((tmptype == "integer-attribute") || (tmptype == "float-attribute")) {
        tmpvalue <- as.numeric(strsplit(tmpvalue," +")[[1]])
      }
      values <- c(values,list(value=tmpvalue))
    }        
  }

  names(values) <- args

  dx <- values$DATASET_DIMENSIONS[1]
  dy <- values$DATASET_DIMENSIONS[2]
  dz <- values$DATASET_DIMENSIONS[3]
  dt <- values$DATASET_RANK[2]
  scale <- values$BRICK_FLOAT_FACS
  size <- file.info(paste(file,".BRIK",sep=""))$size/(dx*dy*dz*dt)

  if (regexpr("MSB",values$BYTEORDER_STRING[1]) != -1) {
    endian <- "big"
  } else {
    endian <- "little"
  }

  if (min(abs(values$DELTA)) != 0) {
    weights <-
      abs(values$DELTA/min(abs(values$DELTA)))
  } else {
    weights <- NULL
  }
  
  if (as.integer(size) == size) {
    conbrik <- file(paste(file,".BRIK",sep=""),"rb")
    myttt<- readBin(conbrik, "int", n=dx*dy*dz*dt*size, size=size, signed=TRUE, endian=endian)
    close(conbrik)
    dim(myttt) <- c(dx,dy,dz,dt)
    for (k in 1:dt) {
      if (scale[k] != 0) {
        cat("scale",k,"with",scale[k],"\n")
        cat(range(myttt[,,,k]),"\n")
        myttt[,,,k] <- scale[k] * myttt[,,,k]
        cat(range(myttt[,,,k]),"\n")
      }
    }
    z <-
      list(ttt=myttt,format="HEAD/BRIK",delta=values$DELTA,origin=values$ORIGIN,orient=values$ORIENT_SPECIFIC,dim=c(dx,dy,dz,dt),weights=weights, header=values)
  } else {
    warning("Error reading file: Could not detect size per voxel\n")
    z <- list(ttt=NULL,format="HEAD/BRIK",delta=NULL,origin=NULL,orient=NULL,dim=NULL,weights=NULL,header=values)    
  }

  mask <- myttt[,,,1]
  mask[mask < quantile(mask,0.75)] <- 0
  mask[mask != 0] <- 1
  z$mask <- mask

  class(z) <- "fmridata"
  attr(z,"file") <- paste(file,".HEAD/BRIK",sep="")
  invisible(z)
}




write.AFNI <- function(file, ttt, label, note="", origin=c(0,0,0), delta=c(4,4,4), idcode="WIAS_noid") {
  ## TODO:
  ## 
  ## create object oriented way!!!!
  
  AFNIheaderpart <- function(type, name, value) {
    a <- "\n"
    a <- paste(a, "type = ", type, "\n", sep="")
    a <- paste(a, "name = ", name, "\n", sep="")
    if (regexpr("string",type) == 1) {
      value <- paste("'", value, "~", sep="")
      a <- paste(a, "count = ", nchar(value) - 1, "\n", sep ="")
      a <- paste(a, value, "\n", sep="")
    } else {
      a <- paste(a, "count = ", length(value), "\n", sep ="")
      j <- 0
      while (j<length(value)) {
        left <- length(value) - j
        if (left>4) left <- 5
        a <- paste(a, paste(value[(j+1):(j+left)],collapse="  "), "\n", sep="  ")
        j <- j+5
      }
    }
    a
  }
  
  conhead <- file(paste(file, ".HEAD", sep=""), "w")
  writeChar(AFNIheaderpart("string-attribute","HISTORY_NOTE",note),conhead,eos=NULL)
  writeChar(AFNIheaderpart("string-attribute","TYPESTRING","3DIM_HEAD_FUNC"),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_STRING",idcode),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_DATE",date()),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","SCENE_DATA",c(0,11,1,-999,-999,-999,-999,-999)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","ORIENT_SPECIFIC",c(0,3,4)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","ORIGIN",origin),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","DELTA",delta),conhead,eos=NULL)  
  minmax <- function(y) {r <- NULL;for (k in 1:dim(y)[4]) {r <- c(r,min(y[,,,k]),max(y[,,,k]))}; r}
  mm <- minmax(ttt)
  writeChar(AFNIheaderpart("float-attribute","BRICK_STATS",mm),conhead,eos=NULL)
  writeChar(AFNIheaderpart("integer-attribute","DATASET_RANK",c(3,dim(ttt)[4],0,0,0,0,0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","DATASET_DIMENSIONS",c(dim(ttt)[1:3],0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","BRICK_TYPES",rep(1,dim(ttt)[4])),conhead,eos=NULL)  

  scale <- rep(0,dim(ttt)[4])
  for (k in 1:dim(ttt)[4]) {
    scale[k] <- max(abs(mm[2*k-1]),abs(mm[2*k]))/32767
    ttt[,,,k] <- ttt[,,,k] / scale[k]
  }

  writeChar(AFNIheaderpart("float-attribute","BRICK_FLOAT_FACS",scale),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BRICK_LABS",paste(label,collapse="~")),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BYTEORDER_STRING","MSB_FIRST"),conhead,eos=NULL)  
  close(conhead)

  conbrik <- file(paste(file, ".BRIK", sep=""), "wb")
  dim(ttt) <- NULL
  writeBin(as.integer(ttt), conbrik,size=2, endian="big")
  close(conbrik)
}


