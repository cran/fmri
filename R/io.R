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
  (dx <- header$dimension[2]) || (dx <- 1)
  (dy <- header$dimension[3]) || (dy <- 1)
  (dz <- header$dimension[4]) || (dz <- 1)
  (dt <- header$dimension[5]) || (dt <- 1)
  endian <- header$endian
  if (header$datatype == 1) { # logical
    what <- "raw"
    signed <- TRUE
    size <- 1
  } else if (header$datatype == 2) { # unsigned char????
    what <- "int"
    signed <- FALSE
    size <- if (header$bitpix) header$bitpix/8 else 2
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

write.ANALYZE.volume <- function(ttt,filename,size) {
  con <- file(paste(filename, ".img", sep=""), "wb")
  dim(ttt) <- NULL
  writeBin(as.integer(ttt), con, size)
  close(con)
}

read.ANALYZE <- function(prefix = c(""), numbered = FALSE, postfix = "", picstart = 1, numbpic = 1) {
  counter <- c(paste("00", 1:9, sep=""), paste("0", 10:99, sep=""),paste(100:999,sep=""));

  prefix <- strsplit(prefix,".img")
  filename <- character(length(prefix))
  if (numbered) {
    for (i in picstart:(picstart+numbpic-1)) {
      filename[i] <- paste(prefix[[1]][1], counter[i], postfix, ".img", sep="")
    }
  } else {
    for (i in 1:length(prefix)) {
      if (length(prefix[[i]]) > 1)
        warning("filename",paste(prefix[[i]],collapse=""),"probably misinterpreted due to extension-like .img in name!")
      filename[i] <- paste(prefix[[i]][1], ".img", sep="")
    }
  }

  if (!is.na(file.info(filename[1])$size)) {
    analyze <- read.ANALYZE.volume(filename[1]);
    ttt <- analyze$ttt
    dt <- dim(ttt)
    cat(".")
    header <- analyze$header;
    
    if ((numbpic > 1) && !numbered) { 
      for (i in 2:numbpic) {
        a <- read.ANALYZE.volume(filename[i])$ttt
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
    
     mask <- array(TRUE,dt[1:3])
     mask[ttt[,,,1] < quantile(ttt[,,,1],0.75)] <- FALSE

    z <- list(ttt=writeBin(as.numeric(ttt),raw(),4),format="ANALYZE",delta=header$pixdim[2:4],origin=NULL,
              orient=NULL,dim=dim(ttt),weights=weights,header=header, mask=mask)
  } else {
    warning(paste("Error: file",filename[1],"does not exist!"))
    z <- list(ttt=NULL,format="ANALYZE",delta=NULL,origin=NULL,orient=NULL,dim=NULL,weights=NULL,header=NULL,mask=NULL)
  }

  class(z) <- "fmridata"

  if (length(filename) > 1) {
    attr(z,"file") <- paste(filename[1], "...",
                            filename[length(filename)],
                            sep="")
  } else {
    attr(z,"file") <- paste(filename[1], sep="")
  }
  
  invisible(z)
}



write.ANALYZE <- function(ttt, header=NULL, filename) {
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

  if (header$datatype == 1) { # logical
    what <- "raw"
    signed <- TRUE
    size <- 1
  } else if (header$datatype == 2) { # unsigned char????
    what <- "int"
    signed <- FALSE
    size <- if (header$bitpix) header$bitpix/8 else 2
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

  write.ANALYZE.header(header,filename)

  write.ANALYZE.volume(ttt, filename,size)
}



read.AFNI <- function(filename) {
  fileparts <- strsplit(filename,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])

  if (ext == "head") {
    filename.head <- filename
    filename.brik <- paste(c(fileparts[-length(fileparts)],"BRIK"),collapse=".")
  } else if (ext == "brik") {
    filename.head <- paste(c(fileparts[-length(fileparts)],"HEAD"),collapse=".")
    filename.brik <- filename
  } else {
    filename.head <- paste(filename,".HEAD",sep="")
    filename.brik <- paste(filename,".BRIK",sep="")
  }
  
  conhead <- file(filename.head,"r")
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
      } else {
        tmpvalue <- sub("~$","",sub("^\'","",tmpvalue))
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
  size <- file.info(filename.brik)$size/(dx*dy*dz*dt)

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
    conbrik <- file(filename.brik,"rb")
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
  mask <- array(TRUE,c(dx,dy,dz))
  mask[myttt[,,,1] < quantile(myttt[,,,1],0.75)] <- FALSE
    z <-
      list(ttt=writeBin(as.numeric(myttt),raw(),4),format="HEAD/BRIK",delta=values$DELTA,origin=values$ORIGIN,orient=values$ORIENT_SPECIFIC,dim=c(dx,dy,dz,dt),weights=weights, header=values,mask=mask)
  } else {
    warning("Error reading file: Could not detect size per voxel\n")
    z <- list(ttt=NULL,format="HEAD/BRIK",delta=NULL,origin=NULL,orient=NULL,dim=NULL,weights=NULL,header=values,mask=NULL)    
  }

  class(z) <- "fmridata"
  attr(z,"file") <- paste(filename,".HEAD/BRIK",sep="")
  invisible(z)
}

write.AFNI <- function(filename, ttt, label=NULL, note=NULL, origin=NULL, delta=NULL, idcode=NULL, header=NULL, taxis = FALSE) {

  afni.header <- list()
  # afni.header$NAME <- c("NAME", "type", category, count, reserved-count, multiple)
  # "NAME"         : Attribute name (string)
  # "type"         : Attribute type (one of "string", "integer", "float")
  # category       : according to AFNI Doc. 0 is mandatory, 1 only for time series
  # count          : number of parameters (used), 0 is unknown, or not implemented
  # reserved-count : number of parameters (reserved), 0 is unknown, or not implemented
  # multiple       : is this a multiple argument like NOTE_NUMBER_001

  # Mandatory Attributes
  afni.header$DATASET_RANK <- c("DATASET_RANK", "integer", 0, 2, 8, FALSE)
  afni.header$DATASET_DIMENSIONS <- c("DATASET_DIMENSION", "integer", 0, 3, 5, FALSE)
  afni.header$TYPESTRING <- c("TYPESTRING", "string", 0, 1, 1, FALSE)
  afni.header$SCENE_DATA <- c("SCENE_DATA", "integer", 0, 3, 8, FALSE)
  afni.header$ORIENT_SPECIFIC <- c("ORIENT_SPECIFIC", "integer", 0, 3, 3, FALSE)
  afni.header$ORIGIN <- c("ORIGIN", "float", 0, 3, 3, FALSE)
  afni.header$DELTA <- c("DELTA", "float", 0, 3, 3, FALSE)

  # Time-dependent Dataset Attributes
  afni.header$TAXIS_NUMS <- c("TAXIS_NUMS", "integer", 1, 3, 8, FALSE)
  afni.header$TAXIS_FLOATS <- c("TAXIS_FLOATS", "float", 1, 5, 8, FALSE)
  afni.header$TAXIS_OFFSETS <- c("TAXIS_OFFSETS", "float", 1, 0, 0, FALSE)

  # Almost Mandatory Attributes
  afni.header$IDCODE_STRING <- c("IDCODE_STRING", "string", 2, 1, 1, FALSE)
  afni.header$IDCODE_DATE <- c("IDCODE_DATE", "string", 2, 1, 1, FALSE)
  afni.header$BYTEORDER_STRING <- c("BYTEORDER_STRING", "string", 2, 1, 1, FALSE)
  afni.header$BRICK_STATS <- c("BRICK_STATS", "float", 2, 0, 0, FALSE)
  afni.header$BRICK_TYPES <- c("BRICK_TYPES", "integer", 2, 0, 0, FALSE)
  afni.header$BRICK_FLOAT_FACS <- c("BRICK_FLOAT_FACS", "float", 2, 0, 0, FALSE)
  afni.header$BRICK_LABS <- c("BRICK_LABS", "string", 2, 0, 0, FALSE)
  afni.header$BRICK_STATAUX <- c("BRICK_STATAUX", "float", 2, 0, 0, FALSE)
  afni.header$STAT_AUX <- c("STAT_AUX", "float", 2, 0, 0, FALSE) # depreciated

  # Note Attributes
  afni.header$HISTORY_NOTE <- c("HISTORY_NOTE", "string", 3, 1, 1, FALSE)
  afni.header$NOTES_COUNT <- c("NOTES_COUNT", "integer", 3, 1, 1, FALSE)
  afni.header$NOTE_NUMBER <- c("NOTE_NUMBER", "string", 3, 0, 0, TRUE) # more than one!!!

  # Registration Attributes
  afni.header$TAGALIGN_MATVEC <- c("TAGALIGN_MATVEC", "float", 4, 12, 12, FALSE)
  afni.header$VOLREG_MATVEC <- c("VOLREG_MATVEC", "float", 4, 12, 12, TRUE) # more than one!!!
  afni.header$VOLREG_ROTCOM <- c("VOLREG_ROTCOM", "string", 4, 1, 1, TRUE) # more than one!!!
  afni.header$VOLREG_CENTER_OLD <- c("VOLREG_CENTER_OLD", "float", 4, 3, 3, FALSE)
  afni.header$VOLREG_CENTER_BASE <- c("VOLREG_CENTER_BASE", "float", 4, 3, 3, FALSE)
  afni.header$VOLREG_ROTPARENT_IDCODE <- c("VOLREG_ROTPARENT_IDCODE", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_ROTPARENT_NAME <- c("VOLREG_ROTPARENT_NAME", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_GRIDPARENT_IDCODE <- c("VOLREG_GRIDPARENT_IDCODE", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_GRIDPARENT_NAME <- c("VOLREG_GRIDPARENT_NAME", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_INPUT_IDCODE <- c("VOLREG_INPUT_IDCODE", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_INPUT_NAME <- c("VOLREG_INPUT_NAME", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_BASE_IDCODE <- c("VOLREG_BASE_IDCODE", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_BASE_NAME <- c("VOLREG_BASE_NAME", "string", 4, 1, 1, FALSE)
  afni.header$VOLREG_ROTCOM_NUM <- c("VOLREG_ROTCOM_NUM", "integer", 4, 1, 1, FALSE)

  # Miscellaneous Attributes
  afni.header$IDCOE_ANAT_PARENT <- c("IDCOE_ANAT_PARENT", "string", 5, 1, 1, FALSE)
  afni.header$TO3D_ZPAD <- c("TO3D_ZPAD", "integer", 5, 3, 3, FALSE)

  # Warping Attributes
  afni.header$IDCOE_WARP_PARENT <- c("IDCOE_WARP_PARENT", "string", 6, 1, 1, FALSE)
  afni.header$WARP_TYPE <- c("WARP_TYPE", "integer", 6, 2, 2, FALSE)
  afni.header$WARP_DATA <- c("WARP_DATA", "float", 6, 0, 0, FALSE)

  # Talairach Markers Attributes
  afni.header$MARKS_XYZ <- c("MARKS_XYZ", "float", 7, 30, 30, FALSE)
  afni.header$MARKS_LAB <- c("MARKS_LAB", "string", 7, 1, 1, FALSE)
  afni.header$MARKS_HELP <- c("MARKS_HELP", "string", 7, 1, 1, FALSE)
  afni.header$MARKS_FLAGS <- c("MARKS_FLAGS", "integer", 7, 2, 2, FALSE)

  # Attributes for User-Defined Tags
  afni.header$TAGSET_NUM <- c("TAGSET_NUM", "integer", 8, 2, 2, FALSE)
  afni.header$TAGSET_FLOATS <- c("TAGSET_FLOATS", "floats", 8, 0, 0, FALSE)
  afni.header$TAGSET_LABELS <- c("TAGSET_LABELS", "string", 8, 1, 1, FALSE)

  # Nearly Useless Attributes
  afni.header$LABEL_1 <- c("LABEL_1", "string", 9, 1, 1, FALSE)
  afni.header$LABEL_2 <- c("LABEL_2", "string", 9, 1, 1, FALSE)
  afni.header$DATASET_NAME <- c("DATASET_NAME", "string", 9, 1, 1, FALSE)
  afni.header$DATASET_KEYWORDS <- c("DATASET_KEYWORDS", "string", 9, 1, 1, FALSE)
  afni.header$BRICK_KEYWORDS <- c("BRICK_KEYWORDS", "string", 9, 1, 1, FALSE)
  
  AFNIheaderpart <- function(name, value, conhead, check=NULL) {
    if (regexpr("_[0-9]*$", name)) if (!is.null(afni.header[[sub("_[0-9]*$","", name)]])) {
      if (as.logical(afni.header[[sub("_[0-9]*$","", name)]][6])) {
        header.entry <- afni.header[[sub("_[0-9]*$","", name)]]
      } else {
        header.entry <- afni.header[[name]]
      }
    } else {
      header.entry <- afni.header[[name]]
    }
    if (!(is.null(header.entry[2]))) {
      # now we know, that entry "name" is valid attribute
      a <- "\n"
      type <- header.entry[2]
      a <- paste(a, "type = ", type, "-attribute\n", sep="")
      a <- paste(a, "name = ", name, "\n", sep="")
      if (is.null(value)) {
        if (header.entry[3] == 0) {
          stop("not found mandatory attribute ",name," in header", call.=FALSE)
        } else if ((header.entry[3] == 1) && taxis) {
          stop("not found mandatory attribute ",name," in header", call.=FALSE)
        } else {
          stop("attempt to write attribute ",name,", which was not found", call.=FALSE)
        }
      }

      if (regexpr("string",type) == 1) {
        if (!is.null(check)) if (!(value %in% check)) 
          warning("value ", value, " for attribute ", name," does not seem to be valid! Please check!", call.=FALSE)
        value <- paste("'", value, "~", sep="")                         # make syntax highlightening work:'
        a <- paste(a, "count = ", nchar(value) - 1, "\n", sep ="")
        a <- paste(a, value, "\n", sep="")
      } else {
        if (header.entry[4] != 0) {
          if (length(value) < header.entry[4]) {
            warning("too few values for ",name,", expecting ", header.entry[4] ," values, filling with zeros!\n", call.=FALSE)
          }
          length(value) <- as.integer(header.entry[5])
          value[is.na(value)] <- 0
          if (!is.null(check)) {
            length(check) <- as.integer(header.entry[5])
            check[is.na(check)] <- 0
            if (!as.logical(prod((value[1:as.integer(header.entry[4])] == check[1:as.integer(header.entry[4])])))) 
              warning("value ", value, " for attribute ", name," does not seem to be valid or match dataset! Please check!", call.=FALSE)
          }
        }
        a <- paste(a, "count = ", length(value), "\n", sep ="")
        j <- 0
        while (j<length(value)) {
          left <- length(value) - j
          if (left>4) left <- 5
          a <- paste(a, paste(value[(j+1):(j+left)],collapse="  "), "\n", sep="  ")
          j <- j+5
        }
      }
      writeChar(a,conhead,eos=NULL)
    } else {
      warning("attempt to write unknown attribute ",name," ...skipping!", call.=FALSE)
    }
  }

  if (!is.null(c(label, note, origin, delta, idcode))) warning("The use of any of the arguments label, note, origin, delta, idcode is depreciated. Please use header and taxis argument, see documentation. THEY WILL VANISH IN SOME FUTURE RELAESE OF THIS SOFTWARE!", call.=FALSE)

  # checking dataset!!!
  if ((length(dim(ttt)) < 3) || (length(dim(ttt)) > 4) || any(dim(ttt)[1:3] == 1)) stop("bad dimension",dim(ttt) ,"for dataset!")
  if (length(dim(ttt)) == 3) dim(ttt) <- c(dim(ttt),1)

  if (is.null(header)) header <- list()

  conhead <- file(paste(filename, ".HEAD", sep=""), "w")
  # write mandatory attributes
  if (is.null(header$DATASET_RANK)) header$DATASET_RANK <- c(3,dim(ttt)[4])
  AFNIheaderpart("DATASET_RANK",header$DATASET_RANK, conhead,c(3,dim(ttt)[4]))
  header$DATASET_RANK <- NULL
  if (is.null(header$DATASET_DIMENSIONS)) header$DATASET_DIMENSIONS <- c(dim(ttt)[1:3])
  AFNIheaderpart("DATASET_DIMENSIONS",header$DATASET_DIMENSIONS, conhead,c(dim(ttt)[1:3]))
  header$DATASET_DIMENSIONS <- NULL
  if (is.null(header$TYPESTRING)) {header$TYPESTRING <- "3DIM_HEAD_FUNC"; warning("TYPESTRING not given, setting to 3DIM_HEAD_FUNC for backward compatibility", call.=FALSE) }
  AFNIheaderpart("TYPESTRING",header$TYPESTRING, conhead,c("3DIM_HEAD_ANAT","3DIM_HEAD_FUNC","3DIM_GEN_ANAT","3DIM_GEN_ANAT"))
  header$TYPESTRING <- NULL
  if (is.null(header$SCENE_DATA)) {header$SCENE_DATA <- c(0,11,1,-999,-999,-999,-999,-999); warning("SCENE_DATA not given, setting to c(0,11,1,-999,-999,-999,-999,-999) for backward compatibility", call.=FALSE) }
  AFNIheaderpart("SCENE_DATA",header$SCENE_DATA, conhead)
  header$SCENE_DATA <- NULL
  if (is.null(header$ORIENT_SPECIFIC)) {header$ORIENT_SPECIFIC <- c(0,3,4); warning("ORIENT_SPECIFIC not given, setting to c(0,3,4) for backward compatibility", call.=FALSE) }
  AFNIheaderpart("ORIENT_SPECIFIC",header$ORIENT_SPECIFIC, conhead)
  header$ORIENT_SPECIFIC <- NULL
  if (!is.null(origin)) header$ORIGIN <- origin
  if (is.null(header$ORIGIN)) { header$ORIGIN <- c(0,0,0); warning("ORIGIN not given, setting to c(0,0,0) for backward compatibility", call.=FALSE) }
  AFNIheaderpart("ORIGIN",header$ORIGIN, conhead)
  header$ORIGIN <- NULL
  if (!is.null(delta)) header$DELTA <- delta
  if (is.null(header$DELTA)) { header$DELTA <- c(4,4,4); warning("DELTA not given, setting to c(4,4,4) for backward compatibility", call.=FALSE) }
  AFNIheaderpart("DELTA",header$DELTA, conhead)
  header$DELTA <- NULL

  # write mandatory attributes for time series
  if (taxis) {
    if (is.null(header$TAXIS_NUMS)) stop("TAXIS_NUMS not given")
    AFNIheaderpart("TAXIS_NUMS",header$TAXIS_NUMS, conhead)
    if (is.null(header$TAXIS_FLOATS)) stop("TAXIS_FLOATS not given")
    AFNIheaderpart("TAXIS_FLOATS",header$TAXIS_FLOATS, conhead)
    header$TAXIS_FLOATS <- NULL
    if (!is.null(header$TAXIS_NUMS[2])) if (header$TAXIS_NUMS[2] != 0) if (is.null(header$TAXIS_OFFSETS)) {
      stop("TAXIS_OFFSETS not given", call.=FALSE) } else { AFNIheaderpart("TAXIS_OFFSETS",header$TAXIS_OFFSETS, conhead); header$TAXIS_OFFSETS <- NULL }
    header$TAXIS_NUMS <- NULL
  } 

  # almost mandatory attributes
  if (!is.null(idcode)) header$IDCODE_STRING <- idcode
  if (is.null(header$IDCODE_STRING)) {header$IDCODE_STRING <- "WIAS_noid"; warning("IDCODE_STRING not given, setting to WIAS_noid", call.=FALSE) }
  AFNIheaderpart("IDCODE_STRING",header$IDCODE_STRING, conhead)
  header$IDCODE_STRING <- NULL
  header$IDCODE_DATE <- date()
  AFNIheaderpart("IDCODE_DATE",header$IDCODE_DATE, conhead)
  header$IDCODE_DATE <- NULL
  if (is.null(header$BYTEORDER_STRING)) {
    endian <- .Platform$endian
    header$BYTEORDER_STRING <- switch(endian, "little" = "LSB_FIRST",
                                       "big" = "MSB_FIRST")
  } else {
    if (!(header$BYTEORDER_STRING%in% c("MSB_FIRST","LSB_FIRST"))) header$BYTEORDER <- "MSB_FIRST"
    endian <- switch(header$BYTEORDER_STRING, "MSB_FIRST" = "big",
                                       "LSB_FIRST" = "little")
  }
  AFNIheaderpart("BYTEORDER_STRING",header$BYTEORDER_STRING, conhead,c("MSB_FIRST","LSB_FIRST"))
  header$BYTEORDER_STRING <- NULL
  header$BRICK_STATS <- apply(ttt,4,range)
  dim(header$BRICK_STATS) <- NULL
  AFNIheaderpart("BRICK_STATS",header$BRICK_STATS, conhead)
  header$BRICK_FLOAT_FACS <- rep(0,dim(ttt)[4])
  if (is.null(header$BRICK_TYPES)) {
    warning("no BRICK_TYPES given, assuming short", call.=FALSE)
    header$BRICK_TYPES <- rep(1,dim(ttt)[4])
  }
  bricktypes <- header$BRICK_TYPES[1]
  if ((bricktypes == 1) && (max(abs(header$BRICK_STATS)) > 32767)) { 
    for (k in 1:dim(ttt)[4]) {
      header$BRICK_FLOAT_FACS[k] <- max(abs(header$BRICK_STATS[2*k-1]),abs(header$BRICK_STATS[2*k]))/32767
      ttt[,,,k] <- ttt[,,,k] / header$BRICK_FLOAT_FACS[k]
    }
  }
  header$BRICK_STATS <- NULL
  AFNIheaderpart("BRICK_TYPES",header$BRICK_TYPES, conhead)
  header$BRICK_TYPES <- NULL
  AFNIheaderpart("BRICK_FLOAT_FACS",header$BRICK_FLOAT_FACS, conhead)
  header$BRICK_FLOAT_FACS <- NULL
  if (!is.null(label)) header$BRICK_LABS <- paste(label,collapse="~")
  if (!is.null(header$BRICK_LABS)) AFNIheaderpart("BRICK_LABS",header$BRICK_LABS, conhead)
  header$BRICK_LABS <- NULL

  # note attributes
  if (!(is.null(note))) header$HISTORY_NOTE <- paste(header$HISTORY_NOTE,note)
  if (is.null(header$HISTORY_NOTE)) header$HISTORY_NOTE <- ""
  AFNIheaderpart("HISTORY_NOTE",header$HISTORY_NOTE, conhead)
  header$HISTORY_NOTE <- NULL

  for (name in names(header)) {
    AFNIheaderpart(name,header[[name]], conhead)
    header[[name]] <- NULL
  }
  close(conhead)

  if (!(bricktypes %in% c(1,3,5))) stop("Sorry, cannot write this BRICK_TYPES.", call.=FALSE)
  conbrik <- file(paste(filename, ".BRIK", sep=""), "wb")
  dim(ttt) <- NULL
  switch(bricktypes, "1" = writeBin(as.integer(ttt), conbrik, size=2, endian=endian),
                     "3" = writeBin(as.numeric(ttt), conbrik, size=4, endian=endian),
                     "5" = writeBin(as.complex(ttt), conbrik, size=16, endian=endian))
  close(conbrik)
}


read.DICOM <- function(filename,includedata=TRUE) {
  read.DICOM.groupelement <- function(con,endian="little") {
    if (endian == "little") {
      paste(paste(rev(readBin(con,"raw",2,1)),collapse=""),paste(rev(readBin(con,"raw",2,1)),collapse=""),sep=",")
    } else {
      paste(paste(readBin(con,"raw",2,1),collapse=""),paste(readBin(con,"raw",2,1),collapse=""),sep=",")
    }    
  }
  
  con <- file(filename,"rb")

  endian <- "little"
  
  headerdetails <- list()
  bytes <- 0

  empty <- paste(readBin(con,"character",128,1),collapse="")
  bytes <- bytes + 128
  if (empty != "") {
    warning("This does not seem to be a DICOM file\n")
  }
  dicom <- readChar(con,4)
  bytes <- bytes + 4
  if (dicom != "DICM") {
    warning("This does not seem to be a DICOM file\n")
  }
  
  groupelement <- 0
  while (TRUE) {
    groupelement <- read.DICOM.groupelement(con)
    bytes <- bytes + 4
    if (groupelement == "7fe0,0010") break
    vr <- readChar(con,2)
    bytes <- bytes + 2

    if (vr %in% c("OB","OW","OF","SQ","UT","UN")) {
      reserved <- readBin(con,"raw",2)
      bytes <- bytes + 2
      length <- readBin(con,"integer",1,4,signed=FALSE,endian=endian)
      bytes <- bytes + 4
      if (length == -1) {
        while (TRUE) {
          groupelement <- read.DICOM.groupelement(con)
          bytes <- bytes + 4
          if (groupelement == "fffe,e0dd") break
          length <- readBin(con,"integer",1,4,signed=FALSE,endian=endian)
          bytes <- bytes + 4
          value <- readBin(con,"raw",length,1)
          bytes <- bytes + length
        }
        length <- readBin(con,"integer",1,4,signed=FALSE,endian=endian)
        bytes <- bytes + 4
        value <- "undecoded sequence"
      }
    } else {
      length <- readBin(con,"integer",1,2,signed=FALSE,endian=endian)
      bytes <- bytes + 2
    }

    if (length != 0) {
      if (vr %in% c("UI","DS","SH","IS")) {
        value <- readChar(con,length)      
        bytes <- bytes + length
      } else if (vr %in% c("US")) {
        total <- 0
        value <- ""
        while (total < length) {
          value <- paste(value,readBin(con,"integer",1,2,signed=FALSE,endian=endian),sep="")
          bytes <- bytes + 2
          total <- total + 2
        }
      } else {
        value <- readBin(con,"raw",length,1)
        bytes <- bytes + length
      }
    }
#    cat(groupelement,vr,length,paste(value,collapse=""),"\n")
    if (groupelement %in% c("0002,0010",
                            "0018,0050",
                            "0020,0010",
                            "0020,0011",
                            "0020,0012",
                            "0020,0013",
                            "0020,0032",
                            "0020,1041",
                            "0028,0002",
                            "0028,0010",
                            "0028,0011",
                            "0028,0030",
                            "0028,0100",
                            "0028,1050",
                            "0028,1051",
                            "0028,1052",
                            "0028,1053")) {
      headerdetails[[groupelement]] <- value
    }
  }

  # belongs to last groupelement "7fe0,0010"
  vr <- readChar(con,2)
  bytes <- bytes + 2
  if (vr %in% c("OB","OW","OF","SQ","UT","UN")) {
    reserved <- readBin(con,"raw",2)
    bytes <- bytes + 2
    length <- readBin(con,"integer",1,4,signed=FALSE,endian=endian)
    bytes <- bytes + 4
  } else {
    length <- readBin(con,"integer",1,2,signed=FALSE,endian=endian)
    bytes <- bytes + 2
  }
  headerdetails[[groupelement]] <- length
#  cat("Bytes for Header:",bytes,"\n")

  close(con)
  
  if (!is.null(headerdetails[["0028,0010"]])) {
    xdim <- as.integer(headerdetails[["0028,0010"]])
  } else {
    xdim <- NULL
  }
  if (!is.null(headerdetails[["0028,0011"]])) {
    ydim <- as.integer(headerdetails[["0028,0011"]])
  } else {
    ydim <- NULL
  }
  if (!is.null(headerdetails[["0028,0100"]])) {
    if (headerdetails[["0028,0100"]] == 16) {
      depth <- 2
    } else {
      depth <- 1
    }
  } else {
    depth <- 1
  }
  if (!is.null(headerdetails[["0028,0030"]])) {
    delta <- as.numeric(strsplit(headerdetails[["0028,0030"]],"\\",fixed=TRUE)[[1]])
    if (!is.null(headerdetails[["0018,0050"]])) {
      delta <- c(delta,as.numeric(headerdetails[["0018,0050"]]))
    }
  } else {
    delta <- NULL
  }
  if (!is.null(headerdetails[["0020,0011"]])) {
    series <- as.integer(headerdetails[["0020,0011"]])
  } else {
    series <- NULL
  }
  if (!is.null(headerdetails[["0020,0013"]])) {
    image <- as.integer(headerdetails[["0020,0013"]])
  } else {
    image <- NULL
  }

  # read again 
  con <- file(filename,"rb")

  header <- readBin(con,"raw",bytes)
  
  if (includedata) {
    if (is.null(xdim) || is.null(ydim)) {
      ttt <- readBin(con,"integer",length/depth,depth,signed=FALSE,endian=endian)
      warning("Cannot assign dimension to image because information was not found!")
    } else {
      ttt <- array(readBin(con,"integer",length/depth,depth,signed=FALSE,endian=endian),c(xdim,ydim))
    }
    
    z <- list(header=header,ttt=writeBin(as.numeric(ttt),raw(),4),format="DICOM",delta=delta,series=series,image=image,dim=c(xdim,ydim))
  } else {
    z <- list(header=header,format="DICOM",delta=delta,series=series,image=image,dim=c(xdim,ydim))
  }
  close(con)
  class(z) <- "fmridata"

  attr(z,"file") <- filename
  invisible(z)
}

read.NIFTI.header <- function(con) {
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
  header$diminfo <- readChar(con,1)
  header$dimension <- readBin(con,"int",8,2,endian=endian)
  header$intentp1 <- readBin(con,"double",1,4,endian=endian)
  header$intentp2 <- readBin(con,"double",1,4,endian=endian)
  header$intentp3 <- readBin(con,"double",1,4,endian=endian)
  header$intentcode <- readBin(con,"int",1,2,endian=endian)
  header$datatype <- readBin(con,"int",1,2,endian=endian)
  header$bitpix <- readBin(con,"int",1,2,endian=endian)
  header$slicestart <- readBin(con,"int",1,2,endian=endian)
  header$pixdim <- readBin(con,"double",8,4,endian=endian)
  header$voxoffset <- readBin(con,"double",1,4,endian=endian)
  header$sclslope <- readBin(con,"double",1,4,endian=endian)
  header$sclinter <- readBin(con,"double",1,4,endian=endian)
  header$sliceend <- readBin(con,"int",1,2,endian=endian)
  header$slicecode <- readChar(con,1)
  header$xyztunits <- readChar(con,1)
  header$calmax <- readBin(con,"double",1,4,endian=endian)
  header$calmin <- readBin(con,"double",1,4,endian=endian)
  header$sliceduration <- readBin(con,"double",1,4,endian=endian)
  header$toffset <- readBin(con,"double",1,4,endian=endian)
  header$glmax <- readBin(con,"int",1,4,endian=endian)
  header$glmin <- readBin(con,"int",1,4,endian=endian)
  header$describ <- readChar(con,80)
  header$auxfile <- readChar(con,24)
  header$qform <- readBin(con,"int",1,2,endian=endian)
  header$sform <- readBin(con,"int",1,2,endian=endian)
  header$quaternb <- readBin(con,"double",1,4,endian=endian)
  header$quaternc <- readBin(con,"double",1,4,endian=endian)
  header$quaternd <- readBin(con,"double",1,4,endian=endian)
  header$qoffsetx <- readBin(con,"double",1,4,endian=endian)
  header$qoffsety <- readBin(con,"double",1,4,endian=endian)
  header$qoffsetz <- readBin(con,"double",1,4,endian=endian)
  header$srowx <- readBin(con,"double",4,4,endian=endian)
  header$srowy <- readBin(con,"double",4,4,endian=endian)
  header$srowz <- readBin(con,"double",4,4,endian=endian)
  header$intentname <- readChar(con,16)
  header$magic <- readChar(con,4)
  
  header
}

read.NIFTI <- function(filename) {
  fileparts <- strsplit(filename,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])

  if (ext == "nii") {
    filename.nii <- filename
    filename.hdr <- paste(c(fileparts[-length(fileparts)],"hdr"),collapse=".")
    filename.img <- paste(c(fileparts[-length(fileparts)],"img"),collapse=".")
  } else if (ext == "hdr") {
    filename.hdr <- filename
    filename.img <- paste(c(fileparts[-length(fileparts)],"img"),collapse=".")
  } else if (ext == "img") {
    filename.hdr <- paste(c(fileparts[-length(fileparts)],"hdr"),collapse=".")
    filename.img <- filename
  } else {
    filename.nii <- paste(filename,".nii",sep="")
    filename.hdr <- paste(filename,".hdr",sep="")
    filename.img <- paste(filename,".img",sep="")
  }

  if ((ext != "hdr") && (ext != "img") && (!is.na(file.info(filename.nii)$size))) {
    con <- file(filename.nii,"rb")
    header <- read.NIFTI.header(con)
    if (!(header$magic == "n+1") && !(header$magic == "ni1")) 
      warning("Hmmm! Dont see the magic NIFTI string! Try to proceed, but maybe some weird results will occur!");
    bytes <- header$voxoffset - 348
    header$extension <- readBin(con,"raw",bytes)
  } else {
    if (is.na(file.info(filename.hdr)$size) | (file.info(filename.hdr)$size < 348))
      stop("Hmmm! This does not seem to be a NIFTI header (hdr/img-pair)! Wrong size or does not exist!");
    con <- file(filename.hdr,"rb")
    header <- read.NIFTI.header(con)
    header$extension <- NULL  
    close(con)
    if (is.na(file.info(filename.img)$size))     
      stop("Hmmm! This does not seem to be a NIFTI header (hdr/img-pair)! img-file not found!");
    con <- file(filename.img,"rb")
  }
    
  dx <- header$dimension[2]
  dy <- header$dimension[3]
  dz <- header$dimension[4]
  dt <- header$dimension[5]
  endian <- header$endian
  if (header$datatype == 1) { # logical
    what <- "raw"
    signed <- TRUE
    size <- 1
  } else if (header$datatype == 2) { # unsigned char????
    what <- "int"
    signed <- FALSE
    size <- if (header$bitpix) header$bitpix/8 else 2
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
  ttt <- readBin(con, what, n=dx*dy*dz*dt*size, size=size, signed=signed, endian=endian) 
  close(con)

  if (min(abs(header$pixdim[2:4])) != 0) {
    weights <-
      abs(header$pixdim[2:4]/min(abs(header$pixdim[2:4])))
  } else {
    weights <- NULL
  }
  dim(ttt) <- c(dx,dy,dz,dt)

  mask <- array(TRUE,c(dx,dy,dz))
  mask[ttt[,,,1] < quantile(ttt[,,,1],0.75)] <- FALSE

  z <- list(ttt=writeBin(as.numeric(ttt),raw(),4),format="NIFTI",delta=header$pixdim[2:4],
                origin=NULL,orient=NULL,dim=header$dimension[2:5],weights=weights,header=header,mask=mask)

  class(z) <- "fmridata"

  invisible(z)
}

extract.data <- function(z,what="data") {
  if (!("fmridata"%in%class(z))) {
    warning("extract.data: data not of class <fmridata>. Try to proceed but strange things may happen")
  }
  if (what=="residuals") {  
      if(!is.null(z$resscale)){
          ttt <- readBin(z$res,"integer",prod(z$dim),2)*z$resscale 
          dim(ttt) <- z$dim
          } else {
          warning("extract.data: No residuals available, returning NULL")
          ttt <- NULL
      }
      } else { 
      ttt <- readBin(z$ttt,"numeric",prod(z$dim),4)
      dim(ttt) <- z$dim
      }
 
  invisible(ttt)
}
