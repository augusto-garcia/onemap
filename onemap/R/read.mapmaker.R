#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: read.mapmaker.R                                               #
# Contains: read.mapmaker, print.f2.onemap, print.bc.onemap           #
#           print.riself.onemap, print.risib.onemap                   #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# Adapted from read.cross.mm (found in the R package qtl)             #
# copyright (c) 2000-6, Karl W Broman                                 #
#                                                                     #
# First version: 09/27/2009                                           #
# Last update:   03/05/2011                                           #
# License: GNU General Public License version 3 (June, 2007) or later #
#                                                                     #
#######################################################################

## Function to read data in MAPMAKER style from input file
read.mapmaker<-function (dir, file) 
{
  ## create file name
  if (missing(file)) 
    stop("Missing file.")
  if (!missing(dir) && dir != "") {
    file <- file.path(dir, file)
  }
  # count lines in rawfile
  n.lines <- length(scan(file, what = character(), skip = 0, 
                         nlines = 0, blank.lines.skip = FALSE,
                         quiet = TRUE, sep = "\n"))
  ## begin reading/parsing the genotype data
  cur.mar <- 0
  cur.phe <- 0
  NEW.symb <- c("1", "2", "3", "4", "5", NA)
  OLD.symb <- c("A", "H", "B", "D", "C", "-")
  flag <- 0
  for (i in 1:n.lines) {
    a <- scan(file, what = character(), skip = i - 1, 
              nlines = 1, blank.lines.skip = TRUE, quiet = TRUE)
    if (length(a) == 0) 
      next
    if (length(grep("#", a[1])) != 0) 
    next
    if (flag == 0) {
      flag <- 1
      if (!is.na(match("intercross", a))) 
        type <- "f2"
      else if (!is.na(match("backcross", a))) 
        type <- "bc"
      else if (!is.na(match("self", a))) 
        type <- "riself"
      else if (!is.na(match("sib", a))) 
        type <- "risib"
      else stop("File indicates invalid cross type: ", 
                a[length(a)], ".")
    }
    else if (flag == 1) {
      flag <- 2
      n.ind <- as.numeric(a[1])
      n.mar <- as.numeric(a[2])
      n.phe <- as.numeric(a[3])
      cat(" --Read the following data:\n")
      cat("\tType of cross:         ", type, "\n")
      cat("\tNumber of individuals: ", n.ind, "\n")
      cat("\tNumber of markers:     ", n.mar, "\n")
      ## if there's a set of "symbols" for non-standard symbols in
      ## the file, use them.
      if (length(a) > 3 && ("symbols" %in% a)) {
        o <- match("symbols", a)
        b <- a[-(1:o)]
        infile.symb <- substring(b, 1, 1)
        std.symb <- substring(b, 3, 3)
        wh <- rep(0, length(std.symb))
        fixed <- rep(0, length(OLD.symb))
        for (j in 1:length(std.symb)) if (std.symb[j] %in% 
                                          OLD.symb) 
          wh[j] <- match(std.symb[j], OLD.symb)
        for (j in 1:length(std.symb)) if (wh[j] != 0) {
          OLD.symb[wh[j]] <- infile.symb[j]
          fixed[wh[j]] <- 1
        }
        temp <- table(OLD.symb)
        if (any(temp > 1)) {
          for (j in names(temp)[temp > 1]) {
            o <- OLD.symb == j & fixed == 0
            if (any(o)) 
              OLD.symb[o] <- paste(OLD.symb[o], "   ")
          }
        }
      }
      marnames <- rep("", n.mar)
      geno <- matrix(0, ncol = n.mar, nrow = n.ind)
      if (n.phe == 0) {
        pheno <- matrix(1:n.ind, ncol = 1)
        phenames <- c("number")
      }
      else {
        pheno <- matrix(0, ncol = n.phe, nrow = n.ind)
        phenames <- rep("", n.phe)
      }
    }
    else {
      if (substring(a[1], 1, 1) == "*") {
        cur.mar <- cur.mar + 1
        cur.row <- 1
        if (cur.mar > n.mar) { ## now reading phenotypes
          cur.phe <- cur.phe + 1
          if (cur.phe > n.phe) 
            next
          phenames[cur.phe] <- substring(a[1], 2)
          if (length(a) > 1) {
            p <- a[-1]
            p[p == "-"] <- NA
            n <- length(p)
            oldna <- is.na(p)
            numerp <- suppressWarnings(as.numeric(p))
            newna <- is.na(numerp)
            wh <- !oldna & newna
            if (any(wh)) {
              droppedasmissing <- unique(p[wh])
              if (length(droppedasmissing) > 1) {
                themessage <- paste("The values", paste("\"", 
                                                        droppedasmissing, "\"", sep = "", collapse = " "))
                themessage <- paste(themessage, " for phenotype \"", 
                                    phenames[cur.phe], "\" were", sep = "")
              }
              else {
                themessage <- paste("The value \"", droppedasmissing, 
                                    "\" ", sep = "")
                themessage <- paste(themessage, " for phenotype \"", 
                                    phenames[cur.phe], "\" was", sep = "")
              }
              themessage <- paste(themessage, "interpreted as missing.")
              warning(themessage)
            }
            pheno[cur.row + (0:(n - 1)), cur.phe] <- numerp
          }
          else n <- 0
          cur.row <- cur.row + n
        }
        else { ## reading genotypes
          marnames[cur.mar] <- substring(a[1], 2)
          if (length(a) > 1) {
            g <- paste(a[-1], collapse = "")
            h <- g <- unlist(strsplit(g, ""))
            for (j in seq(along = NEW.symb)) {
              if (any(h == OLD.symb[j])) 
                g[h == OLD.symb[j]] <- NEW.symb[j]
            }
            n <- length(g)
            if(any(is.na(match(g,NEW.symb))))
              stop(paste("Invalid marker codification. Please check data for marker", marnames[cur.mar]), ".", sep="")           
            geno[cur.row + (0:(n - 1)), cur.mar] <- as.numeric(g)
          }
          else n <- 0
          cur.row <- cur.row + n
        }
      }
      else { ## continuation lines
        if (cur.mar > n.mar) { ## now reading phenotypes
          a[a == "-"] <- NA
          n <- length(a)
          pheno[cur.row + (0:(n - 1)), cur.phe] <- as.numeric(a)
          cur.row <- cur.row + n
        }
        else {
          g <- paste(a, collapse = "")
          h <- g <- unlist(strsplit(g, ""))
          for (j in seq(along = NEW.symb)) {
            if (any(h == OLD.symb[j])) 
              g[h == OLD.symb[j]] <- NEW.symb[j]
          }
          n <- length(g)
          geno[cur.row + (0:(n - 1)), cur.mar] <- as.numeric(g)
          cur.row <- cur.row + n
        }
      }## end continuation line
    }## end non-intro line
  }
  dimnames(geno) <- list(NULL, marnames)
  dimnames(pheno) <- list(NULL, phenames)
  ## done reading the raw file
  
  ## data coding in onemap style
  segr.type<-character(n.mar)
  segr.type.num<-numeric(n.mar)
  ph<-rep(1,n.mar)
  ## maintains the data in MAPMAKER style (maybe it will be removed in a near future)
  geno.mmk<-list(geno=geno, type=type)
  geno.mmk[is.na(geno.mmk)]<-0

  if(type=="f2"){
    cl<-"f2.onemap"
    ##Verifying if there are up to three classes in f2 data, ignoring NAs
    if(any(unlist(lapply(apply(geno, 2, table),length)) > 3))
      stop("check data: there are more than 3 classes for f2")

    ##checking for markers with one class (e.g A A A - - - A - A - - - A)
    ##they are not necessarily monomorphic because we don't know the missing data
    mkt.mono<-NULL
    mkt.mono<-which(apply(geno, 2, function(x) sum(!is.na(unique(x))))<=1)
    if(length(mkt.mono)!=0){
      segr.type[mkt.mono]<-"B3.7"
   #   mkt.mono.names <- paste(sQuote(colnames(geno)[mkt.mono]), collapse = ", ")
   #   msg <- sprintf(ngettext(length(mkt.mono),
   #                           "There is one marker with one class on dataset: %s",
   #                           "There are markers with one class on dataset: %s"), mkt.mono.names)
   #    warning(msg, domain=NA)
    }

    ##more data coding in onemap style (f2 is equivalent to B3.7 and C.8)
    mkt<-apply(geno, 2, function(x) prod(unique(x), na.rm=TRUE))
    segr.type[mkt==2 | mkt==3 | mkt==6]<-"B3.7"
    segr.type[mkt==5 | mkt==12]<-"C.8"
    mkt.wrg<-NULL
    mkt.wrg<-which(segr.type=="")
    if(length(mkt.wrg)!=0){
      mkt.wrg.names <- paste(sQuote(colnames(geno)[mkt.wrg]), collapse = ", ")
      msg <- sprintf(ngettext(length(mkt.wrg),
                              "marker %s has invalid codification",
                              "markers %s have invalid codification"), mkt.wrg.names)
       stop(msg)
    }
    segr.type.num[segr.type=="B3.7"]<-4
    segr.type.num[segr.type=="C.8"]<-5
    ##coding phases: A B H -> 1; C A -> -1; B D -> 1
    ##ph[mkt==2]<-ph[mkt==3]<-ph[mkt==6]<-ph[mkt==12]<-1
    ph[mkt==5]<--1
    ##data codind for markers that segregates in 3:1 fashion (type C.8)
    geno[is.na(geno)]<-0    
    geno[,segr.type.num==5][geno[,segr.type.num==5]==1]<-2
    geno[,segr.type.num==5][geno[,segr.type.num==5]==5]<-1
    geno[,segr.type.num==5][geno [,segr.type.num==5]==4]<-1
    geno[,segr.type.num==5][geno[,segr.type.num==5]==3]<-2
  }
  else if(type=="bc"){
    cl<-"bc.onemap"
    ##Verifying if there are up to two classes in bc data, ignoring NAs
    if(sum(!is.na(unique(as.vector(geno)))) > 2)
      stop("check data: there are more than 2 classes for bc")
    segr.type[]<-"D1.10" ##data coding in onemap style (bc is equivalent to D1.10)
    segr.type.num<-rep(6,ncol(geno))
    ph<-rep(1,n.mar) #phases: all merkers in coupling
    geno[is.na(geno)]<-0 
    geno[geno==3]<-1 #coding for raw data entered as H and B 
  }
  else if(type=="riself" || type=="risib"){
    if (type=="riself") cl<-"riself.onemap"
    else cl<-"risib.onemap"
    ##Verifying if there are up to two classes in ril data, ignoring NAs
    if(sum(!is.na(unique(as.vector(geno)))) > 2)
      stop("check data: there are more than 2 classes for ", type)
    segr.type[]<-"D1.10" ##more data coding in onemap style 
    segr.type.num<-rep(6,ncol(geno))
    ph<-rep(1,n.mar) #phases: all merkers inn coupling
    geno[is.na(geno)]<-0 
    geno[geno==3]<-2 #coding as backcross
  }
  else
    stop("Invalide type of cross")
  structure(list(geno = geno, geno.mmk = geno.mmk, n.ind = n.ind, n.mar = n.mar,
                 segr.type = segr.type, segr.type.num=segr.type.num,phase=ph,
                 input=file, n.phe=n.phe, pheno = pheno),  class = cl)
}

##print method for object class 'f2.onemap'
print.f2.onemap<-function (x, ...){
  ##checking for correct object
  if (any(is.na(match(c("geno","geno.mmk","n.ind","n.mar","segr.type",
                        "segr.type.num","phase","input","n.phe","pheno"), 
                      names(x))))) 
    stop("this is not an object of class 'f2.onemap'")
  mis<-100*sum(x$geno!=0)/length(x$geno) #counting missing data
  ##printing brief summary of the data
  cat("This is an object of class 'f2.onemap'\n")
  cat("    No. individuals:    ", x$n.ind, "\n")
  cat("    No. markers:        ", x$n.mar, "\n")
  cat("    Percent genotyped:  ", round(mis), "\n\n")
  cat("    Number of markers per type:\n")
  ##counting the number of markers with each segregation type
  quant <- table(x$segr.type.num-x$phase)
  names(quant)[which(names(quant)=="3") ] <-"AA : AB : BB -->"
  names(quant)[which(names(quant)=="5") ] <-"AA : AB : BB -->"
  names(quant)[which(names(quant)=="4")]  <-" Not BB : BB -->"
  names(quant)[which(names(quant)=="6")]  <-" Not AA : AA -->" 
  for (i in 1:length(quant)) {
    cat(paste("       ", names(quant)[i],"  ", quant[i], 
              "\n", sep = ""))
  }
  ##checking for phenotipic data
  if(x$n.phe==0) cat("\nThis data contains no phenotypic information\n")
  else cat("\nThis data contains phenotypic information\n")
}

##print method for object class 'bc.onemap'
print.bc.onemap<-function (x, ...) {
  ##checking for correct object
  if (any(is.na(match(c("geno","geno.mmk","n.ind","n.mar","segr.type",
                        "segr.type.num","phase","input","n.phe","pheno"), 
                      names(x))))) 
    stop("this is not an object of class 'bc.onemap'")
  mis<-100*sum(x$geno!=0)/length(x$geno)#counting missing data
  ##printing brief summary of the data
  cat("This is an object of class 'bc.onemap'\n")
  cat("    No. individuals:    ", x$n.ind, "\n")
  cat("    No. markers:        ", x$n.mar, "\n")
  cat("    Percent genotyped:  ", round(mis), "\n\n")
  cat("    Number of markers per type:\n")
  cat(paste("       AA : AB --> ", ncol(x$geno), " marker(s)\n", sep = ""))
  ##checking for phenotipic data
  if(x$n.phe==0) cat("\nThis data contains no phenotypic information\n\n")
  else cat("\nThis data contains phenotypic information\n\n")
}

##print method for object class 'riself.onemap'
print.riself.onemap<-function (x, ...) {
  ##checking for correct object
  if (any(is.na(match(c("geno","geno.mmk","n.ind","n.mar","segr.type",
                        "segr.type.num","phase","input","n.phe","pheno"), 
                      names(x))))) 
    stop("this is not an object of class 'riself.onemap'")
  mis<-100*sum(x$geno!=0)/length(x$geno)#counting missing data
  ##printing brief summary of the data
  cat("This is an object of class 'riself.onemap'\n")
  cat("    No. individuals:    ", x$n.ind, "\n")
  cat("    No. markers:        ", x$n.mar, "\n")
  cat("    Percent genotyped:  ", round(mis), "\n\n")
  cat("    Number of markers per type:\n")
  cat(paste("       AA : BB --> ", ncol(x$geno), " marker(s)\n", sep = ""))
  ##checking for phenotipic data
  if(x$n.phe==0) cat("\nThis data contains no phenotypic information\n\n")
  else cat("\nThis data contains phenotypic information\n\n")
}

##print method for object class 'risib.onemap'
print.risib.onemap<-function (x, ...) {
  ##checking for correct object
  if (any(is.na(match(c("geno","geno.mmk","n.ind","n.mar","segr.type",
                        "segr.type.num","phase","input","n.phe","pheno"), 
                      names(x))))) 
    stop("this is not an object of class 'risib.onemap'")
  mis<-100*sum(x$geno!=0)/length(x$geno)#counting missing data
  ##printing brief summary of the data
  cat("This is an object of class 'risib.onemap'\n")
  cat("    No. individuals:    ", x$n.ind, "\n")
  cat("    No. markers:        ", x$n.mar, "\n")
  cat("    Percent genotyped:  ", round(mis), "\n\n")
  cat("    Number of markers per type:\n")
  cat(paste("       AA : BB --> ", ncol(x$geno), " marker(s)\n", sep = ""))
  ##checking for phenotipic data
  if(x$n.phe==0) cat("\nThis data contains no phenotypic information\n\n")
  else cat("\nThis data contains phenotypic information\n\n")
}


## end of file
