#########################################################################
##                                                                      #
##  Package: onemap                                                     #
##                                                                      #
##  File: read_mapmaker.R                                               #
##  Contains: read_mapmaker                                             #
##                                                                      #
##  Written by Marcelo Mollinari with minor modifications by Cristiane  #
##  Taniguti                                                            #
##  Adapted from read.cross.mm (found in the R package qtl)             #
##  copyright (c) 2000-6, Karl W Broman                                 #
##                                                                      #
##  First version: 09/27/2009                                           #
##  License: GNU General Public License version 3 (June, 2007) or later #
##                                                                      #
#########################################################################

## Function to read data in MAPMAKER style from input file


##' Read data from a Mapmaker raw file
##'
##' Imports data from a Mapmaker raw file.
##'
##' For details about MAPMAKER files see \cite{Lincoln et al.} (1993). The
##' current version supports backcross, F2s and RIL populations. The file
##' can contain phenotypic data, but it will not be used in the analysis.
##'
##' @param dir directory where the input file is located.
##' @param file the name of the input file which contains the data to be read.
##' @param verbose A logical, if TRUE it output progress status
##' information.
##' 
##' @return An object of class \code{onemap}, i.e., a list with the following
##' components: \item{geno}{a matrix with integers indicating the genotypes
##' read for each marker in \code{onemap} fashion. Each column contains data
##' for a marker and each row represents an individual.}
##'
##' \code{MAPMAKER/EXP} fashion, i.e., 1, 2, 3: AA, AB, BB, respectively; 3, 4:
##' BB, not BB, respectively; 1, 5: AA, not AA, respectively. Each column
##' contains data for a marker and each row represents an individual.
##'
##' \item{n.ind}{number of individuals.} \item{n.mar}{number of markers.}
##' \item{segr.type}{a vector with the segregation type of each marker, as
##' \code{strings}. Segregation types were adapted from outcross segregation
##' types, using the same notation. For details see \link{read_onemap}.}
##' \item{segr.type.num}{a vector with the segregation type of each marker,
##' represented in a simplified manner as integers. Segregation types were
##' adapted from outcross segregation types. For details see
##' \link{read_onemap}.} \item{input}{the name of the input file.}
##' \item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
##' values.  Each column contains data for a trait and each row represents an
##' individual. Currently ignored.} \item{error}{matrix containing HMM emission probabilities}
##' 
##' @author Adapted from Karl Broman (package \pkg{qtl}) by Marcelo Mollinari,
##' \email{mmollina@@usp.br}
##' @seealso \code{mapmaker_example_bc} and \code{mapmaker_example_f2} raw files in the
##' package source.
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing genetic
##' linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and reference
##' manual. \emph{A Whitehead Institute for Biomedical Research Technical
##' Report}.
##' @keywords IO
##' @examples
##' \donttest{
##'  map_data <-read_mapmaker(file=system.file("extdata/mapmaker_example_f2.raw", package = "onemap"))
##'  #Checking 'mapmaker_example_f2'
##'  data(mapmaker_example_f2)
##'  names(mapmaker_example_f2)
##' }
##'@export
read_mapmaker<-function (file=NULL, dir=NULL, verbose=TRUE)
{
    ## create file name
    if (is.null(file))
        stop("Missing file.")
    if (!is.null(dir) && dir != "") {
        file <- file.path(dir, file)
    }
    
    ## count lines in rawfile
    n.lines <- length(scan(file, what = character(), skip = 0,
                           nlines = 0, blank.lines.skip = FALSE,
                           quiet = TRUE, sep = "\n"))
    ## begin reading/parsing the genotype data
    cur.mar <- 0
    cur.phe <- 0
    NEW.symb <- c("1", "2", "3", "4", "5", NA)
    OLD.symb <- c("A", "H", "B", "D", "C", "-")
    flag <- 0
    for (i in 1:n.lines)
    {
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
                type <- "backcross"
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
            if(verbose){
                cat(" --Read the following data:\n")
                
                cat("\tType of cross:         ", type, "\n")
                cat("\tNumber of individuals: ", n.ind, "\n")
                cat("\tNumber of markers:     ", n.mar, "\n")
            }
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
    if(type=="f2"){
        cl <- c("onemap", "f2")
        ##checking for markers with one class (e.g A A A - - - A - A - - - A)
        ##they are not necessarily monomorphic because we don't know the missing data
        mkt.mono<-NULL
        mkt.mono<-which(apply(geno, 2, function(x) sum(!is.na(unique(x))))<=1)
        if(length(mkt.mono)!=0){
            segr.type[mkt.mono]<-"A.H.B"
        }
        mkt<-apply(geno, 2, function(x) prod(unique(x), na.rm=TRUE))
        segr.type[mkt==2 | mkt==3 | mkt==6]<-"A.H.B"
        segr.type[mkt==12]<-"D.B"
        segr.type[mkt==5]<-"C.A"
        mkt.rest<-which(segr.type=="")
        for(i in mkt.rest)
        {
            if(any(is.na(match(na.omit(unique(geno[,i])), 1:5))))
            {
                mkt.wrg.names <- paste(sQuote(colnames(geno)[i]), collapse = ", ")
                msg <- sprintf(ngettext(length(mkt.wrg.names),
                                        "marker %s has invalid codification",
                                        "markers %s have invalid codification"), mkt.wrg.names)
                stop(msg)
            }
            else
                segr.type[i]<-"M.X"
        }
        segr.type.num[segr.type=="A.H.B"]<-4
        segr.type.num[segr.type=="C.A"]<-5
        segr.type.num[segr.type=="D.B"]<-5
        segr.type.num[segr.type=="M.X"]<-0
        # Adapting to change in f2 HMM == out HMM
        geno[, segr.type=="C.A"][which(geno[, segr.type=="C.A"] == 1)] <- 2 
        geno[, segr.type=="C.A"][which(geno[, segr.type=="C.A"] == 5)] <- 1 
        geno[, segr.type=="D.B"][which(geno[, segr.type=="D.B"] == 4)] <- 1 
        geno[, segr.type=="D.B"][which(geno[, segr.type=="D.B"] == 3)] <- 2 
        geno[is.na(geno)]<-0
    }
    else if(type=="backcross"){
        cl <- c("onemap", "backcross")
        ##Verifying if there are up to two classes in bc data, ignoring NAs
        if(sum(!is.na(unique(as.vector(geno)))) > 2)
            stop("check data: there are more than 2 classes for backcross")
        segr.type[]<-"A.H"
        segr.type.num<-rep(8,ncol(geno))
        geno[is.na(geno)]<-0
        geno[geno==3]<-1 #coding for raw data entered as H and B
    }
    else if(type=="riself" || type=="risib"){
        if (type=="riself") cl <- c("onemap", "riself")
        else cl <- c("onemap", "risib")
        ##Verifying if there are up to two classes in ril data, ignoring NAs
        if(sum(!is.na(unique(as.vector(geno)))) > 2)
            stop("check data: there are more than 2 classes for ", type)
        segr.type[]<-"A.B"
        segr.type.num<-rep(9,ncol(geno))
        geno[is.na(geno)]<-0
        #geno[geno==3]<-2 #coding as backcross
    }
    else
        stop("Invalid cross type")
    if(n.phe != 0) {
        miss.value.pheno <- apply((apply(pheno, 2,is.na)),2,sum)
        if(verbose){
            cat("\tMissing trait values:      ", "\n")
            for(i in 1:n.phe) {
                cat("\t",formatC(paste(colnames(pheno)[i],":",sep=""),width=max(nchar(paste(colnames(pheno),":",sep="")))), miss.value.pheno[i], "\n")
            }
        }
    }
    
    onemap.obj <- list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                       segr.type = segr.type, segr.type.num=segr.type.num,
                       input=file, n.phe=n.phe, pheno = pheno)
    class(onemap.obj) <- cl
    new.onemap.obj <- create_probs(onemap.obj, global_error = 10^-5)
    
    structure(new.onemap.obj)
}

## end of file
