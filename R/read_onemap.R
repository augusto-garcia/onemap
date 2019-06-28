#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: read_onemap.R                                                 #
# Contains: read_onemap, print.onemap                                 #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2015, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/25/2015                                           #
# Last update: 01/26/2017                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Read data from all types of progenies supported by OneMap
##'
##' Imports data derived from outbred parents (full-sib family) or inbred
##' parents (backcross, F2 intercross and recombinant inbred lines obtained
##' by self- or sib-mating). Creates an object of class \code{onemap}.
##'
##' The file format is similar to that used by \code{MAPMAKER/EXP}
##' (\cite{Lincoln et al.}, 1993). The first line indicates the cross type
##' and is structured as \code{data type \{cross\}}, where \code{cross}
##' must be one of \code{"outcross"}, \code{"f2 intercross"},
##' \code{"f2 backcross"}, \code{"ri self"} or  \code{"ri sib"}. The second line
##' contains five integers: i) the number of individuals; ii) the number of
##' markers; iii) an indicator variable taking the value 1 if there is CHROM
##' information, i.e., if markers are anchored on any reference sequence, and
##' 0 otherwise; iv) a similar 1/0 variable indicating whether there is POS
##' information for markers; and v) the number of phenotypic traits.
##'
##' The next line contains sample IDs, separated by empty spaces or tabs.
##' Addition of this sample ID requirement makes it possible for separate input
##' datasets to be merged.
##'
##' Next comes the genotype data for all markers. Each new marker is initiated
##' with a \dQuote{*} (without the quotes) followed by the marker name, without
##' any space between them. Each marker name is followed by the corresponding
##' segregation type, which may be: \code{"A.1"}, \code{"A.2"}, \code{"A.3"},
##' \code{"A.4"}, \code{"B1.5"}, \code{"B2.6"}, \code{"B3.7"}, \code{"C.8"},
##' \code{"D1.9"}, \code{"D1.10"}, \code{"D1.11"}, \code{"D1.12"},
##' \code{"D1.13"}, \code{"D2.14"}, \code{"D2.15"}, \code{"D2.16"},
##' \code{"D2.17"} or \code{"D2.18"} (without quotes), for full-sibs [see
##' \code{\link[onemap]{marker_type}} and \cite{Wu et al.} (2002) for details].
##' Other cross types have special marker types: \code{"A.H"} for backcrosses;
##' \code{"A.H.B"} for F2 intercrosses; and \code{"A.B"} for recombinant inbred
##' lines.
##'
##' After the segregation type comes the genotype data for the
##' corresponding marker. Depending on the segregation type, genotypes may be
##' denoted by \code{ac}, \code{ad}, \code{bc}, \code{bd}, \code{a}, \code{ba},
##' \code{b}, \code{bc}, \code{ab} and \code{o}, in several possible
##' combinations. To make things easier, we have followed \strong{exactly} the
##' notation used by \cite{Wu et al.} (2002). Allowed values for backcrosses
##' are \code{a} and \code{ab}; for F2 crosses they are \code{a}, \code{ab} and
##' \code{b}; for RILs they may be \code{a} and \code{b}. Genotypes \emph{must}
##' be separated by a space. Missing values are denoted by \code{"-"}.
##'
##' If there is physical information for markers, i.e., if they are anchored at
##' specific positions in reference sequences (usually chromosomes), this is
##' included immediately after the marker data. These lines start with special
##' keywords \code{*CHROM} and \code{*POS} and contain \code{strings} and
##' \code{integers}, respectively, indicating the reference sequence and
##' position for each marker. These also need to be separated by spaces.
##'
##' Finally, if there is phenotypic data, it will be added just after the marker
##' or \code{CHROM}/\code{POS} data. They need to be separated by spaces as
##' well, using the same symbol for missing information.
##'
##' The \code{example} directory in the package distribution contains an
##' example data file to be read with this function. Further instructions can
##' be found at the tutorial distributed along with this package.
##'
##' @param dir directory where the input file is located.
##' @param inputfile the name of the input file which contains the data to be read.
##' @return An object of class \code{onemap}, i.e., a list with the following
##' components: \item{geno}{a matrix with integers indicating the genotypes
##' read for each marker. Each column contains data for a marker and each row
##' represents an individual.} \item{n.ind}{number of individuals.}
##' \item{n.mar}{number of markers.} \item{segr.type}{a vector with the
##' segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
##' vector with the segregation type of each marker, represented in a
##' simplified manner as integers, i.e. 1 corresponds to markers of type
##' \code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
##' to markers of type \code{"B2.6"}; 4 corresponds to markers of type
##' \code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
##' to markers of type \code{"D1"} and 7 corresponds to markers of type
##' \code{"D2"}. Markers for F2 intercrosses are coded as 1; all other crosses
##' are left as \code{NA}.} \item{input}{the name of the input file.}
##' \item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
##' values. Each column contains data for a trait and each row represents an
##' individual.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{combine_onemap}} and the \code{example}
##' directory in the package source.
##' @references Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993)
##' Constructing genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial
##' and reference manual. \emph{A Whitehead Institute for Biomedical Research
##' Technical Report}.
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords IO
##' @examples
##'
##'   \dontrun{
##'     outcr_data <- read_onemap(dir="work_directory", inputfile="data_file.txt")
##'   }
##'@export
read_onemap <- function (inputfile=NULL, dir=NULL) {
  if (is.null(inputfile)){
     stop("missing file")
  }
  if (!is.null(inputfile) && !is.null(dir)) {
     inputfile <- file.path(dir, inputfile)
  }

  f <- file(inputfile, open = "r")
  on.exit(close(f))

  ## Read cross type information
  l <- scan(f, what=character(), nlines = 1,
            blank.lines.skip = TRUE, quiet = TRUE)
  if ((length(l) != 3 && length(l) != 4)  || l[1] != "data" || l[2] != "type") {
    stop("The first line of the input file must conform to: 'data type X'",
         call.= TRUE)
  }
  crosstype <- match.arg(l[3],
                         c("outcross", "f2", "ri"),
                         several.ok= FALSE)
  if (crosstype == "f2") {
    if (length(l) == 3 || (l[4] != "intercross" && l[4] != "backcross")) {
      stop("Unknown cross type.")
    }
    ## "f2" denotes an F2 intercross; "backcross" denotes an F2 backcross
    if (l[4] == "backcross") {
      crosstype <- "backcross"
    }
  }
  else if (crosstype == "ri") {
    if (length(l) == 3 || (l[4] != "self" && l[4] != "sib")) {
      stop("Unknown cross type.")
    }
    ## "riself" denotes RI by selfing; "risib" denotes RI by sib mating
    if (l[4] == "self") {
      crosstype <- "riself"
    }
    else if (l[4] == "sib") {
      crosstype <- "risib"
    }
  }

  ## Read in the number of individuals, markers, CHROM/POS availability and number of phenotypes
  l <- scan(f, what=integer(), nlines = 1,
            blank.lines.skip = TRUE, quiet = TRUE)
  if (length(l) != 5) {
    stop("The second line of the input file must have the following information: 'number of individuals', 'number of markers', 'presence of CHROM data', 'presence of POS data' and 'number of traits'. These numbers must be separated with an empty space.",
         call.= TRUE)
  }
  n.ind <- l[1]
  n.mar <- l[2]
  has_CHROM <- l[3] == 1
  has_POS <- l[4] == 1
  n.phe <- l[5]

  ## Parse the sample IDs
  l <- scan(f, what=character(), nlines = 1,
            blank.lines.skip = TRUE, quiet = TRUE)
  if (length(l) != n.ind) {
    stop("Incomplete or extra sample ID information.", call. = TRUE)
  }
  sample_IDs <- l

  ## Read marker genotype information
  cat(" Working...\n\n")
  l <- matrix(scan(f, what = character(), nlines = n.mar,
                   blank.lines.skip = TRUE, quiet = TRUE),
              n.ind + 2, n.mar)
  if (length(l) != (2 + n.ind) * n.mar) {
    stop("Incomplete or extra genotype information.", call. = TRUE)
  }
  if (length(unique(l[1,])) != length(l[1,])) {
    stop("There are markers with the same name.", call. = TRUE)
  }

  ## Get marker names
  bad_lines <- which(substr(l[1,], 1, 1) != "*")
  if (length(bad_lines)) {
    stop("Lines with genotype information must begin with '*[markername]' and contain one field per individual.", call. = TRUE)
  }
  marnames <- substring(l[1,], 2)

  ## Get marker types
  segr.type <- l[2,]

  ## Get genotype matrix
  geno <- l[-c(1,2),]
  geno[!is.na(geno) & geno == "-"] <- NA
  colnames(geno) <- marnames
  rownames(geno) <- sample_IDs
  colnames(geno) <- marnames
  rownames(geno) <- sample_IDs

  temp.data <- codif_data(geno, segr.type, crosstype)
  geno <- temp.data[[1]]
  segr.type.num <- temp.data[[2]]
  rm(temp.data)

  ## Get CHROM/POS
  if (has_CHROM) {
    l <- scan(f, what = character(), nlines = 1,
              blank.lines.skip = TRUE, quiet = TRUE)
    if (length(l) != n.mar + 1) {
      stop("Incomplete or extra CHROM information.", call. = TRUE)
    }
    if (l[1] != "*CHROM") {
      stop("CHROM information expected but not found.", call. = TRUE)
    }
    l[l == "-"] <- NA
    CHROM <- l[-1]
  }
  else {
    CHROM <- NULL
  }
  if (has_POS) {
    l <- scan(f, what = character(), nlines = 1,
              blank.lines.skip = TRUE, quiet = TRUE)
    if (length(l) != n.mar + 1) {
      stop("Incomplete or extra POS information.", call. = TRUE)
    }
    if (l[1] != "*POS") {
      stop("POS information expected but not found.", call. = TRUE)
    }
    l[l == "-"] <- NA
    ## Check remaining (non-missing) fields for non-integer values
    temp <- as.integer(l[!is.na(l)][-1])
    if (any(is.na(temp))) {
      stop("POS line can only contain integers.", call. = TRUE)
    }
    POS <- as.integer(l[-1])
  }
  else {
    POS <- NULL
  }

  ## Get phenotype data
  if (n.phe) {
    l <- matrix(scan(f, what = character(), nlines = n.phe,
                     blank.lines.skip = TRUE, quiet = TRUE),
                n.ind + 1, n.phe)
    if (length(l) != (1 + n.ind) * n.phe) {
      stop("Incomplete or extra phenotype information.", call. = TRUE)
    }

    bad_lines <- which(substr(l[1,], 1, 1) != "*")
    if (length(bad_lines)) {
      stop("Lines with phenotype information must begin with '*[phenoname]' and contain one field per individual.", call. = TRUE)
    }
    pheno_names <- substring(l[1,], 2)

    pheno <- as.matrix(l[-1,])
    pheno[!is.na(pheno) & pheno == "-"] <- NA
    mode(pheno) <- "numeric"
    colnames(pheno) <- pheno_names
  }
  else {
    pheno <- NULL
  }

  ## Output
  cat(" --Read the following data:\n")
  cat("\tType of cross:          ", crosstype, "\n")
  cat("\tNumber of individuals:  ", n.ind, "\n")
  cat("\tNumber of markers:      ", n.mar, "\n")
  cat("\tChromosome information: ", ifelse(is.null(CHROM), "no", "yes"), "\n")
  cat("\tPosition information:   ", ifelse(is.null(POS), "no", "yes"), "\n")
  cat("\tNumber of traits:       ", n.phe, "\n")
  if(n.phe != 0) {
    miss.value.pheno <- apply((apply(pheno, 2,is.na)),2,sum)
    cat("\tMissing trait values:      ", "\n")
    for(i in 1:n.phe) {
      cat("\t",formatC(paste(colnames(pheno)[i],":",sep=""),width=max(nchar(paste(colnames(pheno),":",sep="")))), miss.value.pheno[i], "\n")
    }
  }

  ## Return "onemap" object
  onemap.obj <- structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                 segr.type = segr.type, segr.type.num = segr.type.num,
                 n.phe = n.phe, pheno = pheno, CHROM = CHROM, POS = POS,
                 input = inputfile),
            class = c("onemap", crosstype))
  new.onemap.obj <- create_probs(onemap.obj)
  return(new.onemap.obj)
}

## Print method for object class 'onemap'
##'@export
##' @method print onemap
print.onemap <- function (x, ...) {
  ## Print a brief summary of the data
  not_miss <- 100*sum(x$geno!=0)/length(x$geno)
  cat("  This is an object of class 'onemap'\n")
  cat("    Type of cross:     ", class(x)[2], "\n")
  cat("    No. individuals:   ", x$n.ind, "\n")
  cat("    No. markers:       ", x$n.mar, "\n")
  cat("    CHROM information: ", ifelse(is.null(x$CHROM), "no", "yes"), "\n")
  cat("    POS information:   ", ifelse(is.null(x$POS), "no", "yes"), "\n")
  cat("    Percent genotyped: ", round(not_miss), "\n\n")

  ## Count the number of markers with each segregation type
  cat("    Segregation types:\n")
  quant <- table(x$segr.type)
  ## F2 intercross
  names(quant)[which(names(quant) == "A.H.B") ] <- "AA : AB : BB -->"
  names(quant)[which(names(quant) == "M.X")]    <- "AA : AB : BB (+ dom)  -->"
  names(quant)[which(names(quant) == "D.B")]    <- " Not BB : BB -->"
  names(quant)[which(names(quant) == "C.A")]    <- " Not AA : AA -->"
  ## Backcross
  names(quant)[which(names(quant) == "A.H")]    <- "     AA : AB -->"
  ## RILs
  names(quant)[which(names(quant) == "A.B")]    <- "     AA : BB -->"
  ## Outcross
  names(quant)[which(names(quant) == "A.1")]    <- "         A.1 -->"
  names(quant)[which(names(quant) == "A.2")]    <- "         A.2 -->"
  names(quant)[which(names(quant) == "A.3")]    <- "         A.3 -->"
  names(quant)[which(names(quant) == "A.4")]    <- "         A.4 -->"
  names(quant)[which(names(quant) == "B1.5")]   <- "        B1.5 -->"
  names(quant)[which(names(quant) == "B2.6")]   <- "        B2.6 -->"
  names(quant)[which(names(quant) == "B3.7")]   <- "        B3.7 -->"
  names(quant)[which(names(quant) == "C.8")]    <- "         C.8 -->"
  names(quant)[which(names(quant) == "D1.9")]   <- "        D1.9 -->"
  names(quant)[which(names(quant) == "D1.10")]  <- "       D1.10 -->"
  names(quant)[which(names(quant) == "D1.11")]  <- "       D1.11 -->"
  names(quant)[which(names(quant) == "D1.12")]  <- "       D1.12 -->"
  names(quant)[which(names(quant) == "D1.13")]  <- "       D1.13 -->"
  names(quant)[which(names(quant) == "D2.14")]  <- "       D2.14 -->"
  names(quant)[which(names(quant) == "D2.15")]  <- "       D2.15 -->"
  names(quant)[which(names(quant) == "D2.16")]  <- "       D2.16 -->"
  names(quant)[which(names(quant) == "D2.17")]  <- "       D2.17 -->"
  names(quant)[which(names(quant) == "D2.18")]  <- "       D2.18 -->"

  for (i in 1:length(quant)) {
    cat(paste("       ", names(quant)[i], "  ", quant[i],
              "\n", sep = ""))
  }

  ## Check for phenotypic data
  cat("\n    No. traits:        ", x$n.phe, "\n")
  if(x$n.phe > 0) {
    miss.value <- apply((apply(x$pheno, 2,is.na)),2,sum)
    cat("    Missing trait values:", "\n")
    for (i in 1:x$n.phe) {
      cat("\t",formatC(paste(colnames(x$pheno)[i],":", sep=""), width= max(nchar(paste(colnames(x$pheno),":",sep="")))), miss.value[i], "\n")
    }
  }
}
## end of file
