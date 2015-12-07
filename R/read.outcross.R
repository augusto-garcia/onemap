#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: read.outcross.R                                               #
# Contains: read.outcross, print.outcross                             #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# Adapted from read.cross.mm (package: R/qtl)                         #
# copyright (c) 2000-6, Karl W Broman                                 #
# First version: 11/07/2007                                           #
#                                                                     #
# On August 29th, 2015, it was modified by Augusto Garcia, by changing#
# the code adding a new feature developed by Luciano da Costa e Silva #
# on his package oneqtl. Help files were also modified.               #
# The modification allows the inclusion of phenotypic traits.         #
#                                                                     #
# Last update: 2015/12/07                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Read data from a full-sib progeny (outcrossing populations)
##'
##' Imports data from a full-sib family derived from the cross of two outbred
##' parents and creates an object of class \code{outcross}.
##'
##' The file format is quite similar to that used by \code{MAPMAKER/EXP}
##' (\cite{Lincoln et al.}, 1993). The first line contains three integers: the
##' number of individuals, the number of markers and the number of traits.
##'
##' Next comes the genotype data for all markers. Each new marker is initiated
##' with a \dQuote{*} (without the quotes) followed by the marker name, without
##' any space between them. Each marker name is followed by the corresponding
##' segregation type, which may be: \code{"A.1"}, \code{"A.2"}, \code{"A.3"},
##' \code{"A.4"}, \code{"B1.5"}, \code{"B2.6"}, \code{"B3.7"}, \code{"C.8"},
##' \code{"D1.9"}, \code{"D1.10"}, \code{"D1.11"}, \code{"D1.12"},
##' \code{"D1.13"}, \code{"D2.14"}, \code{"D2.15"}, \code{"D2.16"},
##' \code{"D2.17"} or \code{"D2.18"} (without quotes) [see
##' \code{\link[onemap]{marker.type}} and \cite{Wu et al.} (2002) for details].
##'
##' After the segregation type comes the genotype data for the
##' corresponding marker. Depending on the segregation type, genotypes may be
##' denoted by \code{ac}, \code{ad}, \code{bc}, \code{bd}, \code{a}, \code{ba},
##' \code{b}, \code{bc}, \code{ab} and \code{o}, in several possible
##' combinations. To make things easier, we have followed \strong{exactly} the
##' notation used by \cite{Wu et al.} (2002). Genotypes \emph{must} be
##' separated by commas. Missing values are denoted by \code{"-"}
##'
##' Finally, if there is phenotypic data, it will be added just after the marker data.
##' They need to be separated by commas as well, using the same symbol for missing information.
##'
##' The \code{example} directory in the package distribution contains an
##' example data file to be read with this function. Further instructions can
##' be found at the tutorial distributed along with the package.
##'
##' @param dir directory where the input file is located.
##' @param file the name of the input file which contains the data to be read.
##' @return An object of class \code{outcross}, i.e., a list with the following
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
##' \code{"D2"}} \item{n.phen}{the number of traits included in the file}
##' \item{pheno}{the name of the phenoytpes} \item{input}{the name of the input file.}
##' @author Adapted from Karl Broman (package \pkg{qtl}) by Gabriel R A
##' Margarido, \email{gramarga@@gmail.com}, later with additions from Luciano C Silva
##' @seealso \code{example} directory in the package source.
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing genetic
##' linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and reference
##' manual. \emph{A Whitehead Institute for Biomedical Research Technical
##' Report}.
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords IO
##' @examples
##'
##'   \dontrun{
##'     outcr_data <-
##' read.outcross(dir="work_directory",file="data_file.txt")
##'   }
##'
read.outcross <- function (dir, file) {
  if (missing(file))
    stop("missing file")
  if (!missing(dir) && dir != "")
    file <- file.path(dir, file)
  n.lines <- length(scan(file, what = character(), skip = 0,
                         nlines = 0, blank.lines.skip = FALSE, quiet = TRUE, sep = "\n"))
  cur.mar <- 0
  cur.pheno <- 0
  n.phen <- 0
  flag <- 0
  for (i in 1:n.lines) {
    a <- scan(file, what = character(), skip = i - 1, nlines = 1,
              blank.lines.skip = TRUE, quiet = TRUE)
    if (length(a) == 0)
      next
    if (length(grep("#", a[1])) != 0)
    next
    if (flag == 0) { #reading first line
      flag <- 1
      if (length(a) != 3)
        stop("The first line of the input file must have the following information: 'number of individuals', 'number of markers', and 'number of traits'. These numbers must be separated with an empty space. For instance, 10 5 0.", call.= TRUE)
      n.ind <- as.numeric(a[1])
      n.mar <- as.numeric(a[2])
      n.phen <- as.numeric(a[3]) #
      cat(" Working...\n\n")
      marnames <- rep("", n.mar)
      geno <- matrix(0, ncol = n.mar, nrow = n.ind)
      segr.type <- character(n.mar)
      if (n.phen == 0) {
        pheno <- numeric(0) #matrix(1:n.ind, ncol = 1)
        phenonames <- character(0) #c("number")
      }
      else {
        pheno <- matrix(0, ncol = n.phen, nrow = n.ind)
        phenonames <- rep("", n.phen)
      }
    } #finishes reading first line in the file (flag==0)
    else {#reading lines of markers and traits
      if (substring(a[1], 1, 1) == "*") {#reading lines of markers and traits that start with "*"
        cur.mar <- cur.mar + 1
        cur.row <- 1
        if (cur.mar > n.mar) {#reading lines of traits that start with "*"
          cur.pheno <- cur.pheno + 1
          if (cur.pheno > n.phen)
            next
          phenonames[cur.pheno] <- substring(a[1], 2)
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
                themessage <- paste(themessage, " for pheno \"",
                                    phenonames[cur.pheno], "\" were", sep = "")
              }
              else {
                themessage <- paste("The value \"", droppedasmissing,
                                    "\" ", sep = "")
                themessage <- paste(themessage, " for pheno \"",
                                    phenonames[cur.pheno], "\" was", sep = "")
              }
              themessage <- paste(themessage, "interpreted as missing.")
              warning(themessage, call.= FALSE)
            }
            pheno[cur.row + (0:(n - 1)), cur.pheno] <- numerp
          }
          else n <- 0
          cur.row <- cur.row + n
        }#finishes reading lines of traits that start with "*" (cur.mar > n.mar)
        else {#reading lines of markers that start with "*" (cur.mar <= n.mar)
          marnames[cur.mar] <- substring(a[1], 2)
          if (length(a) < 2) {
            stop("the segregation type of marker ", marnames[cur.mar],
                 " should be placed next to its name (on the same line)")
          }
          segr.type[cur.mar] <- a[2]
          if (length(a) > 2) {
            g <- paste(a[c(-1, -2)], collapse = "")
            g <- unlist(strsplit(g, ","))
            n <- length(g)
            geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
          }
          else n <- 0
          cur.row <- cur.row + n
        }#finishes reading lines of markers that start with "*" (cur.mar <= n.mar)
      }#finishes reading lines of markers and traits that start with "*"
      else {#reading lines of markers and traits that do not start with "*"
        if (cur.mar > n.mar) {
          a[a == "-"] <- NA
          n <- length(a)
          pheno[cur.row + (0:(n - 1)), cur.pheno] <- as.numeric(a)
          cur.row <- cur.row + n
        }#finishes reading lines of traits that do not start with "*"
        else {
          g <- paste(a, collapse = "")
          g <- unlist(strsplit(g, ","))
          n <- length(g)
          geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
          cur.row <- cur.row + n
        }#finishes reading lines of markers that do not start with "*"
      }#finishes reading lines of markers and traits that do not start with "*"
    }#finishes reading lines of markers and traits
  } #close loop for lines (i)
  colnames(geno) <- marnames
  geno[!is.na(geno) & geno == "-"] <- NA
  temp.data <- codif.data(geno, segr.type)
  geno <- temp.data[[1]]
  segr.type.num <- temp.data[[2]]
  rm(temp.data)
  if(n.phen != 0) {
    colnames(pheno) <- phenonames
    pheno[!is.na(pheno) & pheno == "-"] <- NA
  }
  cat(" --Read the following data:\n")
  cat("\tNumber of individuals: ", n.ind, "\n")
  cat("\tNumber of markers:     ", n.mar, "\n")
  cat("\tNumber of traits:      ", n.phen, "\n")
  if(n.phen != 0) {
    miss.value.pheno <- apply((apply(pheno, 2,is.na)),2,sum)
    cat("\tMissing trait values:      ", "\n")
    for(i in 1:n.phen) {
      cat("\t",formatC(paste(colnames(pheno)[i],":",sep=""),width=max(nchar(paste(colnames(pheno),":",sep="")))), miss.value.pheno[i], "\n")
    }
  }

  structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                 segr.type = segr.type, segr.type.num = segr.type.num,
                 n.phen = n.phen, pheno = pheno,
                 input = file),
            class = "outcross")
}

print.outcross <- function (x, ...) {
  if (any(is.na(match(c("geno", "n.ind", "n.mar", "segr.type"),
                      names(x)))))
    stop("this is not an object of class 'outcross'")
  cat("  This is an object of class 'outcross'\n")
  cat("    No. individuals:   ", x$n.ind, "\n")
  cat("    No. markers:       ", x$n.mar, "\n")
  cat("    Segregation types:\n")
  quant <- cbind(table(x$segr.type))
  for (i in 1:length(quant)) {
    cat(paste("       ", rownames(quant)[i], ":\t", quant[i],
              "\n", sep = ""))
  }
  cat("    No. traits:        ", x$n.phen, "\n")
  if(x$n.phen > 0) {
    miss.value <- apply((apply(x$pheno, 2,is.na)),2,sum)
    cat("    Missing trait values:", "\n")
    for (i in 1:x$n.phen) {
      #cat(paste("       ", colnames(x$pheno)[i], "\n", sep = ""))
      cat("\t",formatC(paste(colnames(x$pheno)[i],":", sep=""), width= max(nchar(paste(colnames(x$pheno),":",sep="")))), miss.value[i], "\n")
    }
  }
}
# end of file
