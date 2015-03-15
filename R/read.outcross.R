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
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to read data from input file


##' Read data from a segregating full-sib population
##' 
##' Imports data from a full-sib family derived from the cross of two outbred
##' parents and creates an object of class \code{outcross}.
##' 
##' The file format is quite similar to that used by \code{MAPMAKER/EXP}
##' (\cite{Lincoln et al.}, 1993). The first line contains two integers: the
##' number of individuals and the number of markers.
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
##' Finally, after the segregation type comes the genotype data for the
##' corresponding marker. Depending on the segregation type, genotypes may be
##' denoted by \code{ac}, \code{ad}, \code{bc}, \code{bd}, \code{a}, \code{ba},
##' \code{b}, \code{bc}, \code{ab} and \code{o}, in several possible
##' combinations. To make things easier, we have followed \strong{exactly} the
##' notation used by \cite{Wu et al.} (2002). Genotypes \emph{must} be
##' separated by commas. Missing values are denoted by \code{"-"}
##' 
##' The \code{example} directory in the package distribution contains an
##' example data file to be read with this function. Further instructions can
##' be found at the tutorial distributed along with this package.
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
##' \code{"D2"}} \item{input}{the name of the input file.}
##' @author Adapted from Karl Broman (package \pkg{qtl}) by Gabriel R A
##' Margarido, \email{gramarga@@gmail.com}
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
read.outcross <- 
function (dir, file) {
  # create file name
  if (missing(file)) 
    stop("missing file")
  if(!missing(dir) && dir != "")
    file <- file.path(dir, file)
  
  # count lines in rawfile
  n.lines <- length(scan(file, what = character(), skip = 0, 
                         nlines = 0, blank.lines.skip = FALSE,
                         quiet = TRUE, sep = "\n"))
  
  # begin reading the genotype data
  cur.mar <- 0
  flag <- 0
  for (i in 1:n.lines) {
    a <- scan(file, what = character(), skip = i - 1, 
              nlines = 1, blank.lines.skip = TRUE, quiet = TRUE)
    if (length(a) == 0) next
    if (length(grep("#", a[1])) != 0) next
    if (flag == 0) {
      # reading first line (basic informatioxn about the data)
      flag <- 1
      n.ind <- as.numeric(a[1])
      n.mar <- as.numeric(a[2])
      cat(" Working...\n\n")
      marnames <- rep("", n.mar)
      geno <- matrix(0, ncol = n.mar, nrow = n.ind)
      segr.type <- character(n.mar)
    }
    else {
      if (substring(a[1], 1, 1) == "*") {
        # new marker
        cur.mar <- cur.mar + 1
        cur.row <- 1
        marnames[cur.mar] <- substring(a[1], 2)
        if (length(a) < 2) {
          stop("the segregation type of marker ", marnames[cur.mar],
               " should be placed next to its name (on the same line)")
        }
        segr.type[cur.mar] <- a[2]
        if (length(a) > 2) {
          # reading genotypes on the line where marker name is
          g <- paste(a[c(-1,-2)], collapse = "")
          g <- unlist(strsplit(g, ","))
          n <- length(g)
          geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
        }
        else n <- 0
        cur.row <- cur.row + n
      }
      else {
        # continuation lines
        g <- paste(a, collapse = "")
        g <- unlist(strsplit(g, ","))
        n <- length(g)
        geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
        cur.row <- cur.row + n
      }  # end continuation lines
    }  # end non-intro line
  }
  # done reading the raw file

  # add marker names to data
  colnames(geno) <- marnames
  # changes -'s (missing data) to NA's
  geno[!is.na(geno) & geno == "-"] <- NA
  
  # recoding data
  temp.data <- codif.data(geno,segr.type)
  geno <- temp.data[[1]]
  segr.type.num <- temp.data[[2]]
  rm(temp.data)
  
  cat(" --Read the following data:\n")
  cat("\tNumber of individuals: ", n.ind, "\n")
  cat("\tNumber of markers:     ", n.mar, "\n")
  
  structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,segr.type = segr.type,
                 segr.type.num=segr.type.num, input=file), class = "outcross")
}



# print method for object class 'outcross'
print.outcross <-
function(x,...) {
  # checking for correct object
  if (any(is.na(match(c("geno", "n.ind", "n.mar", "segr.type"),
                      names(x))))) 
    stop("this is not an object of class 'outcross'")

  # printing brief summary of the data
  cat("This is an object of class 'outcross'\n")
  cat("    No. individuals:   ", x$n.ind, "\n")
  cat("    No. markers:       ", x$n.mar, "\n")
  cat("    Segregation types:\n")
  
  # counting the number of markers with each segregation type
  quant <- cbind(table(x$segr.type))
  for(i in 1:length(quant)){
    cat(paste("       ", rownames(quant)[i], ":\t", quant[i],
              "\n", sep=""))
  }
}

# end of file
