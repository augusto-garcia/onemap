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
