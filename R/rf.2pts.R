#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: rf.2pts.R                                                     ##
## Contains: rf.2pts, print.rf.2pts                                    ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  ##
## copyright (c) 2007-15, Gabriel R A Margarido and Marcelo Mollinari  ##
##                                                                     ##
## First version: 11/07/2007                                           ##
## Last update: 12/03/2015                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     #
#######################################################################

## Function to perform two-point analyses for all markers in a data set


##' Two-point analysis between genetic markers
##' 
##' Performs the two-point (pairwise) analysis proposed by \cite{Wu et al.
##' (2002)} between all pairs of markers.
##' 
##' For \code{n} markers, there are \deqn{\frac{n(n-1)}{2}}{n*(n-1)/2} pairs of
##' markers to be analyzed. Therefore, completion of the two-point analyses can
##' take a long time.
##' 
##' @aliases rf.2pts print.rf.2pts
##' @param input.obj an object of class \code{outcross}, \code{bc.onemap},
##' \code{f2.onemap}, \code{riself.onemap} or \code{risib.onemap}.
##' @param LOD minimum LOD Score to declare linkage (defaults to \code{3}).
##' @param max.rf maximum recombination fraction to declare linkage (defaults
##' to \code{0.50}).
##' @param verbose logical. If \code{TRUE}, current progress is shown; if
##' \code{FALSE}, no output is produced.
##' @param x an object of class \code{rf.2pts}.
##' @param mrk1,mrk2 optionally, two markers can be specified. If so, detailed
##' results of the two-point analysis will be printed for this pair. Both
##' arguments can be numeric or character strings indicating the numbers/names
##' corresponding to any markers in the input file.
##' @param \dots further arguments, passed to other methods. Currently ignored.
##' @return An object of class \code{rf.2pts}, which is a list containing the
##' following components:  \item{n.mar}{total number of markers.} \item{LOD}{minimum LOD Score to declare
##' linkage.} \item{max.rf}{maximum recombination fraction to declare linkage.}
##' \item{input}{the name of the input file.} \item{analysis}{an array with the
##' complete results of the two-point analysis for each pair of markers.} 
##' @note The thresholds used for \code{LOD} and \code{max.rf} will be used in
##' subsequent analyses, but can be overriden.
##' @author Gabriel R A Margarido \email{gramarga@gmail.com} and Marcelo Mollinari \email{mmollina@usp.br}
##' @references Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002)
##' Simultaneous maximum likelihood estimation of linkage and linkage phases in
##' outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords utilities
##' @examples
##' 
##'   data(example.out)
##' 
##'   twopts <- rf.2pts(example.out,LOD=3,max.rf=0.5) # perform two-point analyses
##'   twopts
##' 
##'   print(twopts,"M1","M2") # detailed results for markers 1 and 2
##' 
rf.2pts <- function(input.obj, LOD=3, max.rf=0.50, verbose = TRUE) {
    ## checking for correct object
    if(!any(class(input.obj)=="outcross"||class(input.obj)=="f2.onemap" || class(input.obj)=="bc.onemap" || class(input.obj)=="riself.onemap" || class(input.obj)=="risib.onemap")) stop(deparse(substitute(input.obj))," is not an object of class 'outcross', 'bc.onemap', 'f2.onemap', 'riself.onemap' or 'risib.onemap'")
    if (input.obj$n.mar<2) stop("there must be at least two markers to proceed with analysis")
    ## creating variables (result storage and progress output)
    if(class(input.obj)=="outcross")
        r<-est_rf_out(geno = input.obj$geno, seg_type = input.obj$segr.type.num, nind = input.obj$n.ind, verbose = verbose)
    else if(class(input.obj)=="f2.onemap")
        r<-est_rf_f2(geno = input.obj$geno, seg_type = input.obj$segr.type.num, nind = input.obj$n.ind, verbose = verbose)
    else if(class(input.obj)=="bc.onemap")
        r<-est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=0, verbose = verbose)
    else if(class(input.obj)=="riself.onemap")
        r<-est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=1, verbose = verbose)
    else if(class(input.obj)=="risib.onemap")
        r<-est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=2, verbose = verbose)
    structure(list(data.name=as.character(sys.call())[2], n.mar=input.obj$n.mar, LOD=LOD, max.rf=max.rf, input=input.obj$input, analysis=r), class = c("rf.2pts", class(input.obj)))
}

## print method for object class 'rf.2pts'
print.rf.2pts <- function(x, mrk1=NULL, mrk2=NULL,...) {
    ## checking for correct object
    if(!any(class(x)=="rf.2pts")) stop(deparse(substitute(x))," is not an object of class 'rf.2pts'")
    
    if (is.null(mrk1) || is.null(mrk2)) {
        ## printing a brief summary
        cat("  This is an object of class 'rf.2pts'\n")
        cat("\n  Criteria: LOD =", x$LOD, ", Maximum recombination fraction =",
            x$max.rf, "\n")
        cat("\n  This object is too complex to print\n")
        cat("  Type 'print(object,mrk1=marker,mrk2=marker)' to see the analysis for two markers\n")
        cat("    mrk1 and mrk2 can be the names or numbers of both markers\n")
    }
  else {
      ## printing detailed results for two markers
      ## checking if markers exist and converting character to numeric
      if(any(class(x)=="f2.onemap") || any(class(x)=="bc.onemap") || any(class(x)=="risib.onemap") || any(class(x)=="riself.onemap"))
      {
          if (is.character(mrk1) && is.character(mrk2)) {
              mrk1name<-mrk1
              mrk2name<-mrk2
              mrk1<-match(mrk1, colnames(x$analysis))
              mrk2<-match(mrk2, colnames(x$analysis))
              if (is.na(mrk1)) stop("marker ", mrk1name, " not found")
              if (is.na(mrk2)) stop("marker ", mrk2name, " not found")
          }
          else if(is.numeric(mrk1) && is.numeric(mrk2))
          {
              if(mrk1 > nrow(x$analysis)) stop("marker ", mrk1, " not found: marker number out of bounds")
              if(mrk2 > nrow(x$analysis)) stop("marker ", mrk2, " not found: marker number out of bounds")
          }
          else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
          cat("  Results of the 2-point analysis for markers:", colnames(x$analysis)[mrk1],
              "and", colnames(x$analysis)[mrk2], "\n")
          
          ## results found
          cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
              x$max.rf, "\n\n")
          ## do not print anything for the same marker
          if(mrk1 == mrk2) stop("mrk1 and mrk2 are the same")
          ## results of the two-point analysis
          if (mrk1 > mrk2){
              r<-x$analysis[mrk1,mrk2]
              LOD<-x$analysis[mrk2,mrk1]
              output<-c(r, LOD)
              names(output)<-c("rf","LOD")
              print(output)
          }
          else
          {
              r<-x$analysis[mrk2,mrk1]
              LOD<-x$analysis[mrk1,mrk2]
              output<-c(r, LOD)
              names(output)<-c("rf","LOD")
              print(output)
          }
      }
    else if(any(class(x)=="outcross"))
    {
        if (is.character(mrk1) && is.character(mrk2)) {
            mrk1name<-mrk1
            mrk2name<-mrk2
            mrk1<-match(mrk1, colnames(x$analysis[[1]]))
            mrk2<-match(mrk2, colnames(x$analysis[[1]]))
            if (is.na(mrk1)) stop("marker ", mrk1name, " not found")
            if (is.na(mrk2)) stop("marker ", mrk2name, " not found")
        }
        else if(is.numeric(mrk1) && is.numeric(mrk2))
        {
            if(mrk1 > nrow(x$analysis[[1]])) stop("marker ", mrk1, " not found: marker number out of bounds")
            if(mrk2 > nrow(x$analysis[[1]])) stop("marker ", mrk2, " not found: marker number out of bounds")
        }
        else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
        cat("  Results of the 2-point analysis for markers:", colnames(x$analysis[[1]])[mrk1],
            "and", colnames(x$analysis[[1]])[mrk2], "\n")
        ## results found
        cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
            x$max.rf, "\n\n")
        ## do not print anything for the same marker
        if(mrk1 == mrk2) stop("mrk1 and mrk2 are the same")
        ## results of the two-point analysis
        output<-NULL
        if (mrk1 > mrk2){
            for(i in 1:4)
                output<-rbind(output,c(x$analysis[[i]][mrk1,mrk2],
                                       x$analysis[[i]][mrk2,mrk1]))
            dimnames(output)<-list(c("CC", "CR", "RC", "RR"), c("rf","LOD"))
            print(output)
        }
        else
        {
            for(i in 1:4)
                output<-rbind(output, c(x$analysis[[i]][mrk2,mrk1],
                                        x$analysis[[i]][mrk1,mrk2]))
            dimnames(output)<-list(c("CC", "CR", "RC", "RR"), c("rf","LOD"))
            print(output)
        }
    }
  }
}

##get twopt information for a given pair of markers
get_twopt_info<-function(twopt, small, big)
{
    if(any(class(twopt)=="outcross"))
        return(t(sapply(twopt$analysis, function(x,i,j) c(x[j,i], x[i,j]), i=small, j=big)))
    else
        return(matrix(c(twopt$analysis[big,small], twopt$analysis[small,big]), ncol=2))
}




## end of file


