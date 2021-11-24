#########################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: rf_2pts.R                                                     ##
## Contains: rf_2pts, print.rf_2pts                                    ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  ##
## with minor changes by Cristiane Taniguti
## copyright (c) 2007-15, Gabriel R A Margarido and Marcelo Mollinari  ##
##                                                                     ##
## First version: 11/07/2007                                           ##
## Last update: 07/06/2017                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#########################################################################

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
##'@importFrom methods is
##'
##' @aliases rf_2pts
##' @param input.obj an object of class \code{onemap}.
##' @param LOD minimum LOD Score to declare linkage (defaults to \code{3}).
##' @param max.rf maximum recombination fraction to declare linkage (defaults
##' to \code{0.50}).
##' @param verbose logical. If \code{TRUE}, current progress is shown; if
##' \code{FALSE}, no output is produced.
##' @param rm_mks logical. If \code{TRUE} the algorithm will remove the markers for which it found numerical 
##' problems to calculates the recombination fraction. The numerical problems can happens because of excess of 
##' missing data or segregation deviation.
##' 
##' @return An object of class \code{rf_2pts}, which is a list containing the
##' following components:  \item{n.mar}{total number of markers.} \item{LOD}{minimum LOD Score to declare
##' linkage.} \item{max.rf}{maximum recombination fraction to declare linkage.}
##' \item{input}{the name of the input file.} \item{analysis}{an array with the
##' complete results of the two-point analysis for each pair of markers.}
##' 
##' @note The thresholds used for \code{LOD} and \code{max.rf} will be used in
##' subsequent analyses, but can be overriden.
##' @author Gabriel R A Margarido \email{gramarga@@gmail.com} and Marcelo Mollinari \email{mmollina@@usp.br}
##' @references Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002)
##' Simultaneous maximum likelihood estimation of linkage and linkage phases in
##' outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords utilities
##' @examples
##' 
##'   data(onemap_example_out)
##'
##'   twopts <- rf_2pts(onemap_example_out,LOD=3,max.rf=0.5) # perform two-point analyses
##'   twopts
##'
##'   print(twopts,c("M1","M2")) # detailed results for markers 1 and 2
##'   
##'   
##'@export
rf_2pts <- function(input.obj, LOD=3, max.rf=0.50, verbose = TRUE, rm_mks = FALSE) {
  ## checking for correct object
  if(!is(input.obj, c("onemap", "outcross", "f2", "backcross", "riself", "risib")))
    stop(deparse(substitute(input.obj))," is not an object of class 'onemap'.")
  if (input.obj$n.mar<2) stop("there must be at least two markers to proceed with analysis")
  ## creating variables (result storage and progress output)
  if(is(input.obj,"outcross"))
    r<- est_rf_out(geno = input.obj$geno, seg_type = input.obj$segr.type.num, nind = input.obj$n.ind, verbose = verbose)
  else if(is(input.obj, "f2"))
    r<-est_rf_out(geno = input.obj$geno, seg_type = input.obj$segr.type.num, nind = input.obj$n.ind, verbose = verbose)
  else if(is(input.obj,"backcross"))
    r<- est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=0, verbose = verbose)
  else if(is(input.obj,"riself"))
    r<-est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=1, verbose = verbose)
  else if(is(input.obj, "risib"))
    r<-est_rf_bc(geno = input.obj$geno, nind = input.obj$n.ind, type=2, verbose = verbose)
  
  # The recombination matrix should not have NA or NaN
  if(anyNA(r, recursive = T)) {
    if(rm_mks){
      if(is.list(r)){
        mks_rm <- unlist(sapply(r, function(x) rownames(which(is.na(x), arr.ind = T))))
      } else {
        which(apply(r, 2, function(x) all(is.na(x))))
        mks_rm <- unique(rownames(which(is.na(r),arr.ind = T)))
      }
      
      if(length(mks_rm) == input.obj$n.mar){
        stop("Check if all markers have at least one genotype information or if they have segregation pattern deviation. 
            We suggest filter_missing function to avoid excess of missing data and test_segregation.")
      }
      
      message("Recombination fraction for ", length(mks_rm), " markers could not be estimated. They were removed from analysis. Check if these markers have at least one genotype information or if they have segregation pattern deviation. We suggest filter_missing function to avoid excess of missing data and test_segregation.")
      
      mks_rm_num <- which(colnames(input.obj$geno) %in% mks_rm)
      mks_num <- which(!colnames(input.obj$geno) %in% mks_rm)
      input.obj$geno <- input.obj$geno[,-mks_rm_num]
      input.obj$segr.type.num <- input.obj$segr.type.num[-mks_rm_num]
      input.obj$segr.type <- input.obj$segr.type[-mks_rm_num]
      input.obj$CHROM <- input.obj$CHROM[-mks_rm_num]
      input.obj$POS <-  input.obj$POS[-mks_rm_num]
      input.obj$n.mar <- length(input.obj$segr.type)
      input.obj$error <- input.obj$error[mks_num + rep(c(0:(input.obj$n.ind-1))*input.obj$n.mar, each=length(mks_num)),]
      
      twopts <- rf_2pts(input.obj, LOD=LOD, max.rf=max.rf, verbose = verbose)
      return(twopts)
    } else 
      message("We could not estimate all recombination fraction. Check if these markers have at least one genotype information or if they have segregation pattern deviation. We suggest filter_missing function to avoid excess of missing data, test_segregation and rm_mks argument.")
  }
  
  structure(list(data.name= input.obj, n.mar=input.obj$n.mar, LOD=LOD, max.rf=max.rf,
                 input=input.obj$input, CHROM = input.obj$CHROM, POS= input.obj$POS, analysis=r),
            class = c("rf_2pts", class(input.obj)[2]))
}

##' Print method for object class 'rf_2pts'
##'
##' It shows the linkage groups as well as the unlinked markers.
##'
##' @aliases print.rf_2pts
##'
##' @param x an object of class \code{rf_2pts}.
##'
##' @param mrk a vector containing a pair of markers, so detailed
##'     results of the two-point analysis will be printed for them.
##'     Can be numeric or character strings indicating the
##'     numbers/names corresponding to any markers in the input file.
##'
##' @param ... further arguments, passed to other methods. Currently ignored.
##'
##' @return \code{NULL}
##' @keywords internal
##' @method print rf_2pts
##' @export
print.rf_2pts <- function(x, mrk=NULL,...) {
  ## checking for correct object
  if(!is(x, "rf_2pts"))
    stop(deparse(substitute(x))," is not an object of class 'rf_2pts'")
  if (any(is.null(mrk))) {
    ## printing a brief summary
    cat("  This is an object of class 'rf_2pts'\n")
    cat("\n  Criteria: LOD =", x$LOD, ", Maximum recombination fraction =",
        x$max.rf, "\n")
    cat("\n  This object is too complex to print\n")
    cat("  Type 'print(object, c(mrk1=marker, mrk2=marker))' to see\n")
    cat("    the analysis for two markers\n")
    cat("    mrk1 and mrk2 can be the names or numbers of both markers\n")
  }
  else {
    ## printing detailed results for two markers
    ## checking if markers exist and converting character to numeric
    if(length(mrk)!=2)
      stop(deparse(substitute(mrk))," must be a pair of markers")
    if(is(x, "backcross") || is(x, "risib") || is(x, "riself"))
    {
      if (is.character(mrk[1]) && is.character(mrk[2])) {
        mrk1name<-mrk[1]
        mrk2name<-mrk[2]
        mrk[1]<-match(mrk[1], colnames(x$analysis))
        mrk[2]<-match(mrk[2], colnames(x$analysis))
        mrk<-as.numeric(mrk)
        if (is.na(mrk[1])) stop("marker ", mrk1name, " not found")
        if (is.na(mrk[2])) stop("marker ", mrk2name, " not found")
      }
      else if(is.numeric(mrk[1]) && is.numeric(mrk[2]))
      {
        if(mrk[1] > nrow(x$analysis)) stop("marker ", mrk[1], " not found: marker number out of bounds")
        if(mrk[2] > nrow(x$analysis)) stop("marker ", mrk[2], " not found: marker number out of bounds")
      }
      else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
      cat("  Results of the 2-point analysis for markers:", colnames(x$analysis)[mrk[1]],
          "and", colnames(x$analysis)[mrk[2]], "\n")
      
      ## results found
      cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
          x$max.rf, "\n\n")
      ## do not print anything for the same marker
      if(mrk[1] == mrk[2]) stop("mrk1 and mrk2 are the same")
      ## results of the two-point analysis
      if (mrk[1] > mrk[2]){
        r<-x$analysis[mrk[1],mrk[2]]
        LOD<-x$analysis[mrk[2],mrk[1]]
      }
      else
      {
        r<-x$analysis[mrk[2],mrk[1]]
        LOD<-x$analysis[mrk[1],mrk[2]]
      }
      output<-c(r, LOD)
      names(output)<-c("rf","LOD")
      print(output)
    }
    else if(is(x, "outcross") || is(x,"f2"))
    {
      if (is.character(mrk[1]) && is.character(mrk[2])) {
        mrk1name<-mrk[1]
        mrk2name<-mrk[2]
        mrk[1]<-match(mrk[1], colnames(x$analysis[[1]]))
        mrk[2]<-match(mrk[2], colnames(x$analysis[[1]]))
        mrk<-as.numeric(mrk)
        if (is.na(mrk[1])) stop("marker ", mrk1name, " not found")
        if (is.na(mrk[2])) stop("marker ", mrk2name, " not found")
      }
      else if(is.numeric(mrk[1]) && is.numeric(mrk[2]))
      {
        if(mrk[1] > nrow(x$analysis[[1]])) stop("marker ", mrk[1], " not found: marker number out of bounds")
        if(mrk[2] > nrow(x$analysis[[1]])) stop("marker ", mrk[2], " not found: marker number out of bounds")
      }
      else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
      cat("  Results of the 2-point analysis for markers:", colnames(x$analysis[[1]])[mrk[1]],
          "and", colnames(x$analysis[[1]])[mrk[2]], "\n")
      ## results found
      cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
          x$max.rf, "\n\n")
      ## do not print anything for the same marker
      if(mrk[1] == mrk[2]) stop("mrk1 and mrk2 are the same")
      ## results of the two-point analysis
      output<-NULL
      if (mrk[1] > mrk[2]){
        for(i in 1:4)
          output<-rbind(output,c(x$analysis[[i]][mrk[1],mrk[2]],
                                 x$analysis[[i]][mrk[2],mrk[1]]))
      }
      else
      {
        for(i in 1:4)
          output<-rbind(output, c(x$analysis[[i]][mrk[2],mrk[1]],
                                  x$analysis[[i]][mrk[1],mrk[2]]))
      }
      dimnames(output)<-list(c("CC", "CR", "RC", "RR"), c("rf","LOD"))
      print(output)
    }
  }
}

##get twopt information for a given pair of markers
get_twopt_info<-function(twopt, small, big)
{
  if(is(twopt, "outcross"))
    return(t(sapply(twopt$analysis, function(x,i,j) c(x[j,i], x[i,j]), i=small, j=big)))
  else
    return(matrix(c(twopt$analysis[big,small], twopt$analysis[small,big]), ncol=2))
}

## end of file


