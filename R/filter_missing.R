#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: filter_missing.R                                              ##
## Contains: filter_missing                                            ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
## copyright (c) 2007-9, Cristiane Taniguti                            ##
##                                                                     ##
## First version: 22/11/2019                                           ## 
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

## Function filter markers by missing data

##' Filter markers according with a missing data threshold
##'
##' @param onemap.obj an object of class \code{onemap}.
##' @param threshold a numeric from 0 to 1 to define the threshold of missing data allowed
##' @return Returns an object of class \code{onemap}
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
##' @examples
##'
##'   data(onemap_example_out)
##'   filt_obj <- filter_missing(onemap_example_out, threshold=0.25)
##'   
##'@export
filter_missing <- function(onemap.obj=NULL, threshold= 0.25){
  if(class(onemap.obj)[1]!="onemap"){
    stop("onemap.obj should be of class onemap")
  }
  perc.mis <- apply(onemap.obj$geno, 2, function(x) sum(x == 0)/length(x))
  idx <- which(!perc.mis > threshold)
  
  new.onemap.obj <- onemap.obj
  new.onemap.obj$geno <- onemap.obj$geno[,idx]
  new.onemap.obj$n.mar <- length(idx)
  new.onemap.obj$segr.type <- onemap.obj$segr.type[idx]
  new.onemap.obj$segr.type.num <- onemap.obj$segr.type.num[idx]
  new.onemap.obj$CHROM <- onemap.obj$CHROM[idx]
  new.onemap.obj$POS <- onemap.obj$POS[idx]
  new.onemap.obj$error <- onemap.obj$error[idx + rep(c(0:(onemap.obj$n.ind-1))*onemap.obj$n.mar, each=length(idx)),]
  cat("Number of markers removed from the onemap object: ", length(which(perc.mis > threshold)))
  return(new.onemap.obj)
}