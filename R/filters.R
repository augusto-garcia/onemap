#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: filters.R                                                     ##
## Contains: filter_missing filter_by_prob                             ##
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
  if(!is(onemap.obj,"onemap")){
    stop("onemap.obj should be of class onemap\n")
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
  cat("Number of markers removed from the onemap object: ", length(which(perc.mis > threshold)), "\n")
  return(new.onemap.obj)
}


##' Function filter genotypes by genotype probability
##'
##' @param onemap.obj an object of class \code{onemap}.
##' @param threshold a numeric from 0 to 1 to define the threshold for 
##' the probability of the called genotype (highest probability)
##' @return Returns an object of class \code{onemap}
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
##' @examples
##'
##'   data(onemap_example_out)
##'   filt_obj <- filter_prob(onemap_example_out, threshold=0.8)
##' 
##' 
##' @importFrom reshape2 melt dcast
##' 
##' @export
filter_prob <- function(onemap.obj=NULL, threshold= 0.8){
  idx <- apply(onemap.obj$error, 1, which.max)
  rm <- which(onemap.obj$error[cbind(seq_along(idx), idx)] < threshold)
  onemap.obj$error[rm,] <- 1
  cat(paste(length(rm), "genotypes were converted to missing data."))
  
  geno_melt <- melt(onemap.obj$geno)
  geno_melt[rm,3] <- 0
  geno <- dcast(geno_melt, Var1 ~ Var2)
  rownames(geno) <- geno$Var1
  geno <- geno[,-1]
  
  onemap.obj$geno <- as.matrix(geno)
  
  return(onemap.obj)
}
