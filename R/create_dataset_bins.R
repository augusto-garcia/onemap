#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: create_dataset_bins.R                                         #
# Contains: select_data_bins add_redundants                           #
#                                                                     #
# Written by Marcelo Mollinari with minor changes by Cristiane        #
# Taniguti                                                            #
# copyright (c) 2015, Marcelo Mollinari                               #
#                                                                     #
# First version: 09/2015                                              #
# Last update: 08/2021                                                #
# License: GNU General Public License version 3                       #
#                                                                     #
#######################################################################

#' New dataset based on bins
#'
#' Creates a new dataset based on \code{onemap_bin} object
#'
#' Given a \code{onemap_bin} object,
#' creates a new data set where the redundant markers are
#' collapsed into bins and represented by the marker with the lower
#' amount of missing data among those on the bin.
#'
#' @aliases create_data_bins
#' @param input.obj an object of class \code{onemap}.
#' @param bins an object of class \code{onemap_bin}.
#'
#' @return an object of class \code{onemap}.
#' @author Marcelo Mollinari, \email{mmollina@@usp.br}
#' @seealso \code{\link[onemap]{find_bins}}
#' @keywords bins dimension reduction
#' @examples
##'   data("onemap_example_f2")
##'   (bins<-find_bins(onemap_example_f2, exact=FALSE))
##'   onemap_bins <- create_data_bins(onemap_example_f2, bins)
#'@export
create_data_bins <- function(input.obj, bins)
{
  ## checking for correct object
  if(!is(input.obj,"onemap"))
    stop(deparse(substitute(input.obj))," is not an object of class 'onemap'")

  if(!is(bins, "onemap_bin"))
    stop(deparse(substitute(bins))," is not an object of class 'onemap_bin'")

  if (input.obj$n.mar<2) stop("there must be at least two markers to proceed with analysis")

  nm<-names(input.obj)
  dat.temp<-structure(vector(mode="list", length(nm)), class=class(input.obj))
  names(dat.temp)<-nm
  wrk<-match(names(bins$bins), colnames(input.obj$geno))
  dat.temp$geno<-input.obj$geno[,wrk]
  dat.temp$n.ind<-nrow(dat.temp$geno)
  dat.temp$n.mar<-ncol(dat.temp$geno)
  dat.temp$segr.type<-input.obj$segr.type[wrk]
  dat.temp$segr.type.num<-input.obj$segr.type.num[wrk]
  #dat.temp$phase<-input.obj$phase[wrk]
  dat.temp$n.phe<-input.obj$n.phe
  dat.temp$pheno<-input.obj$pheno
  dat.temp$CHROM <- input.obj$CHROM[wrk]
  dat.temp$POS <- input.obj$POS[wrk]
  dat.temp$error <- input.obj$error[wrk + rep(c(0:(input.obj$n.ind-1))*input.obj$n.mar, each=length(wrk)),]
 return(dat.temp)
}

#' Add the redundant markers removed by create_data_bins function
#' 
#' @param sequence object of class \code{sequence}
#' @param onemap.obj object of class \code{onemap.obj} before redundant markers were removed
#' @param bins object of class \code{onemap_bin}
#' 
#' @export
add_redundants <- function(sequence, onemap.obj, bins){
  
  if(!is(sequence, c("sequence"))) stop("Input object must be of class sequence")
  if(!is(onemap.obj, c("onemap"))) stop("Input object must be of class onemap")
  if(!is(bins, c("onemap_bin"))) stop("Input object must be of class onemap_bin")
  
  idx <- match(colnames(sequence$data.name$geno)[sequence$seq.num], names(bins[[1]]))
  sizes <- sapply(bins[[1]][idx], function(x) dim(x)[1])
  
  mks <- sapply(bins[[1]][idx], rownames)
  mks <- do.call(c, mks)
  mks.num <- match(mks, colnames(onemap.obj$geno))
  
  new.seq.rf <- as.list(cumsum(c(0,sequence$seq.rf)))
  
  new.phases <- as.list(c(NA, sequence$seq.phases))
  
  for(i in 1:length(new.seq.rf)){
    new.seq.rf[[i]] <- rep(new.seq.rf[[i]], each = sizes[i])
    new.phases[[i]] <- rep(new.phases[[i]], each = sizes[i])
  }
  
  new.seq.rf <- do.call(c, new.seq.rf)
  new.phases <- do.call(c, new.phases)[-1]
  new.seq.rf <- diff(new.seq.rf)
  new_sequence <- sequence
  new_sequence$seq.num <- mks.num
  new_sequence$seq.phases <- new.phases
  new_sequence$seq.rf <- new.seq.rf
  new_sequence$data.name <- onemap.obj
  new_sequence$probs <- "with redundants"
  return(new_sequence)  
}


