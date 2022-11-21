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
##' individual.} \item{error}{matrix containing HMM emission probabilities}
#' 
#' 
#' @author Marcelo Mollinari, \email{mmollina@@usp.br}
#' @seealso \code{\link[onemap]{find_bins}}
#' @keywords bins dimension reduction
#' @examples
#' 
##'   data("onemap_example_f2")
##'   (bins<-find_bins(onemap_example_f2, exact=FALSE))
##'   onemap_bins <- create_data_bins(onemap_example_f2, bins)
#' 
#'@export
create_data_bins <- function(input.obj, bins)
{
  ## checking for correct object
  if(!inherits(input.obj,"onemap"))
    stop(deparse(substitute(input.obj))," is not an object of class 'onemap'")

  if(!inherits(bins, "onemap_bin"))
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
##' @return New sequence object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{object of class \code{onemap} with the raw
##' data.} \item{twopt}{object of class \code{rf_2pts} with the
##' 2-point analyses.} 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#' 
#' @seealso \code{\link[onemap]{find_bins}}
#' 
#' @keywords redundants bins
#' 
#' @export
add_redundants <- function(sequence, onemap.obj, bins){
  
  if(!inherits(sequence, c("sequence"))) stop("Input object must be of class sequence")
  if(!inherits(onemap.obj, c("onemap"))) stop("Input object must be of class onemap")
  if(!inherits(bins, c("onemap_bin"))) stop("Input object must be of class onemap_bin")
  
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


