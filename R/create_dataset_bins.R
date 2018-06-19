#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: create_dataset_bins.R                                         #
# Contains: select_data_bins                                          #
#                                                                     #
# Written by Marcelo Mollinari with minor changes by Cristiane        #
# Taniguti                                                            #
# copyright (c) 2015, Marcelo Mollinari                               #
#                                                                     #
# First version: 09/2015                                              #
# Last update: 07/05/2017                                             #
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
#'  \dontrun{
#'   load(url("https://github.com/mmollina/data/raw/master/fake_big_data_f2.RData"))
#'   fake.big.data.f2
#'   (bins <- find_bins(fake.big.data.f2, exact=FALSE))
#'   (new.data <- create_data_bins(fake.big.data.f2, bins))}
#'@export
create_data_bins <- function(input.obj, bins)
{
  ## checking for correct object
  if(class(input.obj)[1] != "onemap")
    stop(deparse(substitute(input.obj))," is not an object of class 'onemap'")

  if(is.na(match("onemap_bin", class(bins))))
    stop(deparse(substitute(bins))," is not an object of class 'onemap_bin'")

  if (input.obj$n.mar<2) stop("there must be at least two markers to proceed with analysis")

  nm<-names(input.obj)
  dat.temp<-structure(vector(mode="list", length(nm)), class=class(input.obj))
  names(dat.temp)<-nm
  wrk<-match(names(bins$bins), colnames(input.obj$geno))
  dat.temp$geno<-input.obj$geno[,wrk]
  if(class(input.obj)[2] != "outcross"){
    dat.temp$geno.mmk<-list(geno=dat.temp$geno, type=gsub("\\..*","",class(input.obj)[2]))
    dat.temp$geno.mmk$geno[dat.temp$geno.mmk$geno==0]<-NA
  }
  dat.temp$n.ind<-nrow(dat.temp$geno)
  dat.temp$n.mar<-ncol(dat.temp$geno)
  dat.temp$segr.type<-input.obj$segr.type[wrk]
  dat.temp$segr.type.num<-input.obj$segr.type.num[wrk]
  #dat.temp$phase<-input.obj$phase[wrk]
  dat.temp$n.phe<-input.obj$n.phe
  dat.temp$pheno<-input.obj$pheno
  dat.temp$CHROM <- input.obj$CHROM[wrk]
  dat.temp$POS <- input.obj$POS[wrk]
 return(dat.temp)
}

## end of file


