#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: utils.R                                                       #
# Contains: acum                                                      #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

acum <- function(w) {
  if (w<0) stop("'w' should be equal to or higher than zero")
  
  # the famous gaussian sum from 1 to w
  w*(w+1)/2
}


#' Extract marker number by name
#' 
#' @param sequence object of class or sequence
#' @param mk_type vector of character with marker type to be selected
#' 
#' @return new sequence object with selected marker type
#' @export
#' 
seq_by_type <- function(sequence, mk_type){
  if(!is(sequence, c("sequence"))) stop("Input object must be of class sequence")
  if(length(mk_type) > 1) pat <- paste0(mk_type, collapse = "|") else pat <- mk_type
  type <- sequence$seq.num[grep(pat, sequence$data.name$segr.type[sequence$seq.num])]
  new.seq <- make_seq(sequence$twopt, type)
  return(new.seq)
}

#' Repeat HMM if map find unlinked maker
#'
#' @param input.seq object of class sequence
#' @param size The center size around which an optimum is to be searched
#' @param overlap The desired overlap between batches
#' @param phase_cores The number of parallel processes to use when estimating
#' the phase of a marker. (Should be no more than 4)
#' @param tol tolerance for the C routine, i.e., the value used to evaluate
#' convergence.
#' 
#' @export
map_avoid_unlinked <- function(input.seq, 
                               size = NULL, 
                               overlap = NULL,
                               round = 5,
                               phase_cores = 1, 
                               tol = 1e-05,
                               max.gap=F){
  #TODO: error checks...
  if(phase_cores == 1){
    map_df <- map(input.seq, rm_unlinked = T)
  } else {
    if(is.null(size) | is.null(overlap)){
      stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you should also define `size` and `overlap` arguments. See ?map_avoid_unlinked and ?pick_batch_sizes")
    } else {
    map_df <- map_overlapping_batches(input.seq = input.seq,
                                      size = size, overlap = overlap,
                                      phase_cores = phase_cores, 
                                      tol=tol, rm_unlinked = T, max.gap = max.gap)
    }
  }
  while(is(map_df, "integer")){
    seq_true <- make_seq(input.seq$twopt, map_df)
    if(is.null(size) & is.null(overlap) & phase_cores == 1){
      map_df <- map(input.seq = seq_true, rm_unlinked = T, tol=tol)
    }else{
      map_df <- map_overlapping_batches(input.seq = seq_true,
                                        size = size, overlap = overlap, 
                                        phase_cores = phase_cores, 
                                        tol=tol, rm_unlinked = T)
    }
  }
  return(map_df)
}


