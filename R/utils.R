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
  new.seq <- make_seq(sequence$data.name, type)
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
                               phase_cores = 1, 
                               tol = 1e-05){
  #TODO: error checks...
  map_df <- onemap:::map_save_ram(input.seq, rm_unlinked = T, 
                                  size = size, 
                                  overlap = overlap, 
                                  tol=tol, 
                                  phase_cores = phase_cores)
  
  while(is(map_df, "integer")){
    seq_true <- make_seq(input.seq$twopt, map_df)
    map_df <- map_save_ram(input.seq = seq_true, 
                           rm_unlinked = T, 
                           tol=tol, 
                           size = size, 
                           overlap = overlap, 
                           phase_cores = phase_cores)
  }
  return(map_df)
}

# Split 2pts object by mks
split_2pts <- function(twopts.obj, mks){
  split.dat <- split_onemap(onemap.obj = twopts.obj$data.name, mks)
  twopts.obj$data.name <- split.dat
  twopts.obj$n.mar <- length(mks)
  twopts.obj$CHROM <- twopts.obj$CHROM[mks]
  twopts.obj$POS <- twopts.obj$POS[mks]
  if(is(twopts.obj$data.name, c("outcross","f2"))){
    new.twopts <- rep(list(matrix(0,nrow = length(mks), ncol = length(mks))),4)
    for(j in 1:(length(mks)-1)) {
      for(i in (j+1):length(mks)) {
        k<-sort(c(mks[i], mks[j]))
        for(w in 1:4){
          r.temp<-twopts.obj$analysis[[w]][k[2], k[1]]
          new.twopts[[w]][i,j]<-r.temp
          LOD.temp<-twopts.obj$analysis[[w]][k[1], k[2]]
          new.twopts[[w]][j,i]<-LOD.temp
          colnames(new.twopts[[w]]) <- rownames(new.twopts[[w]]) <- colnames(split.dat$geno)
        }
      }
    }
    names(new.twopts) <- c("CC", "CR", "RC", "RR")
  } else {
    new.twopts <- matrix(0, nrow = length(mks), ncol = length(mks))
    for(i in 1:(length(mks)-1)) {
      for(j in (i+1):length(mks)) {
        k<-sort(c(mks[i], mks[j]))
        r.temp<-twopts.obj$analysis[k[2], k[1]]
        new.twopts[i,j]<-r.temp
        LOD.temp<-twopts.obj$analysis[k[1], k[2]]
        new.twopts[j,i]<-LOD.temp
      }
    }
    colnames(new.twopts) <- rownames(new.twopts) <- colnames(split.dat$geno)
  }
  twopts.obj$analysis <- new.twopts
  return(twopts.obj)
}

# perform map with backgroups onemap object and twopts only with sequence markers information
# it save space in ram memory - very useful if dealing with many markers in total dataset
map_save_ram <- function(input.seq,
                         tol=10E-5, 
                         verbose=FALSE, 
                         rm_unlinked=FALSE, 
                         phase_cores = 1, 
                         size = NULL, 
                         overlap = NULL){
  
  input.seq.tot <- input.seq
  if(length(input.seq$seq.num) < input.seq.tot$data.name$n.mar){
    split.twopts <- onemap:::split_2pts(twopts.obj = input.seq$twopt, mks = input.seq$seq.num) 
    input.seq <- make_seq(split.twopts, "all")
  }
  if(phase_cores == 1){
    return.map <- map(input.seq, tol = tol, 
                      verbose = verbose, 
                      rm_unlinked = rm_unlinked, 
                      phase_cores = phase_cores)
  } else {
    if(is.null(size) | is.null(overlap)){
      stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you should also define `size` and `overlap` arguments. See ?map_avoid_unlinked and ?pick_batch_sizes")
    } else {
      return.map <- map_overlapping_batches(input.seq = input.seq,
                                            size = size, overlap = overlap, 
                                            phase_cores = phase_cores, 
                                            tol=tol, rm_unlinked = rm_unlinked)
    }
  }
  if(length(input.seq$seq.num) < input.seq.tot$data.name$n.mar){
    if(!is(return.map, "integer")){ # When rm_unlinked == F
      return.map$seq.num <- input.seq.tot$seq.num
      return.map$data.name <- input.seq.tot$data.name
      return.map$twopt <- input.seq.tot$twopt
    }
  }
  return(return.map)
}


#'Remove inviduals from the onemap object
#'
#'@param onemap.obj object of class onemap
#'@param rm.ind vector of charaters with individuals names
#'
#'@export
remove_inds <- function(onemap.obj, rm.ind){
  if(!is(onemap.obj, "onemap")) stop("Input must to be of onemap class \n")
  if(!(length(which(rownames(onemap.obj$geno) %in% rm.ind)) >0)) stop("We could not find any of these individuals in the dataset \n")
  
  rm.ind <- c("II_3_08", "II_1_37", "I_4_62", "I_4_28", "I_4_21", "I_3_72", "I_3_70")
  new.onemap.obj <- onemap.obj
  new.onemap.obj$geno <- onemap.obj$geno[-which(rownames(onemap.obj$geno) %in% rm.ind),]
  new.onemap.obj$n.ind <- onemap.obj$n.ind - length(rm.ind)
  for(i in 1:length(rm.ind)){
    rm.idx <- grep(paste0("_",rm.ind[i],"$"), rownames(new.onemap.obj$error))
    new.onemap.obj$error <- new.onemap.obj$error[-rm.idx,]
  }
  return(new.onemap.obj)
}
