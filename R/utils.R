#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: utils.R                                                       #
# Contains: acum seq_by_type map_avoid_unlinked split_2pts            #
# map_save_ram remove_inds sort_by_pos empty_onemap_obj               #
# try_seq_by_seq add_redundants rm_dupli_mks                          #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Cristiane Taniguti #
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
                               phase_cores = 1, 
                               tol = 1e-05){
  #TODO: error checks...
  map_df <- map_save_ram(input.seq, rm_unlinked = T, 
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
  input.seq_ram <- input.seq
  if(length(input.seq$seq.num) < input.seq.tot$data.name$n.mar){
    split.twopts <- split_2pts(twopts.obj = input.seq$twopt, mks = input.seq$seq.num) 
    input.seq_ram <- make_seq(split.twopts, "all")
  }
  if(phase_cores == 1){
    return.map <- map(input.seq_ram, tol = tol, 
                      verbose = verbose, 
                      rm_unlinked = rm_unlinked, 
                      phase_cores = phase_cores)
  } else {
    if(is.null(size) | is.null(overlap)){
      stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you should also define `size` and `overlap` arguments. See ?map_avoid_unlinked and ?pick_batch_sizes")
    } else {
      return.map <- map_overlapping_batches(input.seq = input.seq_ram,
                                            size = size, overlap = overlap, 
                                            phase_cores = phase_cores, 
                                            tol=tol, rm_unlinked = rm_unlinked)
    }
  }
  if(length(input.seq_ram$seq.num) < input.seq.tot$data.name$n.mar){
    if(!is(return.map, "integer")){ # When rm_unlinked == F
      return.map$seq.num <- input.seq.tot$seq.num
      return.map$data.name <- input.seq.tot$data.name
      return.map$twopt <- input.seq.tot$twopt
    } else {
      remain <- colnames(input.seq_ram$data.name$geno)[return.map]
      old <- colnames(input.seq.tot$data.name$geno)[input.seq.tot$seq.num]
      return.map <- input.seq.tot$seq.num[old %in% remain] 
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

#' Sort markers in onemap object by their position in reference genome
#' 
#' @param onemap.obj object of class onemap
#' 
#' @export
sort_by_pos <- function(onemap.obj){
  
  idx <- order(onemap.obj$CHROM, onemap.obj$POS)
  
  new.obj <- onemap.obj
  new.obj$geno <- onemap.obj$geno[,idx]
  new.obj$segr.type <- onemap.obj$segr.type[idx]
  new.obj$segr.type.num <- onemap.obj$segr.type.num[idx]
  new.obj$CHROM <- onemap.obj$CHROM[idx]
  new.obj$POS <- onemap.obj$POS[idx]
  new.obj$error <- onemap.obj$error[idx + rep(c(0:(onemap.obj$n.ind-1))*onemap.obj$n.mar, each=length(idx)),]
  return(new.obj)
}

# Produce empty object to avoid code break
empty_onemap_obj <- function(vcf, P1, P2, cross){
  legacy_crosses <- setNames(c("outcross", "f2", "backcross", "riself", "risib"), 
                             c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"))
  
  geno <- matrix(0, ncol = 0, nrow = length(colnames(vcf@gt)[-c(1, P1, P2)]))
  
  rownames(geno) <- colnames(vcf@gt)[-c(1, P1, P2)]    
  onemap.obj <- structure(list(geno= geno,
                               n.ind = dim(geno)[2],
                               n.mar = 0,
                               segr.type = logical(),
                               segr.type.num = as.numeric(),
                               n.phe = 0,
                               pheno = NULL,
                               CHROM = logical(),
                               POS = logical(),
                               input = "vcfR.object"),
                          class=c("onemap",legacy_crosses[cross]))
  return(onemap.obj)
}


try_seq_by_seq <- function(sequence, markers, cM.thr= 10, lod.thr=-10){
  
  seq_now <- sequence
  for(i in 1:length(markers)){
    try_edit <- try_seq(seq_now, markers[i])  
    pos <- which(try_edit$LOD == 0)[1]
    new_map <- make_seq(try_edit, pos)
    size_new <- cumsum(kosambi(new_map$seq.rf))[length(new_map$seq.rf)]
    size_old <- cumsum(kosambi(seq_now$seq.rf))[length(seq_now$seq.rf)]
    lod_new <- new_map$seq.like
    lod_old <- seq_now$seq.like
    diff_size <- size_new - size_old
    diff_lod <- lod_new - lod_old
    if(diff_size < cM.thr & diff_lod > lod.thr){
      seq_now <- new_map
      cat("Marker", markers[i], "was included \n")
    } 
  }
  return(seq_now)
}

add_redundants <- function(sequence, onemap.obj, bins){
  
  idx <- match(colnames(sequence$data.name$geno)[sequence$seq.num], names(bins[[1]]))
  
  sizes <- sapply(bins[[1]][idx], function(x) dim(x)[1])
  
  mks <- sapply(bins[[1]][idx], rownames)
  mks <- do.call(c, mks)
  mks.num <- match(mks, colnames(onemap.obj$geno))
  
  new.seq.rf <- as.list(cumsum(c(0,sequence$seq.rf)))
  
  for(i in 1:length(new.seq.rf)){
    new.seq.rf[[i]] <- rep(new.seq.rf[[i]], each = sizes[i])
  }
  
  new.seq.rf <- do.call(c, new.seq.rf)
  new.seq.rf <- diff(new.seq.rf)
  new_sequence <- sequence
  new_sequence$seq.num <- mks.num
  new_sequence$seq.rf <- new.seq.rf
  new_sequence$data.name <- onemap.obj
  new_sequence$probs <- "with redundants"
  return(new_sequence)  
}

# Removing duplicated markers, kept the one with less missing data
rm_dupli_mks <- function(onemap.obj){
  MKS <- colnames(onemap.obj$geno)
  GT_matrix <- t(onemap.obj$geno)
  n.mk <- length(MKS)
  dupli <- MKS[duplicated(MKS)]
  if(length(dupli)>0){
    n.rm.mks <- length(dupli)
    dupli <- unique(dupli)
    warning(paste("There are duplicated markers IDs:", paste(MKS[duplicated(MKS)], collapse = " "), "\nOnly the one with less missing data was kept."))
    for(w in 1:length(dupli)){
      temp_GT <- GT_matrix[MKS==dupli[w],]
      mis_count <- apply(temp_GT, 1, function(x) sum(x==0))
      discard <- temp_GT[-which.min(mis_count),]
      if(is(discard, "matrix")){
        for(j in 1:dim(discard)[1]){
          idx <- which(apply(GT_matrix, 1, function(x) all(x == discard[j,])))
          idx <- idx[MKS[idx] == dupli[w]][1]
          GT_matrix <- GT_matrix[-idx,]
          mk.type <- mk.type[-idx]
          mk.type.num <- mk.type.num[-idx] 
          onemap.obj$CHROM <- onemap.obj$CHROM[-idx]
          onemap.obj$POS <- onemap.obj$POS[-idx]
          MKS <- MKS[-idx]
        }
      } else {
        idx <- which(apply(GT_matrix, 1, function(x) all(x == discard)))
        idx <- idx[MKS[idx] == dupli[w]][1]
        GT_matrix <- GT_matrix[-idx,]
        onemap.obj$segr.type <- onemap.obj$segr.type[-idx]
        onemap.obj$segr.type.num <- onemap.obj$segr.type.num[-idx] 
        onemap.obj$CHROM <- onemap.obj$CHROM[-idx]
        onemap.obj$POS <- onemap.obj$POS[-idx]
        MKS <- MKS[-idx]
      }
    }
    onemap.obj$n.mar <- n.mk - n.rm.mks
    onemap.obj$geno <- t(GT_matrix)
  } 
  return(onemap.obj)
}

#' Onemap object sanity check 
#' 
#' Based on MAPpoly check_data_sanity function by Marcelo Mollinari
#' 
#' @param x an object of class \code{onemap}
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @examples 
#' 
#' data(onemap_example_bc)
#' check_data(onemap_example_bc)
#' 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@usp.br}
#' 
#' @export
check_data <- function(x){
  test <- logical(24L)
  names(test) <- 1:24
  
  test[1] <- any(is.na(x$geno))    
  test[2] <- any(is.na(x$error))
  test[3] <- !all(dim(x$geno) == c(x$n.ind, x$n.mar))
  test[4] <- !dim(x$error)[1] == prod(dim(x$geno))
  test[5] <- if(!is.null(x$CHROM)) length(x$CHROM) != x$n.mar else FALSE
  test[6] <- if(!is.null(x$POS)) length(x$POS) != x$n.mar else FALSE
  test[7] <- if(is(x, "f2")) {
    !all(unique(x$segr.type) %in% c("A.H.B", "D.B", "C.A"))
  } else if(is(x, "outcross")){
    !all(unique(x$segr.type) %in% c("A.1", "A.2", "A.3", "A.4", "B1.5", 
                                    "B2.6", "B3.7", "C.8", "D1.9", "D1.10", 
                                    "D1.11", "D1.12", "D1.13", "D2.14", 
                                    "D2.15", "D2.16", "D2.17", "D2.18"))
  } else if(is(x, "backcross")){
    !all(unique(x$segr.type) %in% c("A.H"))
  } else if(is(x, "riself") | is(x, "risib")){
    !all(unique(x$segr.type) %in% c("A.B"))
  }
  
  test[8] <- if(is(x, "f2")) {
    !all(unique(x$segr.type.num) %in% c(4,6,7))
  } else if(is(x, "outcross")){
    !all(unique(x$segr.type.num) %in% 1:7)
  } else if(is(x, "backcross")){
    !all(unique(x$segr.type.num) %in% 8)
  } else if(is(x, "riself") | is(x, "risib")){
    !all(unique(x$segr.type.num) %in% 9)
  }
  
  if(any(test))
    return(test)
  else 
    return(0)
}

#' Twopts object sanity check 
#' 
#' Based on MAPpoly check_data_sanity function by Marcelo Mollinari
#' 
#' @param x an object of class \code{onemap}
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @examples 
#' 
#' data(onemap_example_bc)
#' twopts <- rf_2pts(onemap_example_bc)
#' check_twopts(twopts)
#' 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@usp.br}
#' 
#' @export
check_twopts <- function(x){
  test <- logical(4L)
  names(test) <- 1:4
  
  test[1] <- if(check_data(x$data.name) == 0) FALSE else TRUE
  if(is(x$data.name, "outcross") | is(x$data.name, "f2")){
    test[2] <- !is(x$analysis, "list")
    test[3] <- all(dim(x$analysis[[1]]) != rep(x$data.name$n.mar,2))
    test[4] <- any(sapply(x$analysis, function(x) any(is.na(x))))
  } else {
    test[2] <- !is(x$analysis, "matrix")
    test[3] <- all(dim(x$analysis) != rep(x$data.name$n.mar,2))
    test[4] <- any(is.na(x$analysis))
  }
  
  if(any(test))
    return(test)
  else 
    return(0)
}

