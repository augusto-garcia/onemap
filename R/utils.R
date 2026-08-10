#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: utils.R                                                       #
# Contains: acum seq_by_type map_avoid_unlinked split_2pts            #
# map_save_ram remove_inds sort_by_pos empty_onemap_obj               #
# try_seq_by_seq add_redundants rm_dupli_mks ord_by_geno              #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Cristiane Taniguti #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#' Perform gaussian sum
#' 
#' @param w vector of numbers
#' 
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
##' @return New sequence object of class \code{sequence} with selected marker type, 
##' which is a list containing the
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
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' @seealso \code{\link[onemap]{make_seq}}
##' 
#' @export
seq_by_type <- function(sequence, mk_type){
  if(!inherits(sequence, c("sequence"))) stop("Input object must be of class sequence")
  if(length(mk_type) > 1) pat <- paste0(mk_type, collapse = "|") else pat <- mk_type
  type <- sequence$seq.num[grep(pat, sequence$data.name$segr.type[sequence$seq.num])]
  new.seq <- make_seq(sequence$twopt, type)
  return(new.seq)
}

#' Split rf_2pts object by markers
#' 
#' @param twopts.obj object of class rf_2pts
#' @param mks markers names (vector of characters) or number (vector of integers) to be removed and added to a new rf_2pts object
#' 
##' @return An object of class \code{rf_2pts} with only the selected markers, which is a list containing the
##' following components:  \item{n.mar}{total number of markers.} \item{LOD}{minimum LOD Score to declare
##' linkage.} \item{max.rf}{maximum recombination fraction to declare linkage.}
#' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' 
#' @export
split_2pts <- function(twopts.obj, mks){
  split.dat <- split_onemap(onemap.obj = twopts.obj$data.name, mks)
  twopts.obj$data.name <- split.dat
  twopts.obj$n.mar <- length(mks)
  twopts.obj$CHROM <- twopts.obj$CHROM[mks]
  twopts.obj$POS <- twopts.obj$POS[mks]
  if(inherits(twopts.obj$data.name, c("outcross","f2"))){
    new.twopts <- lapply(twopts.obj$analysis, function(x){
      temp <- matrix(0,nrow = length(mks), ncol = length(mks))
      k <- matrix(c(rep(mks[1:(length(mks))], each = length(mks)), 
                    rep(mks[1:(length(mks))], length(mks))), ncol = 2)
      k <- k[-which(k[,1] == k[,2]),]
      k <- t(apply(k, 1, sort))
      k <- k[-which(duplicated(k)),]
      LOD.temp<-x[k[,c(1,2)]]
      temp[lower.tri((temp))]<-LOD.temp
      temp <- t(temp) 
      r.temp<-x[k[,c(2,1)]]
      temp[lower.tri(temp)]<-r.temp
      colnames(temp) <- rownames(temp) <- colnames(split.dat$geno)
      return(temp)
    })
    names(new.twopts) <- c("CC", "CR", "RC", "RR")
  } else {
    new.twopts <- matrix(0,nrow = length(mks), ncol = length(mks))
    k <- matrix(c(rep(mks[1:(length(mks))], each = length(mks)), 
                  rep(mks[1:(length(mks))], length(mks))), ncol = 2)
    k <- k[-which(k[,1] == k[,2]),]
    k <- t(apply(k, 1, sort))
    k <- k[-which(duplicated(k)),]
    LOD.temp<- twopts.obj$analysis[k[,c(1,2)]]
    new.twopts[lower.tri((new.twopts))] <- LOD.temp
    new.twopts <- t(new.twopts) 
    r.temp<- twopts.obj$analysis[k[,c(2,1)]]
    new.twopts[lower.tri(new.twopts)] <- r.temp
    colnames(new.twopts) <- rownames(new.twopts) <- colnames(split.dat$geno)
  }
  twopts.obj$analysis <- new.twopts
  return(twopts.obj)
}


#' Remove individuals from the onemap object
#'
#' @param onemap.obj object of class onemap
#' @param rm.ind vector of characters with individuals names
#' @param list.seqs list of objects of class sequence
#'
#' @return An object of class \code{onemap} without the selected individuals
#' if onemap object is used as input, or a list of objects of class \code{sequence}
#' without the selected individuals if a list of sequences objects is use as input
#'
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#'
#'@export
remove_inds <- function(onemap.obj=NULL, rm.ind=NULL, list.seqs = NULL){
  if(!is.null(onemap.obj)){
    if(!inherits(onemap.obj, "onemap")) stop("Input must to be of onemap class \n")
    if(!(length(which(rownames(onemap.obj$geno) %in% rm.ind)) >0)) stop("We could not find any of these individuals in the dataset \n")
    
    new.onemap.obj <- onemap.obj
    new.onemap.obj$geno <- onemap.obj$geno[-which(rownames(onemap.obj$geno) %in% rm.ind),]
    new.onemap.obj$n.ind <- onemap.obj$n.ind - length(rm.ind)
    for(i in 1:length(rm.ind)){
      rm.idx <- grep(paste0("_",rm.ind[i],"$"), rownames(new.onemap.obj$error))
      new.onemap.obj$error <- new.onemap.obj$error[-rm.idx,]
    }
    return(new.onemap.obj)
  } else if(!is.null(list.seqs)){
    new.onemap.obj <- list.seqs[[1]]$data.name
    if(!(length(which(rownames(new.onemap.obj$geno) %in% rm.ind)) >0)) stop("We could not find any of these individuals in the dataset \n")
    
    new.onemap.obj$geno <- new.onemap.obj$geno[-which(rownames(new.onemap.obj$geno) %in% rm.ind),]
    new.onemap.obj$n.ind <- new.onemap.obj$n.ind - length(rm.ind)
    for(i in 1:length(rm.ind)){
      rm.idx <- grep(paste0("_",rm.ind[i],"$"), rownames(new.onemap.obj$error))
      new.onemap.obj$error <- new.onemap.obj$error[-rm.idx,]
    }
    
    new.list.seqs <- list.seqs
    for(i in 1:length(list.seqs)){
      new.list.seqs[[i]]$data.name <- new.onemap.obj
    }
    return(new.list.seqs)
  } else {
    stop("Please, indicate an onemap object or a list of sequences using onemap.obj and list.seqs arguments.")
  }
}

#' Sort markers in onemap object by their position in reference genome
#' 
#' @param onemap.obj object of class onemap
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
##' individual.}
#' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#' 
#' @export
sort_by_pos <- function(onemap.obj){
  if(!inherits(onemap.obj, "onemap")) stop("Input must to be of onemap class \n")
  
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

#' Produce empty object to avoid code break. Function for internal purpose.
#'  
#' @param vcf object of class vcfR
#' @param P1 character with parent 1 ID
#' @param P2 character with parent 2 ID
#' @param cross type of cross. Must be one of: \code{"outcross"} for full-sibs;
#' \code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
#' \code{"ri self"} for recombinant inbred lines by self-mating; or
#' \code{"ri sib"} for recombinant inbred lines by sib-mating.
#' 
##' @return An empty object of class \code{onemap}, i.e., a list with the following
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
##' individual.}
#' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' 
#' @export
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

#' Remove duplicated markers keeping the one with less missing data
#'  
#' @param onemap.obj object of class \code{onemap}
#'  
##' @return An empty object of class \code{onemap}, i.e., a list with the following
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
##' individual.}
#' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' 
#'  
#' @export
rm_dupli_mks <- function(onemap.obj){
  
  if(!inherits(onemap.obj, c("onemap"))) stop("Input object must be of class onemap")
  
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
      if(inherits(discard, "matrix")){
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
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
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
  test[7] <- if(inherits(x, "f2")) {
    !all(unique(x$segr.type) %in% c("A.H.B", "D.B", "C.A"))
  } else if(inherits(x, "outcross")){
    !all(unique(x$segr.type) %in% c("A.1", "A.2", "A.3", "A.4", "B1.5", 
                                    "B2.6", "B3.7", "C.8", "D1.9", "D1.10", 
                                    "D1.11", "D1.12", "D1.13", "D2.14", 
                                    "D2.15", "D2.16", "D2.17", "D2.18"))
  } else if(inherits(x, "backcross")){
    !all(unique(x$segr.type) %in% c("A.H"))
  } else if(inherits(x, "riself") | inherits(x, "risib")){
    !all(unique(x$segr.type) %in% c("A.B"))
  }
  
  test[8] <- if(inherits(x, "f2")) {
    !all(unique(x$segr.type.num) %in% c(4,5))
  } else if(inherits(x, "outcross")){
    !all(unique(x$segr.type.num) %in% 1:7)
  } else if(inherits(x, "backcross")){
    !all(unique(x$segr.type.num) %in% 8)
  } else if(inherits(x, "riself") | inherits(x, "risib")){
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
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#' 
#' @export
check_twopts <- function(x){
  test <- logical(4L)
  names(test) <- 1:4
  
  test[1] <- if(check_data(x$data.name) == 0) FALSE else TRUE
  if(inherits(x$data.name, "outcross") | inherits(x$data.name, "f2")){
    test[2] <- !inherits(x$analysis, "list")
    test[3] <- all(dim(x$analysis[[1]]) != rep(x$data.name$n.mar,2))
    test[4] <- any(sapply(x$analysis, function(x) any(is.na(x))))
  } else {
    test[2] <- !inherits(x$analysis, "matrix")
    test[3] <- all(dim(x$analysis) != rep(x$data.name$n.mar,2))
    test[4] <- any(is.na(x$analysis))
  }
  
  if(any(test))
    return(test)
  else 
    return(0)
}


#' Filter markers based on 2pts distance
#' 
#' @param input.seq object of class sequence with ordered markers
#' 
#' @param max.gap maximum gap measured in kosambi centimorgans allowed between adjacent markers. 
#' Markers that presents the defined distance between both adjacent neighbors will be removed.
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
#' @export
filter_2pts_gaps <- function(input.seq, max.gap=10){
  
  if(!inherits(input.seq, "sequence"))
    stop("input.seq object must be of class sequence.")
  
  ## extracting data
  if(inherits(input.seq$data.name, "outcross") | inherits(input.seq$data.name, "f2"))
  {
    ## making a list with necessary information
    n.mrk <- length(input.seq$seq.num)
    LOD <- lapply(input.seq$twopt$analysis,
                  function(x, w){
                    m<-matrix(0, length(w), length(w))
                    for(i in 1:(length(w)-1)){
                      for(j in (i+1):length(w)){
                        z<-sort(c(w[i],w[j]))
                        m[j,i]<-m[i,j]<-x[z[1], z[2]]
                      }
                    }
                    return(m)
                  }, input.seq$seq.num
    )
    mat<-t(get_mat_rf_out(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
  } else {
    ## making a list with necessary information
    n.mrk <- length(input.seq$seq.num) 
    LOD<-matrix(0, length(input.seq$seq.num), length(input.seq$seq.num))
    for(i in 1:(length(input.seq$seq.num)-1)){
      for(j in (i+1):length(input.seq$seq.num)){
        z<-sort(c(input.seq$seq.num[i],input.seq$seq.num[j]))
        LOD[j,i]<-LOD[i,j]<-input.seq$twopt$analysis[z[1], z[2]]
      }
    }
    mat<-t(get_mat_rf_in(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
  }
  
  dist.LOD <- dist.rf <- vector()
  for(i in 1:dim(mat)[1]-1){
    dist.LOD <- c(dist.LOD, mat[i, i+1])
    dist.rf <- c(dist.rf, round(kosambi(mat[i+1, i]),2))
  }
  
  idx <- which(dist.rf > max.gap)
  rm.seq <- vector()
  for(i in 1:length(idx)){
    if(idx[i] == 1){
      rm.seq <- c(rm.seq, 1)
    } else if(idx[i] == length(dist.rf) -1){
      rm.seq <- c(rm.seq, idx[i] + 1)
      if(idx[i-1] == idx[i] -1)
        rm.seq <- c(rm.seq, idx[i])
    } else if(i-1 != 0){
      if(idx[i-1] == idx[i] -1) {
        rm.seq <- c(rm.seq, idx[i])
      }
    }
  }
  
  if(length(rm.seq) > 0) new.seq <- make_seq(input.seq$twopt, input.seq$seq.num[-rm.seq]) else new.seq <- input.seq
  
  return(new.seq)
}


##' Filter markers according with a two-points recombination fraction and LOD threshold. Adapted from MAPpoly.
##'
##' @param input.seq an object of class \code{onemap}.
#' @param thresh.LOD.rf LOD score threshold for recombination fraction (default = 5)
#' @param thresh.rf threshold for recombination fractions (default = 0.15)
#' @param probs indicates the probability corresponding to the filtering quantiles. (default = c(0.05, 1))
##'
##' @return An object of class \code{sequence}, which is a list containing the
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
##'
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' @examples
##'
##'  data("vcf_example_out")
##'  twopts <- rf_2pts(vcf_example_out)
##'  seq1 <- make_seq(twopts, which(vcf_example_out$CHROM == "1"))
##' filt_seq <- rf_snp_filter_onemap(seq1, 20, 0.5, c(0.5,1))
##'
##'@export
rf_snp_filter_onemap <- function(input.seq, thresh.LOD.rf = 5, thresh.rf = 0.15, probs = c(0.05,1)){
  if(inherits(input.seq$data.name, "outcross") | inherits(input.seq$data.name, "f2")){
    rf.mat <- get_mat_rf_out(input.seq, LOD = T)
  } else {
    rf.mat <- get_mat_rf_in(input.seq, LOD = T)
  }
  rf.mat[lower.tri(rf.mat)][rf.mat[lower.tri(rf.mat)]  <= thresh.LOD.rf] <- NA
  rf.mat[upper.tri(rf.mat)][rf.mat[upper.tri(rf.mat)]  >= thresh.rf] <- NA
  x <- apply(rf.mat, 1, function(x) sum(!is.na(x)))
  th <- quantile(x, probs = probs)
  rem <- c(which(x < th[1]), which(x > th[2]))
  ids <- names(which(x >= th[1] & x <= th[2]))
  new.seq <- make_seq(input.seq$twopt, match(ids, colnames(input.seq$data.name$geno)))
  return(new.seq)
}

##'   Creates a new sequence by adding markers.
##'
##'   Creates a new sequence by adding markers from a predetermined
##'   one. The markers are added in the end of the sequence.
##'
##' @param input.seq an object of class \code{sequence}.
##'
##' @param mrks a vector containing the markers to be added from the
##'     \code{sequence}.
##'
##' @return An object of class \code{sequence}, which is a list
##'     containing the following components:
##'
##' \item{seq.num}{a \code{vector} containing the (ordered) indices of
##'     markers in the sequence, according to the input file.}
##'
##' \item{seq.phases}{a \code{vector} with the linkage phases between
##'     markers in the sequence, in corresponding positions. \code{-1}
##'     means that there are no defined linkage phases.}
##'
##' \item{seq.rf}{a \code{vector} with the recombination fractions
##'     between markers in the sequence. \code{-1} means that there
##'     are no estimated recombination fractions.}
##'
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##'     \item{data.name}{name of the object of class \code{onemap}
##'     with the raw data.}
##'
##' \item{twopt}{name of the object of class \code{rf_2pts} with the
##'     2-point analyses.}
##'
##'  @author Marcelo Mollinari, \email{mmollina@@usp.br}
##'
##' @seealso \code{\link[onemap]{drop_marker}}
##'
##' @examples
##' data(onemap_example_out)
##' twopt <- rf_2pts(onemap_example_out)
##' all_mark <- make_seq(twopt,"all")
##' groups <- group(all_mark)
##' (LG1 <- make_seq(groups,1))
##' (LG.aug<-add_marker(LG1, c(4,7)))
##'
##' @export
add_marker<-function(input.seq, mrks)
{
  if (!inherits(input.seq,"sequence"))
    stop(sQuote(deparse(substitute(input.seq))), " is not an object of class 'sequence'")
  seq.num<-c(input.seq$seq.num,mrks)
  return(make_seq(input.seq$twopt,seq.num, twopt=input.seq$twopt))
}


##' Keep in the onemap and twopts object only markers in the sequences
##' 
##' @param list.sequences a list of objects 'sequence'
##' 
##' @return a list of objects 'sequences' with internal onemap and twopts objects reduced
##' 
##' @author Cristiane Taniguti
##' 
##' @export
keep_only_selected_mks <- function(list.sequences= NULL){
  if(!inherits(list.sequences, "list")) stop("Object is not a list")
  if(!all(sapply(list.sequences, function(x) inherits(x, "sequence")))) stop("One or more of the list components is/are not of class sequence")
  mk.numbers <- sapply(list.sequences, function(x) x$seq.num)
  mk.names <- sapply(list.sequences, function(x) colnames(x$data.name$geno)[x$seq.num])
  mk.numbers <- unlist(mk.numbers)
  
  new_onemap <- split_onemap(list.sequences[[1]]$data.name, unique(mk.numbers))
  new_twopts <- rf_2pts(new_onemap, verbose = FALSE)
  
  new_seqs <- list()
  for(i in 1:length(mk.names)){
    new_seqs[[i]] <- make_seq(new_twopts, match(mk.names[[i]], colnames(new_twopts$data.name$geno)))
  }
  return(new_seqs)
}


#' Order the markers in a sequence using the genomic position
#' 
#' @param input.seq object of class `sequence`
#' 
#' @return An object of class \code{sequence}
#' 
#' @author Cristiane Taniguti
#' 
#' @export
ord_by_geno <- function(input.seq){
  
  if(!(inherits(input.seq,c("sequence")))) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(is.null(input.seq$data.name$CHROM) | is.null(input.seq$data.name$POS)) stop("Reference genome chromosome and position information are not available for this dataset.")
  
  ord.seq <- input.seq$seq.num[order(input.seq$data.name$CHROM[input.seq$seq.num],input.seq$data.name$POS[input.seq$seq.num])]
  new.seq <- make_seq(input.seq$twopt, ord.seq)
  return(new.seq)
}

