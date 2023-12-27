#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: filters.R                                                     ##
## Contains: filter_missing filter_prob edit_order                     ##
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
##' @param by character defining if `markers` or `individuals` should be filtered
##' @param verbose A logical, if TRUE it output progress status
##' information.
##' 
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
##' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu} 
##' @examples
##' 
##'   data(onemap_example_out)
##'   filt_obj <- filter_missing(onemap_example_out, threshold=0.25)
##'  
##'@export
filter_missing <- function(onemap.obj=NULL, threshold= 0.25, by = "markers", verbose = TRUE){
  if(!inherits(onemap.obj,"onemap")){
    stop("onemap.obj should be of class onemap\n")
  }
  
  if(by == "markers"){
    perc.mis <- apply(onemap.obj$geno, 2, function(x) sum(x == 0)/length(x))
    idx <- which(!perc.mis > threshold)
    
    new.onemap.obj <- onemap.obj
    new.onemap.obj$geno <- onemap.obj$geno[,idx]
    new.onemap.obj$n.mar <- length(idx)
    new.onemap.obj$segr.type <- onemap.obj$segr.type[idx]
    new.onemap.obj$segr.type.num <- onemap.obj$segr.type.num[idx]
    if(!is.null(onemap.obj$CHROM)) new.onemap.obj$CHROM <- onemap.obj$CHROM[idx]
    if(!is.null(onemap.obj$POS)) new.onemap.obj$POS <- onemap.obj$POS[idx]
    if(!is.null(onemap.obj$ref_alt_alleles)) new.onemap.obj$ref_alt_alleles <- onemap.obj$ref_alt_alleles[idx,]
    new.onemap.obj$error <- onemap.obj$error[idx + rep(c(0:(onemap.obj$n.ind-1))*onemap.obj$n.mar, each=length(idx)),]
    if(verbose) cat("Number of markers removed from the onemap object: ", length(which(perc.mis > threshold)), "\n")
  } else if (by == "individuals"){
    perc.mis <- apply(onemap.obj$geno, 1, function(x) sum(x == 0)/length(x))
    idx <- which(!perc.mis > threshold)
    new.onemap.obj <- onemap.obj
    new.onemap.obj$geno <- onemap.obj$geno[idx,]
    new.onemap.obj$n.ind <- length(idx)
    new.onemap.obj$error <- onemap.obj$error[1:onemap.obj$n.mar + rep((idx-1)*onemap.obj$n.mar, each=onemap.obj$n.mar),]
    if(verbose) cat("Number of indiduals removed from the onemap object: ", length(which(perc.mis > threshold)), "\n")
  } else {
    stop("Input for argument by is not defined. Please choose between `markers` or `individuals` options.")
  }
  return(new.onemap.obj)
}


##' Function filter genotypes by genotype probability
##'
##' @param onemap.obj an object of class \code{onemap}.
##' @param threshold a numeric from 0 to 1 to define the threshold for 
##' the probability of the called genotype (highest probability)
##' @param verbose If \code{TRUE}, print tracing information.
##' 
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
##' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu} 
##' @examples
##' \donttest{
##'   data(onemap_example_out)
##'   filt_obj <- filter_prob(onemap_example_out, threshold=0.8)
##'  }
##' @importFrom reshape2 melt dcast
##' 
##' @export
filter_prob <- function(onemap.obj=NULL, threshold= 0.8, verbose=TRUE){
  idx <- apply(onemap.obj$error, 1, which.max)
  rm <- which(onemap.obj$error[cbind(seq_along(idx), idx)] < threshold)
  onemap.obj$error[rm,] <- 1
  if(verbose) cat(paste(length(rm), "genotypes were converted to missing data."))
  
  geno_melt <- melt(onemap.obj$geno)
  geno_melt[rm,3] <- 0
  geno <- dcast(geno_melt, Var1 ~ Var2)
  rownames(geno) <- geno$Var1
  geno <- geno[,-1]
  
  onemap.obj$geno <- as.matrix(geno)
  
  return(onemap.obj)
}


#' Edit sequence ordered by reference genome positions 
#' comparing to another set order
#' 
#' @param input.seq object of class sequence with alternative order (not genomic order)
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
#' @export
edit_order_onemap <- function(input.seq){
  
  if (!inherits(input.seq, "sequence")) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'sequence'")
  }
  
  if(length(unique(input.seq$data.name$CHROM[input.seq$seq.num])) > 1) stop("There are markers from more than one chromosome in the sequence.")
  
  get_weird <- data.frame(x = 1:length(input.seq$seq.num),
                          y = input.seq$data.name$POS[input.seq$seq.num])
  
  rownames(get_weird) <- colnames(input.seq$data.name$geno)[input.seq$seq.num]
  get_weird <- get_weird[order(get_weird$y),]
  plot(get_weird$x, get_weird$y, xlab="alternative order", ylab = "Genome position")
  
  inverted <- removed <- vector()
  if(interactive()){
    ANSWER <- "Y"
    while(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == ""){
      plot(get_weird$x, get_weird$y, xlab="sequence order", ylab = "Genome position")
      mks.to.remove <- gatepoints::fhs(get_weird, mark = TRUE)
      if(length(which(rownames(get_weird) %in% mks.to.remove)) > 0){
        ANSWER2 <- readline("Enter 'invert/remove' to proceed with the edition: ")
        if(ANSWER2 == "invert"){
          inverted <- c(inverted, as.vector(mks.to.remove))
          repl <- get_weird[rev(which(rownames(get_weird) %in% as.vector(mks.to.remove))),]
          get_weird[which(rownames(get_weird) %in% as.vector(mks.to.remove)),2] <- repl[,2]
          rownames(get_weird)[which(rownames(get_weird) %in% as.vector(mks.to.remove))] <- rownames(repl)
        } else {
          removed <- c(removed, as.vector(mks.to.remove))
          get_weird <- get_weird[-which(rownames(get_weird) %in% mks.to.remove),]
        }
      }
      ANSWER <- readline("Enter 'Y/n' to proceed with interactive edition or quit: ")
    }
    plot(get_weird$x, get_weird$y, xlab="sequence order", ylab = "Genome position")
  }
  
  return(structure(list(edited_order = rownames(get_weird),
                        removed = removed,
                        inverted = inverted,
                        data.name = input.seq$data.name,
                        twopts = input.seq$twopt), class = "onemap.edit.order"))
}

