#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: create_probs.R                                                ##
## Contains: create_probs                                              ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
## copyright (c) 2019, Cristiane Taniguti                              ##
##                                                                     ##
## First version: 28/06/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################
## This function builds matrix with genotypes probabilites for emission function in hmm approach

##' Build genotype probabilities matrix for hmm
##'
##' The genotypes probabilities can be calculated considering a global error (default method)
##' or considering a genotype error probability for each genotype. Furthermore, user can provide 
##' directly the genotype probability matrix.
##' 
##' The genotype probability matrix has number of individuals x number of markers rows and
##' four columns (or two if considering backcross or RILs populations), one for each possible genotype
##' of the population. This format follows the one proposed by MAPpoly.
##' 
##' The genotype probabilities come from SNP calling methods. If you do not have them, you can use a global
##' error or a error value for each genotype. The OneMap until 2.1 version have only the global error option.
##' 
##' @param input.obj object of class onemap or onemap sequence
##' @param global_error a integer specifing the global error value
##' @param genotypes_errors a matrix with dimensions (number of individuals) x (number of markers) with genotypes errors values
##' @param genotypes_probs a matrix with dimensions (number of individuals)*(number of markers) x possible genotypes 
##' (i.e., a ab ba b) with four columns for f2 and outcrossing populations, and two for backcross and RILs).
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
##' @author Cristiane Taniguti \email{chtaniguti@@tamu.edu} 
##' @seealso \code{\link[onemap]{make_seq}}
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' @examples
##' 
##'   data(onemap_example_out)
##'   new.data <- create_probs(onemap_example_out, global_error = 10^-5)
##'   
##'@importFrom reshape2 melt
##'
##'@export
create_probs <- function(input.obj = NULL, 
                         global_error = NULL, 
                         genotypes_errors = NULL, 
                         genotypes_probs = NULL){
  
  if(!(inherits(input.obj,c("onemap", "sequence")))){
    stop("input.obj should be of class onemap or sequence")
  }
  
  if(inherits(input.obj, "sequence")) {
    seq.obj <- input.obj
    input.obj <- input.obj$data.name
    flag <- TRUE
  } else {
    flag <- FALSE
  }
  
  # Empty object
  if(input.obj$n.mar == 0){
    warning("It is a empty onemap object. Nothing will be done.")
    return(input.obj)
  }
  
  if(all(is.null(c(global_error, genotypes_errors, genotypes_probs)))){
    global_error <- 10^-5
  }
  
  crosstype <- class(input.obj)[2]
  
  probs <- melt(t(input.obj$geno))
  probs$type <- rep(input.obj$segr.type.num, input.obj$n.ind)
  
  if(!is.null(global_error) | !is.null(genotypes_errors)){
    
    if(!is.null(global_error)) {
      error <- rep(global_error, length(probs$value))
    } else {
      # checks
      if(!all(colnames(input.obj$geno)%in%colnames(genotypes_errors))){
        stop("Not all markers in onemap object have corresponding genotype errors in matrix")
      }
      
      if(!all(colnames(genotypes_errors)%in%colnames(input.obj$geno))){
        stop("There are more markers in errors matrix than in onemap object")
      }
      
      if(!all(rownames(input.obj$geno)%in%rownames(genotypes_errors))){
        stop("Not all individuals in onemap object have corresponding genotype errors in matrix")
      }
      
      if(!all(rownames(input.obj$geno)%in%rownames(genotypes_errors))){
        stop("There are more individuals in errors matrix than in onemap object")
      }
      
      error <- melt(t(genotypes_errors))
      error <- error$value
    }
    
    # Check if is working with f2
    if(crosstype == "outcross" | crosstype == "f2"){
      prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      idx <- which(is.na(error))
      prob[idx,] <- 1
      # A
      idx <- which(probs$value == 1  & probs$type == 1)
      prob[idx,] <- c(1-error[idx], rep(error[idx]/3,3))
      idx <- which(probs$value == 2  & probs$type == 1)
      prob[idx,] <- c(error[idx]/3, 1-error[idx], rep(error[idx]/3,2))
      idx <- which(probs$value == 3  & probs$type == 1)
      prob[idx,] <- c(rep(error[idx]/3,2), 1-error[idx], error[idx]/3)
      idx <- which(probs$value == 4  & probs$type == 1)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # B1
      idx <- which(probs$value == 1  & probs$type == 2)
      prob[idx,] <- c(rep(1-error[idx],2), rep(error[idx],2))
      idx <- which(probs$value == 2  & probs$type == 2)
      prob[idx,] <- c(rep(error[idx]/3,2), 1-error[idx], error[idx]/3)
      idx <- which(probs$value == 3  & probs$type == 2)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # B2
      idx <- which(probs$value == 1  & probs$type == 3)
      prob[idx,] <- c(1-error[idx], error[idx], 1-error[idx], error[idx])
      idx <- which(probs$value == 2  & probs$type == 3)
      prob[idx,] <- c(error[idx]/3, 1-error[idx], rep(error[idx]/3,2))
      idx <- which(probs$value == 3  & probs$type == 3)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # B3.7
      idx <- which(probs$value == 1  & probs$type == 4)
      prob[idx,] <- c(1-error[idx], rep(error[idx]/3,3))
      idx <- which(probs$value == 2  & probs$type == 4)
      prob[idx,] <- c(error[idx], rep(1-error[idx],2), error[idx])
      idx <- which(probs$value == 3  & probs$type == 4)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # C
      idx <- which(probs$value == 1  & probs$type == 5)
      prob[idx,] <- c(rep((1-error[idx])/3,3), error[idx])
      idx <- which(probs$value == 2  & probs$type == 5)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # D1
      idx <- which(probs$value == 1  & probs$type == 6)
      prob[idx,] <- c(rep(1-error[idx],2), rep(error[idx],2))
      idx <- which(probs$value == 2  & probs$type == 6)
      prob[idx,] <- c(rep(error[idx],2), rep(1-error[idx],2))
      idx <- which(probs$value == 3  & probs$type == 6)
      prob[idx,] <- 1
      
      # D2
      idx <- which(probs$value == 1  & probs$type == 7)
      prob[idx,] <- c(1-error[idx], error[idx], 1-error[idx], error[idx])
      idx <- which(probs$value == 2  & probs$type == 7)
      prob[idx,] <- c(error[idx], 1-error[idx], error[idx], 1-error[idx])
      idx <- which(probs$value == 3  & probs$type == 7)
      prob[idx,] <- 1
      
      # } else if(crosstype == "f2"){
      #   prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      #   idx <- which(probs$value == 0)
      #   prob[idx,] <- 1
      #   idx <- which(probs$value == 1)
      #   prob[idx,] <- c(1- error[idx], rep(error[idx]/3,3))
      #   idx <- which(probs$value == 2)
      #   prob[idx,] <- c(error[idx]/2, rep(1-error[idx],2), error[idx]/2)
      #   idx <- which(probs$value == 3)
      #   prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      #   idx <- which(probs$value == 4)
      #   prob[idx,] <- c(rep(1-error[idx]/3,3), error[idx])
      #   idx <- which(probs$value == 5)
      #   prob[idx,] <- c(error[idx], rep(1-error[idx]/3,3))
    } else if(crosstype == "backcross" | crosstype == "riself" | crosstype == "risib"){
      prob <- matrix(NA, nrow=length(probs$value), ncol = 2)
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      idx <- which(probs$value == 1)
      prob[idx,] <- c(1- error[idx], error[idx])
      idx <- which(probs$value == 2)
      prob[idx,] <- c(error[idx], 1-error[idx])
      idx <- which(probs$value == 3)
      prob[idx,] <- c(error[idx], 1-error[idx])
    } 
    
  }
  
  if(!is.null(genotypes_probs)){
    
    # Only for biallelic markers codominant markers
    if(crosstype == "outcross" | crosstype == "f2"){
      # Sometimes the 1 and 3 are inverted
      # D2.15 when P2 are heterozygote receives 1 instead of 3
      # D1.10 when P1 are heterozygote receives 1 instead of 3
      prob.temp <- matrix(NA, nrow=length(probs$value), ncol = 3)
      
      prob.temp[,2] <- genotypes_probs[,2]
      het.idx <- which(probs$value == 2)
      
      # The homozygotes can be inverted in heterozygotes genotypes
      hom1.idx <- which(probs$value == 1)
      # If it has only one marker the apply breaks
      if(length(hom1.idx) > 1){
        sub <- genotypes_probs[hom1.idx,-2]
      } else {
        sub <- t(as.matrix(genotypes_probs[hom1.idx,-2]))
      }
      
      for_het <- table(apply(sub,1, which.max)) 
      for_het <- as.numeric(which.max(for_het))
      
      if(for_het == 2){
        prob.temp[het.idx,3] <- genotypes_probs[het.idx,1]
        prob.temp[het.idx,1] <- genotypes_probs[het.idx,3]
      } else {
        prob.temp[het.idx,1] <- genotypes_probs[het.idx,1]
        prob.temp[het.idx,3] <- genotypes_probs[het.idx,3]
      }
      
      # We consider that the genotype number is correct and the probability need to be the max at the genotype column
      # We change the position of the wrong probabilities
      hom1.idx.prob <- which(apply(sub,1, which.max) == 1)
      prob.temp[hom1.idx[hom1.idx.prob],1] <- genotypes_probs[hom1.idx[hom1.idx.prob], 1]
      prob.temp[hom1.idx[hom1.idx.prob],3] <- genotypes_probs[hom1.idx[hom1.idx.prob], 3]
      hom1.idx.prob <- which(apply(sub,1, which.max) == 2)
      prob.temp[hom1.idx[hom1.idx.prob],1] <- genotypes_probs[hom1.idx[hom1.idx.prob], 3]
      prob.temp[hom1.idx[hom1.idx.prob],3] <- genotypes_probs[hom1.idx[hom1.idx.prob], 1]
      
      hom3.idx <- which(probs$value == 3)
      if(length(hom3.idx) > 1){
        sub <- genotypes_probs[hom3.idx,-2]
      } else {
        sub <- t(as.matrix(genotypes_probs[hom3.idx,-2]))
      }
      
      hom3.idx.prob <- which(apply(sub,1, which.max) == 1)
      prob.temp[hom3.idx[hom3.idx.prob],3] <- genotypes_probs[hom3.idx[hom3.idx.prob], 1]
      prob.temp[hom3.idx[hom3.idx.prob],1] <- genotypes_probs[hom3.idx[hom3.idx.prob], 3]
      hom3.idx.prob <- which(apply(sub,1, which.max) == 2)
      prob.temp[hom3.idx[hom3.idx.prob],3] <- genotypes_probs[hom3.idx[hom3.idx.prob], 3]
      prob.temp[hom3.idx[hom3.idx.prob],1] <- genotypes_probs[hom3.idx[hom3.idx.prob], 1]
      
      prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      
      # B3.7
      idx <- which(probs$type == 4)
      prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,2], prob.temp[idx,2], prob.temp[idx,3])
      
      # D1 
      idx <- which(probs$type == 6)
      prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,1], prob.temp[idx,2], prob.temp[idx,2])
      
      # D2
      idx <- which(probs$type == 7)
      prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,2], prob.temp[idx,1], prob.temp[idx,2])
      
      # Missing data -- all genotypes 0 will receive 1 for all possible genotypes probabilities
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      
    } else {
      prob.temp <- matrix(NA, nrow=length(probs$value), ncol = 3)
      # Sometimes the 1 and 3 are inverted
      prob.temp[,2] <- genotypes_probs[,2]
      het.idx <- which(probs$value == 2)
      
      hom1.idx <- which(probs$value == 1)
      # If it has only one marker the apply breaks
      if(length(hom1.idx) > 1){
        sub <- genotypes_probs[hom1.idx,-2]
      } else {
        sub <- t(as.matrix(genotypes_probs[hom1.idx,-2]))
      }
      
      hom1.idx.prob <- unique(apply(genotypes_probs[hom1.idx,],1, which.max))
      prob.temp[hom1.idx,1] <- genotypes_probs[hom1.idx, hom1.idx.prob]
      if(hom1.idx.prob == 3){
        prob.temp[hom1.idx,3] <- genotypes_probs[hom1.idx, 1]
        prob.temp[het.idx,3] <- genotypes_probs[het.idx,1]
        prob.temp[het.idx,1] <- genotypes_probs[het.idx,3]
      } else {
        prob.temp[hom1.idx,3] <- genotypes_probs[hom1.idx, 3]
        prob.temp[het.idx,1] <- genotypes_probs[het.idx,1]
        prob.temp[het.idx,3] <- genotypes_probs[het.idx,3]
      }
      
      hom3.idx <- which(probs$value == 3)
      # If it has only one marker the apply breaks
      if(length(hom3.idx) > 1){
        sub <- genotypes_probs[hom3.idx,-2]
      } else {
        sub <- t(as.matrix(genotypes_probs[hom3.idx,-2]))
      }
      
      hom3.idx.prob <- unique(apply(genotypes_probs[hom3.idx,],1, which.max))
      prob.temp[hom3.idx,3] <- genotypes_probs[hom3.idx, hom3.idx.prob]
      if(hom3.idx.prob == 3){
        prob.temp[hom3.idx,1] <- genotypes_probs[hom3.idx, 1]
      } else {
        prob.temp[hom3.idx,1] <- genotypes_probs[hom3.idx, 3]
      }
      # if(crosstype == "f2"){
      #   prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      #   prob <- cbind(prob.temp[,1], prob.temp[,2], prob.temp[,2], prob.temp[,3])
      #   
      #   idx <- which(probs$value == 0)
      #   prob[idx,] <- 1
      #   
      # } else {
      if(crosstype == "backcross" | crosstype == "riself" | crosstype == "risib") {
        prob <- matrix(NA, nrow=length(probs$value), ncol = 2)
        
        idx <- which(probs$value == 1)
        prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,2])
        idx <- which(probs$value == 2)
        prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,2])
        idx <- which(probs$value == 3)
        prob[idx,] <- cbind(prob.temp[idx,1], prob.temp[idx,3])
        
        idx <- which(probs$value == 0)
        prob[idx,] <- 1
      }
    }
    # The values of estimated genotypes must sum 1
  }
  
  rownames(prob) <- paste0(probs$Var1, "_", probs$Var2)
  
  input.obj$error<- t(apply(prob, 1, function(x) if(all(is.na(x))) rep(1, 4) else x))
  
  if(flag) {
    seq.obj$data.name <- input.obj
    seq.obj$twopt$data.name <- input.obj
    input.obj <- seq.obj
  }
  
  return(input.obj)
}