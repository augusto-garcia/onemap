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

##'Converts probabilities to onemap format
##'
##'@param onemap.obj define
##'@param global_error define
##'@param genotypes_errors define
##'@param genotypes_probs define
##'
##'@importFrom reshape2 melt
##'
##'@export
create_probs <- function(onemap.obj = NULL, 
                         global_error = NULL, 
                         genotypes_errors = NULL, 
                         genotypes_probs = NULL){
  
  if(class(onemap.obj)[1]!="onemap"){
    stop("onemap.obj should be of class onemap")
  }
  
  # Empty object
  if(onemap.obj$n.mar == 0){
    warning("It is a empty onemap object. Nothing will be done.")
    return(onemap.obj)
  }
  
  if(all(is.null(c(global_error, genotypes_errors, genotypes_probs)))){
    global_error <- 10^-5
  }
  
  crosstype <- class(onemap.obj)[2]
  
  probs <- melt(t(onemap.obj$geno))
  probs$type <- rep(onemap.obj$segr.type.num, onemap.obj$n.ind)
  
  if(!is.null(global_error) | !is.null(genotypes_errors)){
    
    if(!is.null(global_error)) {
      error <- rep(global_error, length(probs$value))
    } else {
      # checks
      if(!all(colnames(onemap.obj$geno)%in%colnames(genotypes_errors))){
        stop("Not all markers in onemap object have corresponding genotype errors in matrix")
      }
      
      if(!all(colnames(genotypes_errors)%in%colnames(onemap.obj$geno))){
        stop("There are more markers in errors matrix than in onemap object")
      }
      
      if(!all(rownames(onemap.obj$geno)%in%rownames(genotypes_errors))){
        stop("Not all individuals in onemap object have corresponding genotype errors in matrix")
      }
      
      if(!all(rownames(onemap.obj$geno)%in%rownames(genotypes_errors))){
        stop("There are more individuals in errors matrix than in onemap object")
      }
      
      error <- melt(t(genotypes_errors))
      error <- error$value
    }
    
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
      prob[idx,] <- c(rep(1-error[idx],2), rep(error[idx]/2,2))
      idx <- which(probs$value == 2  & probs$type == 2)
      prob[idx,] <- c(rep(error[idx]/3,2), 1-error[idx], error[idx]/3)
      idx <- which(probs$value == 3  & probs$type == 2)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # B2
      idx <- which(probs$value == 1  & probs$type == 3)
      prob[idx,] <- c(1-error[idx], error[idx]/2, 1-error[idx], error[idx]/2)
      idx <- which(probs$value == 2  & probs$type == 3)
      prob[idx,] <- c(error[idx]/3, 1-error[idx], rep(error[idx]/3,2))
      idx <- which(probs$value == 3  & probs$type == 3)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # B3.7
      idx <- which(probs$value == 1  & probs$type == 4)
      prob[idx,] <- c(1-error[idx], rep(error[idx]/3,3))
      idx <- which(probs$value == 2  & probs$type == 4)
      prob[idx,] <- c(error[idx]/2, rep(1-error[idx],2), error[idx]/2)
      idx <- which(probs$value == 3  & probs$type == 4)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # C
      idx <- which(probs$value == 1  & probs$type == 5)
      prob[idx,] <- c(rep(1-error[idx],3), error[idx])
      idx <- which(probs$value == 2  & probs$type == 5)
      prob[idx,] <- c(rep(error[idx]/3,3), 1-error[idx])
      
      # D1
      idx <- which(probs$value == 1  & probs$type == 6)
      prob[idx,] <- c(rep(1-error[idx],2), rep(error[idx]/2,2))
      idx <- which(probs$value == 2  & probs$type == 6)
      prob[idx,] <- c(rep(error[idx]/2,2), rep(1-error[idx],2))
      idx <- which(probs$value == 3  & probs$type == 6)
      prob[idx,] <- 1
      
      # D2
      idx <- which(probs$value == 1  & probs$type == 7)
      prob[idx,] <- c(1-error[idx], error[idx]/2, 1-error[idx], error[idx]/2)
      idx <- which(probs$value == 2  & probs$type == 7)
      prob[idx,] <- c(error[idx]/2, 1-error[idx], error[idx]/2, 1-error[idx])
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
      for_het <- table(apply(genotypes_probs[hom1.idx,-2],1, which.max)) 
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
      hom1.idx.prob <- which(apply(genotypes_probs[hom1.idx,-2],1, which.max) == 1)
      prob.temp[hom1.idx[hom1.idx.prob],1] <- genotypes_probs[hom1.idx[hom1.idx.prob], 1]
      prob.temp[hom1.idx[hom1.idx.prob],3] <- genotypes_probs[hom1.idx[hom1.idx.prob], 3]
      hom1.idx.prob <- which(apply(genotypes_probs[hom1.idx,-2],1, which.max) == 2)
      prob.temp[hom1.idx[hom1.idx.prob],1] <- genotypes_probs[hom1.idx[hom1.idx.prob], 3]
      prob.temp[hom1.idx[hom1.idx.prob],3] <- genotypes_probs[hom1.idx[hom1.idx.prob], 1]
      
      hom3.idx <- which(probs$value == 3)
      hom3.idx.prob <- which(apply(genotypes_probs[hom3.idx,-2],1, which.max) == 1)
      prob.temp[hom3.idx[hom3.idx.prob],3] <- genotypes_probs[hom3.idx[hom3.idx.prob], 1]
      prob.temp[hom3.idx[hom3.idx.prob],1] <- genotypes_probs[hom3.idx[hom3.idx.prob], 3]
      hom3.idx.prob <- which(apply(genotypes_probs[hom3.idx,-2],1, which.max) == 2)
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
  }
  
  rownames(prob) <- paste0(probs$Var1, "_", probs$Var2)
  
  onemap.obj$error <- prob
  
  return(onemap.obj)
}


