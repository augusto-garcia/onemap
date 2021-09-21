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
##'@param input.obj object of class onemap or onemap sequence
##'@param global_error single value to be considered as error probability in HMM emission function
##'@param genotypes_errors matrix individuals x markers with error values for each marker
##'@param genotypes_probs table containing the probability distribution for each combination of marker × individual. 
##'Each line on this table represents the combination of one marker with one individual, and the respective probabilities.
##'The table should contain four three columns (prob(AA), prob(AB) and prob(BB)) and individuals x markers rows.
##'
##'@importFrom reshape2 melt
##'
##'@export
create_probs <- function(input.obj = NULL, 
                         global_error = NULL, 
                         genotypes_errors = NULL, 
                         genotypes_probs = NULL){
  
  if(!(is(input.obj,"onemap") | is(input.obj,"sequence"))){
    stop("input.obj should be of class onemap or sequence")
  }
  
  if(is(input.obj, "sequence")) {
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
  
  input.obj$error <- prob
  
  if(flag) {
    seq.obj$data.name <- input.obj
    seq.obj$twopt$data.name <- input.obj
    input.obj <- seq.obj
  }
  
  return(input.obj)
}