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
##' @param onemap.obj an object of class \code{onemap}.
##' @param prob a integer specifing the global error value, or a matrix with dimensions
##' (number of marker) x (number of markers) with genotypes errors values, or a matrix
##' with dimensions (number of individuals)*(number of markers) x possible genotypes (i.e., a ab ba b)
##' with four columns for f2 and outcrossing populations, and two for backcross and RILs).
##' 
##' @return An object of class \code{onemap} with the built matrix at prob component of the list
##' @author Cristiane Taniguti \email{chtaniguti@@usp.br} 
##' @seealso \code{\link[onemap]{make_seq}}
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' @examples
##'
##'   data(onemap_example_out)
##'   new.data <- create_probs(onemap_example_out, prob = 10^-5)
##'   
##'   
##'@importFrom reshape2 melt
##'@export
create_probs <- function(onemap.obj = NULL, 
                         global_error = NULL, 
                         genotypes_errors = NULL, 
                         genotypes_probs = NULL){
  
  if(class(onemap.obj)[1]!="onemap"){
    stop("onemap.obj should be of class onemap")
  }
  
  if(all(is.null(c(global_error, genotypes_errors, genotypes_probs)))){
    global_error <- 10^-5
  }
  
  crosstype <- class(onemap.obj)[2]
  
  if(!is.null(global_error)){
    probs <- reshape2::melt(t(onemap.obj$geno))
    probs$type <- rep(onemap.obj$segr.type.num, onemap.obj$n.ind)
    
    # Global error according to observated data
    
    if(crosstype == "outcross"){
      prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      # A
      idx <- which(probs$value == 1  & probs$type == 1)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error/3, 3*length(idx)))
      idx <- which(probs$value == 2  & probs$type == 1)
      prob[idx,] <- c(rep(global_error/3, length(idx)), rep(1-global_error, length(idx)), rep(global_error/3, 2*length(idx)))
      idx <- which(probs$value == 3  & probs$type == 1)
      prob[idx,] <- c(rep(global_error/3, 2*length(idx)), rep(1-global_error, length(idx)), rep(global_error/3, length(idx)))
      idx <- which(probs$value == 4  & probs$type == 1)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      
      # B1
      idx <- which(probs$value == 1  & probs$type == 2)
      prob[idx,] <- c(rep(1- global_error, 2*length(idx)), rep(global_error/2, 2*length(idx)))
      idx <- which(probs$value == 2  & probs$type == 2)
      prob[idx,] <- c(rep(global_error/3, 2*length(idx)), rep(1-global_error, length(idx)), rep(global_error/3, length(idx)))
      idx <- which(probs$value == 3  & probs$type == 2)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      
      # B2
      idx <- which(probs$value == 1  & probs$type == 3)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error/2, length(idx)), rep(1- global_error, length(idx)), rep(global_error/2, length(idx)))
      idx <- which(probs$value == 2  & probs$type == 3)
      prob[idx,] <- c(rep(global_error/3, length(idx)), rep(1-global_error, length(idx)), rep(global_error/3, 2*length(idx)))
      idx <- which(probs$value == 3  & probs$type == 3)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      
      # B3.7
      idx <- which(probs$value == 1  & probs$type == 4)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error/3, 3*length(idx)))
      idx <- which(probs$value == 2  & probs$type == 4)
      prob[idx,] <- c(rep(global_error/2, length(idx)), rep(1-global_error, 2*length(idx)), rep(global_error/2, length(idx)))
      idx <- which(probs$value == 3  & probs$type == 4)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      
      # C
      idx <- which(probs$value == 1  & probs$type == 5)
      prob[idx,] <- c(rep(1- global_error, 3*length(idx)), rep(global_error, length(idx)))
      idx <- which(probs$value == 2  & probs$type == 5)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      
      # D1
      idx <- which(probs$value == 1  & probs$type == 6)
      prob[idx,] <- c(rep(1- global_error, 2*length(idx)), rep(global_error/2, 2*length(idx)))
      idx <- which(probs$value == 2  & probs$type == 6)
      prob[idx,] <- c(rep(global_error/2, 2*length(idx)), rep(1-global_error, 2*length(idx)))
      
      # D2
      idx <- which(probs$value == 1  & probs$type == 7)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error/2, length(idx)), rep(1- global_error, length(idx)), rep(global_error/2, length(idx)))
      idx <- which(probs$value == 2  & probs$type == 7)
      prob[idx,] <- c(rep(global_error/2, length(idx)), rep(1-global_error, length(idx)), rep(global_error/2, length(idx)), rep(1-global_error, length(idx)))
      
    } else if(crosstype == "f2"){
      prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      idx <- which(probs$value == 1)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error/3, 3*length(idx)))
      idx <- which(probs$value == 2)
      prob[idx,] <- c(rep(global_error/2, length(idx)), rep(1-global_error, 2*length(idx)), rep(global_error/2, length(idx)))
      idx <- which(probs$value == 3)
      prob[idx,] <- c(rep(global_error/3, 3*length(idx)), rep(1-global_error, length(idx)))
      idx <- which(probs$value == 4)
      prob[idx,] <- c(rep(1- global_error/3, 3*length(idx)), rep(global_error, length(idx)))
      idx <- which(probs$value == 5)
      prob[idx,] <- c(rep(global_error, length(idx)), rep(1-global_error/3, 3*length(idx)))
    } else if(crosstype == "backcross" | crosstype == "riself" | crosstype == "risib"){
      prob <- matrix(NA, nrow=length(probs$value), ncol = 2)
      idx <- which(probs$value == 0)
      prob[idx,] <- 1
      idx <- which(probs$value == 1)
      prob[idx,] <- c(rep(1- global_error, length(idx)), rep(global_error, length(idx)))
      idx <- which(probs$value == 2)
      prob[idx,] <- c(rep(global_error, length(idx)), rep(1-global_error, length(idx)))
      idx <- which(probs$value == 3)
      prob[idx,] <- c(rep(global_error, length(idx)), rep(1-global_error, length(idx)))
    } 
    rownames(prob) <- paste0(probs$Var1, "_", probs$Var2)
    onemap.obj$error <- prob
  }
  
  # Fazer para erros por genótipo e para probabilidades
  
  return(onemap.obj)
}
