#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## Contains: parmap and map_avoid_unlinked                             ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
## copyright (c) 2019, Cristiane Taniguti                              ##
##                                                                     ##
## First version: 11/11/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

globalVariables(c("var"))

##' A parallelized version of map function
##' 
##' Based on the strategy propoused by Schiffthaler et al. the markers of a
##' linkage group with pre-defined order are divided in user defined batches 
##' with overlap markers between them. Each batch runs in a different core
##' which increase the speed of the genetic distances estimation. The distances 
##' and phases of the previous batch are considered to joint the batches. Also, 
##' it is possible to store the differences between the genetic distances 
##' between the overlaped markers as a measure of the process quality. In 
##' other words, if these differences are to too high, you can considered that
##' divide the group in batch did not compromise the distance estimation.
##'  
##' @param input.seq an object of class \code{sequence}.
##' @param cores an integer defining the number of cores to be used and also the
##' number of batches generated
##' @param overlap number of markers overlapping between batches
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @param avoid_link_errors logical. If TRUE markers which do not reach the linkage
##' criteria are removed of the sequence and the distances are automatically reestimated.
##' If FALSE an error stops the algorithm it find these markers.
##' @param export_diff If TRUE also returns (in the first level of the list) the differences
##' in genetic distances between overlaped markers
##' @param verbose A logical, if TRUE its output progress status
##' information.
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
##' \item{data.name}{name of the object of class \code{onemap} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
##' 2-point analyses.}
##' 
##' @author Cristiane Taniguti \email{chtaniguti@@tamu.edu} 
##' @seealso \code{\link[onemap]{map}}
##' 
##' @references Schiffthaler, B., Bernhardsson, C., Ingvarsson, P. K., & Street, 
##' N. R. (2017). BatchMap: A parallel implementation of the OneMap R package 
##' for fast computation of F1 linkage maps in outcrossing species. PLoS ONE, 
##' 12(12), 1–12. https://doi.org/10.1371/journal.pone.0189256
##' 
##' @import parallel
##'  
##' @export
parmap <- function(input.seq=NULL, 
                   cores=3, 
                   overlap=4, 
                   tol=10E-5, 
                   avoid_link_errors = TRUE,
                   export_diff = FALSE, verbose=TRUE){
  
  twopts <- input.seq$twopt
  
  interv <- length(input.seq$seq.num)/cores
  
  seqs <- 1:length(input.seq$seq.num)
  
  list_idx <- list(1:cores)
  
  init <- 1
  end <- interv
  for(i in 1:cores){
    if(i != cores){
      list_idx[[i]] <- init:(end+(overlap-1))
      init <- end
      end <- end + interv
    } else{
      list_idx[[i]] <- init:end
    }
  }
  
  list_seq <- lapply(list_idx, function(x) make_seq(twopts, input.seq$seq.num[x]))
  
  
  clust <- makeCluster(cores)
  clusterExport(clust, c("map_avoid_unlinked"))
  if(avoid_link_errors){
    new.maps <- parLapply(clust, list_seq, function(x) map_avoid_unlinked(x, tol=tol))      
  } else {
    new.maps <- parLapply(clust, list_seq, function(x) map(x, tol=tol))
  }
  
  stopCluster(clust)
  
  joint.map <- new.maps[[1]]
  new.seq.num <- new.seq.rf <- new.seq.phases <- vector()
  diff1 <- vector()
  for(i in 1:(length(new.maps)-1)){
    idx.end <- (length(new.maps[[i]]$seq.num) - overlap+1):(length(new.maps[[i]]$seq.num))
    new.seq.num <- c(new.seq.num, new.maps[[i]]$seq.num[-idx.end])
    new.seq.rf <- c(new.seq.rf, new.maps[[i]]$seq.rf[-idx.end[-length(idx.end)]])
    new.seq.phases <- c(new.seq.phases, new.maps[[i]]$seq.phases[-idx.end[-length(idx.end)]])
    
    if(i == (length(new.maps)-1)){
      new.seq.num <- c(new.seq.num, new.maps[[i+1]]$seq.num)
      new.seq.rf <- c(new.seq.rf, new.maps[[i+1]]$seq.rf)
      new.seq.phases <- c(new.seq.phases, new.maps[[i+1]]$seq.phases)
    }
    
    end <- new.maps[[i]]$seq.rf[idx.end[-length(idx.end)]]
    init <- new.maps[[i+1]]$seq.rf[1:overlap-1]
    diff1 <- c(diff1,end - init)
  }
  
  if(verbose) cat("The overlap markers have mean ", mean(diff1), " of  recombination fraction differences, and variance of ", var(diff1), "\n")
  
  joint.map$seq.num <- new.seq.num
  joint.map$seq.phases <- new.seq.phases
  joint.map$seq.rf <- new.seq.rf
  
  if(export_diff){
    return(list(diff1,joint.map))
  } else {
    return(joint.map)
  }
}



