#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## Contains: parmap and avoid_unlinked                                 ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
## copyright (c) 2019, Cristiane Taniguti                              ##
##                                                                     ##
## First version: 11/11/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' A parallelized version of map function
##' 
##'  Based on the strategy propoused by Schiffthaler et al. the markers of a
##'  linkage group with pre-defined order are divided in user defined batches 
##'  with overlap markers between them. Each batch runs in a different core
##'  which increase the speed of the genetic distances estimation. The distances 
##'  and phases of the previous batch are considered to joint the batches. Also, 
##'  it is possible to store the differences between the genetic distances 
##'  between the overlaped markers as a measure of the process quality. In 
##'  other words, if these differences are to too high, you can considered that
##'  divide the group in batch did not compromise the distance estimation.
##'  
##'  @param input.seq an object of class \code{sequence}.
##'  @param cores an integer defining the number of cores to be used and also the
##'  number of batches generated
##'  @param overlap number of markers overlapping between batches
##'  @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##'  @param avoid_link_errors logical. If TRUE markers which do not reach the linkage
##'  criteria are removed of the sequence and the distances are automatically reestimated.
##'   If FALSE an error stops the algorithm it find these markers.
##'   
##'  @return
##'  @author Cristiane Taniguti \email{chtaniguti@@usp.br} 
##'  @seealso \code{\link[onemap]{map}}
##'  @references Schiffthaler, B., Bernhardsson, C., Ingvarsson, P. K., & Street, 
##'  N. R. (2017). BatchMap: A parallel implementation of the OneMap R package 
##'  for fast computation of F1 linkage maps in outcrossing species. PLoS ONE, 
##'  12(12), 1–12. https://doi.org/10.1371/journal.pone.0189256
##'     

parmap <- function(input.seq=NULL, 
                   cores=3, 
                   overlap=4, 
                   tol=10E-5, 
                   avoid_link_errors = TRUE){
  
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
  clusterExport(clust, c("avoid_unlinked"))
  if(avoid_link_errors){
    new.maps <- parLapply(clust, list_seq, function(x) avoid_unlinked(x, tol=tol))      
  } else {
    new.maps <- parLapply(clust, list_seq, function(x) onemap::map(x, tol=tol))
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
  
  cat("The overlap markers have mean ", mean(diff1), " of  recombination fraction diff1erences, and variance of ", var(diff1), "\n")
  
  joint.map$seq.num <- new.seq.num
  joint.map$seq.phases <- new.seq.phases
  joint.map$seq.rf <- new.seq.rf
  
  return(list(diff1,joint.map))
}

avoid_unlinked <- function(input.seq, tol){
  map_df <- onemap::map(input.seq, mds.seq = T)
  while(class(map_df) == "integer"){
    seq_true <- onemap::make_seq(input.seq$twopt, map_df)
    map_df <- onemap::map(input.seq = seq_true, mds.seq = T)
  }
  return(map_df)
}

