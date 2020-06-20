#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: record.R                                                      ##
## Contains: record                                                    ##
##                                                                     ##
## Written by Marcelo Mollinari                                        ##
## copyright (c) 2007-9, Marcelo Mollinari                             ##
##                                                                     ##
## First version: 11/29/2009                                           ##
## Last update: 12/2015                                                ##
## The detailed description of the algorithm was removed by Augusto    ##
## Garcia, on 2015/07/25, since it was not compiling in new R versions ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################



##' Recombination Counting and Ordering
##'
##' Implements the marker ordering algorithm \emph{Recombination Counting and
##' Ordering} (\cite{Van Os et al., 2005}).
##'
##' \emph{Recombination Counting and Ordering} (\emph{RECORD}) is an algorithm
##' for marker ordering in linkage groups. It is not an exhaustive search
##' method and, therefore, is not computationally intensive. However, it does
##' not guarantee that the best order is always found. The only requirement is
##' a matrix with recombination fractions between markers.
##'
##' After determining the order with \emph{RECORD}, the final map is
##' constructed using the multipoint approach (function
##' \code{\link[onemap]{map}}).
##'
##' @param input.seq an object of class \code{sequence}.
##' @param times integer. Number of replicates of the RECORD procedure.
##' @param LOD minimum LOD-Score threshold used when constructing the pairwise
##' recombination fraction matrix.
##' @param max.rf maximum recombination fraction threshold used as the LOD
##' value above.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param phase_cores The number of parallel processes to use when estimating
##' the phase of a marker. (Should be no more than 4)
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
#' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
#' if \code{TRUE} one of the markers is removed and record is performed again.
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
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{make_seq}} and \code{\link[onemap]{map}}
##' @references Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
##' A. A. F. (2009) Evaluation of algorithms used to order markers on genetics
##' maps. \emph{Heredity} 103: 494-502.
##'
##' Van Os, H., Stam, P., Visser, R.G.F. and Van Eck, H.J. (2005) RECORD: a
##' novel method for ordering loci on a genetic linkage map. \emph{Theoretical
##' and Applied Genetics} 112: 30-40.
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'   ##outcross example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.rec <- record(LG1)
##'
##'   ##F2 example
##'   data(onemap_example_f2)
##'   twopt <- rf_2pts(onemap_example_f2)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.rec <- record(LG1)
##'   LG1.rec
##' }
##'@export
record<-function(input.seq, times=10, LOD=0, max.rf=0.5, tol=10E-5, 
                 rm_unlinked = TRUE,
                 size = NULL, 
                 overlap = NULL, 
                 phase_cores = 1){
  ## checking for correct object
  if(!is(input.seq,"sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
  n.mrk <- length(input.seq$seq.num)
  
  ## create reconmbination fraction matrix
  
  if(is(input.seq$twopt,"outcross") || is(input.seq$twopt,"f2"))
    r<-get_mat_rf_out(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
  else
    r<-get_mat_rf_in(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
  r[is.na(r)]<-0.5
  diag(r)<-0
  
  ##RECORD algorithm
  X<-r*input.seq$data.name$n.ind ## Obtaining X multiplying the MLE of the recombination
  ## fraction by the number of individuals
  
  COUNT<-function(X, sequence){ ## See eq. 1 on the paper (Van Os et al., 2005)
    return(sum(diag(X[sequence[-length(sequence)],sequence[-1]]),na.rm=TRUE))
  }
  ## For two markers
  if(n.mrk==2)
    return(map(make_seq(input.seq$twopt,input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))
  
  ## For three markers (calculation of 3 possible orders: 3!/2)
  else if(n.mrk==3) {
    all.perm<-perm_pars(1:3)
    m.old<-Inf
    for(k in 1:nrow(all.perm)){
      m.new<-COUNT(X, all.perm[k,])
      if(m.new < m.old)
        result.new <- all.perm[k,]; m.old<-m.new
    }
  }
  
  ## For more than three markers (RECORD algorithm itself)
  else{
    marks<-c(1:n.mrk)
    result.new<-sample(marks)## randomize markers
    for(k in 1:times){ ## loop for replicates the RECORD procedure
      result<-sample(marks, 2)
      for(i in 2:(n.mrk-2)){
        next.mark<-sample(c(1:n.mrk)[-result],1)
        partial<-rep(NA,2*length(result)+1)
        partial[seq(2,(length(partial)-1), by = 2)]<-result
        temp<-10e1000
        for(j in seq(1,(length(partial)), by = 2)){
          partial.temp<-partial
          partial.temp[j]<-next.mark
          temp.new<-COUNT(X, partial.temp[is.na(partial.temp)==FALSE])
          if(temp.new<temp) {
            result<-partial.temp[is.na(partial.temp)==FALSE]
            temp<-temp.new
          }
        }
      }
      next.mark<-c(1:n.mrk)[-result]
      partial<-rep(NA,2*length(result)+1)
      partial[seq(2,(length(partial)-1), by = 2)]<-result
      temp<-10e1000
      for(j in seq(1,(length(partial)), by = 2)){
        partial.temp<-partial
        partial.temp[j]<-next.mark
        temp.new<-COUNT(X, partial.temp[ is.na(partial.temp)==FALSE])
        if(temp.new<temp) {
          result<-partial.temp[ is.na(partial.temp)==FALSE]
          temp<-temp.new
        }
      }
      for(i in 1:(n.mrk-1)){
        for(j in 1:(n.mrk-i)){
          y<-result[j:(j+i)]
          if(j==1){
            if(i!=n.mrk-1){
              perm.order<-c(rev(result[1:(1+i)]),result[(i+2):length(result)])
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp) {
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
          else{
            if(j!=(n.mrk-i)){
              perm.order<-c(result[1:(j-1)],rev(result[j:(j+i)]),result[(j+i+1):length(result)])
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }
            else{
              perm.order<-c(result[1:(j-1)],rev(result[j:length(result)]))
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
        }
      }
      if(COUNT(X,result.new) > COUNT(X,result)){
        result.new<-result
        ##print(COUNT(X,result.new))
      }
    }
  }
  ## end of RECORD algorithm
  cat("\norder obtained using RECORD algorithm:\n\n", input.seq$seq.num[avoid_reverse(result.new)], "\n\ncalculating multipoint map using tol", tol, ".\n\n")
  
  if(phase_cores == 1){
    record_map <- map(make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(result.new)],
                               twopt=input.seq$twopt), 
                      tol=tol,
                      rm_unlinked = rm_unlinked)
    
  } else{
    if(is.null(size) | is.null(overlap)){
      stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you must also define `size` and `overlap` arguments.")
    } else {
      record_map <- map_overlapping_batches(make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(result.new)],
                                                     twopt=input.seq$twopt), 
                                            tol=tol,
                                            size = size, overlap = overlap, 
                                            phase_cores = phase_cores,
                                            rm_unlinked = rm_unlinked)
    }
  }
  
  if(!is.list(record_map)) {
    new.seq <- make_seq(input.seq$twopt, record_map)
    record_map <- record(input.seq = new.seq, 
                         LOD=LOD, 
                         max.rf=max.rf, tol=tol, 
                         rm_unlinked= rm_unlinked,
                         size = size, 
                         overlap = overlap, 
                         phase_cores = phase_cores)
  }
  return(record_map)
}

##end of file


















