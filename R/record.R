#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: record.R                                                      #
# Contains: record                                                    #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2007-9, Marcelo Mollinari                             #
#                                                                     #
# First version: 11/29/2009                                           #
# Last update: 01/21/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
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
##' a matrix with recombination fractions between markers. Next is an adapted
##' excerpt from \cite{Mollinari et al (2009)} describing the \emph{RECORD}
##' algorithm:
##' 
##' Based on the expected number of recombination events, an \eqn{S}
##' matrix is constructed, \eqn{S = [S_{M_{i}M_{j}}]_{m \times }{S =
##' [S_MiMj]_(m x m)}\eqn{ m}}{S = [S_MiMj]_(m x m)} (for \eqn{M_{i} =
##' M_{j}}{M_i = M_j}, \eqn{S_{M_{i}M_{j}} = 0}{S_MiMj = 0}), where \eqn{m} is
##' the number of markers. The procedure to obtain \eqn{S} is based on the
##' expected number of crossovers between marker pairs, conditioned by the
##' observation of the markers' phenotype. The optimization criterion
##' \emph{COUNT} for a sequence of \eqn{m} markers may be calculated by
##' \eqn{COUNT = \sum_{i=1}^{m-1} S_{M_{i}M_{i+1}}}{COUNT = Sum_{i=1}^{m-1}
##' S_MiMi+1}, where smaller \emph{COUNT} values correspond to better orders.
##' Map building is carried out by randomly taking two markers and positioning
##' a third one at the beginning, at the end and between them. The marker is
##' fixed at the position that gives a smaller value of \emph{COUNT}.
##' Similarly, the remaining markers are positioned at pre-established orders
##' until completion of the map. Subsequently, a search for smaller values of
##' \emph{COUNT} is performed, inverting the position on the map of
##' subsequences of size \eqn{m' = 2, \ldots, 20}{m' = 2, ..., 20}. If the map
##' resulting from the inverted positions presents a \emph{COUNT} value smaller
##' than the previous one, it is kept. The procedure is repeated \code{times}
##' times and the sequence presenting the smallest \emph{COUNT} value is
##' chosen.
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
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{outcross} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.}
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{make.seq}} and \code{\link[onemap]{map}}
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
##'   data(example.out)
##'   twopt <- rf.2pts(example.out)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rec <- record(LG1)
##' 
##'   ##F2 example
##'   data(fake.f2.onemap)
##'   twopt <- rf.2pts(fake.f2.onemap)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rec <- record(LG1)
##'   LG1.rec
##' }
##' 
record<-function(input.seq, times=10, LOD=0, max.rf=0.5, tol=10E-5){
  ## checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
  markers <- length(input.seq$seq.num)

  ## create recombination fraction matrix 
  r <- matrix(NA,markers,markers)
  for(i in 1:(markers-1)) {
    for(j in (i+1):markers) {
      big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[j])
      small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[j])
      temp <- get(input.seq$twopt)$analysis[acum(big-2)+small,,]
      ## check if any assignment meets the criteria
      relevant <- which(temp[,2] > (max(temp[,2])-0.005)) # maximum LOD scores
      phases <- relevant[which((temp[relevant,1] <= max.rf & temp[relevant,2] >= LOD))]
      if(length(phases) == 0) r[i,j] <- r[j,i] <- 0.5
      else r[i,j] <- r[j,i] <- temp[phases[1],1]
    }
  }
  diag(r)<-0
  
  ##RECORD algorithm
  n.mark<-ncol(r)
  X<-r*get(input.seq$data.name, pos=1)$n.ind ## Obtaining X multiplying the MLE of the recombination
                                      ## fraction by the number of individuals
  
  COUNT<-function(X, sequence){ ## See eq. 1 on the paper (Van Os et al., 2005)
    return(sum(diag(X[sequence[-length(sequence)],sequence[-1]]),na.rm=TRUE))
  }

  ## For two markers
  if(n.mark==2)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))
  
  ## For three markers (calculation of 3 possible orders: 3!/2)
  else if(n.mark==3) {
    all.perm<-perm.pars(1:3)
    m.old<-Inf
    for(k in 1:nrow(all.perm)){
      m.new<-COUNT(X, all.perm[k,])
      if(m.new < m.old)
        result.new <- all.perm[k,]; m.old<-m.new
    }
  }
  
  ## For more than three markers (RECORD algorithm itself)
  else{
    marks<-c(1:n.mark)
    result.new<-sample(marks)## randomize markers
    for(k in 1:times){ ## loop for replicates the RECORD procedure
      result<-sample(marks, 2)
      for(i in 2:(n.mark-2)){
        next.mark<-sample(c(1:n.mark)[-result],1)
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
      next.mark<-c(1:n.mark)[-result]
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
      for(i in 1:(n.mark-1)){
        for(j in 1:(n.mark-i)){
          y<-result[j:(j+i)]
          if(j==1){
            if(i!=n.mark-1){
              perm.order<-c(rev(result[1:(1+i)]),result[(i+2):length(result)])
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp) {
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
          else{
            if(j!=(n.mark-i)){
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
 cat("\norder obtained using RECORD algorithm:\n\n", input.seq$seq.num[avoid.reverse(result.new)], "\n\ncalculating multipoint map using tol", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(result.new)],twopt=input.seq$twopt), tol=tol)
}

##end of file
  
















  
