#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: ug.R                                                          ##
## Contains: ug                                                        ##
##                                                                     ##
## Written by Marcelo Mollinari                                        ##
## copyright (c) 2007-9, Marcelo Mollinari                             ##
##                                                                     ##
## First version: 11/29/2009                                           ##
## Description modified by Augusto Garcia on 2015/07/25                ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################



##' Unidirectional Growth
##'
##' Implements the marker ordering algorithm \emph{Unidirectional Growth}
##' (\cite{Tan & Fu, 2006}).
##'
##' \emph{Unidirectional Growth} (\emph{UG}) is an algorithm for marker
##' ordering in linkage groups. It is not an exhaustive search method and,
##' therefore, is not computationally intensive. However, it does not guarantee
##' that the best order is always found. The only requirement is a matrix with
##' recombination fractions between markers.
##'
##' After determining the order with \emph{UG}, the final map is constructed
##' using the multipoint approach (function \code{\link[onemap]{map}}).
##'
##' @param input.seq an object of class \code{sequence}.
##' @param LOD minimum LOD-Score threshold used when constructing the pairwise
##' recombination fraction matrix.
##' @param max.rf maximum recombination fraction threshold used as the LOD
##' value above.
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
#' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
#' if \code{TRUE} one of the markers is removed and ug is performed again.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param phase_cores The number of parallel processes to use when estimating
##' the phase of a marker. (Should be no more than 4)
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
#' @param hmm logical defining if the HMM must be applied to estimate multipoint
#' genetic distances
##' @param verbose A logical, if TRUE it output progress status
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
##' \item{data.name}{object of class \code{onemap} with the raw
##' data.} \item{twopt}{object of class \code{rf_2pts} with the
##' 2-point analyses.}
##' 
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{make_seq}}, \code{\link[onemap]{map}}
##' @references Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
##' A. A. F. (2009) Evaluation of algorithms used to order markers on genetics
##' maps. \emph{Heredity} 103: 494-502.
##'
##' Tan, Y. and Fu, Y. (2006) A novel method for estimating linkage maps.
##' \emph{Genetics} 173: 2383-2390.
##' @keywords utilities
##' @examples
##'
##' \donttest{
##'   #outcross example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.ug <- ug(LG1)
##'
##'   #F2 example
##'   data(mapmaker_example_f2)
##'   twopt <- rf_2pts(mapmaker_example_f2)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.ug <- ug(LG1)
##'   LG1.ug
##' }
##'@export
ug<-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5, 
             rm_unlinked = TRUE,
             size = NULL, 
             overlap = NULL, 
             phase_cores = 1, hmm=TRUE, parallelization.type = "PSOCK", verbose=TRUE)
{
  ## checking for correct object
  if(!inherits(input.seq,"sequence"))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  n.mrk <- length(input.seq$seq.num)
  
  ## create reconmbination fraction matrix
  if(inherits(input.seq$twopt,"outcross") || inherits(input.seq$twopt,"f2"))
    r<-get_mat_rf_out(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
  else
    r<-get_mat_rf_in(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
  r[is.na(r)]<-0.5
  diag(r)<-0
  
  ## For two markers
  if(n.mrk==2)
    return(map(make_seq(input.seq$twopt,input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))
  
  ## UG algorithm
  ## The equation numbers below refer to the ones in the article (Tan and Fu, 2006)
  
  ##eq 1 and eq 2
  calc.T<-function(d,i,j) 2*d[i,j]-(sum(d[i,-i], na.rm = TRUE)+sum(d[j,-j], na.rm = TRUE))
  
  ##Function to pick the position that has the minimum value (considering ties)
  pick.pos.min<-function(x){
    if (length(which(min(x,na.rm=TRUE)==x))>1) return(sample((which(min(x,na.rm=TRUE)==x)),1))
    else return(which(min(x,na.rm=TRUE)==x))
  }
  lambda<-matrix(1,n.mrk,n.mrk) # see the last paragraph on page 2385
  
  ##Function to calculate distance between loci i and j
  ##eq.13
  r.to.d<-function(r,i,j,lambda){
    N<-0;q<-0
    for(k in 1:n.mrk){
      if(r[i,j]>r[i,k] && r[i,j]>r[j,k]){
        q<-q+(r[i,k]*r[j,k])
        N<-N+1
      }
    }
    if(N==0) N<-1
    return(r[i,j]+(2*lambda[i,j]/N*q))
  }
  
  ##Calculating d
  d<-matrix(NA,(2*n.mrk-1),(2*n.mrk-1))
  for(i in 1:n.mrk){
    for(j in i:n.mrk){
      d[i,j]<-d[j,i]<-r.to.d(r,i,j,lambda)
    }
  }
  
  ##Calculating T
  T<-matrix(NA,n.mrk,n.mrk)
  for(i in 1:n.mrk){
    for(j in i:n.mrk){
      T[i,j]<-calc.T(d,j,i)
    }
    diag(T)<-NA
  }
  
  ##verification of eq 3
  ##linkage<-c(T[1,2],T[19,20],T[20,21])
  ##min(linkage)== min(T, na.rm=T)
  
  
  ##UG algorithm itself
  
  ##step 1
  ##Finding the smallest T-value
  partial<-which(min(T,na.rm=TRUE)==T, arr.ind = TRUE)
  if(length(partial)>2) partial<-partial[sample(c(1:nrow(partial)),1),]
  new.locus<-n.mrk+1
  
  din1<-function(x,y,d,n.mrk){
    d.new.v<-rep(NA, n.mrk)
    for(i in 1:n.mrk){
      if((d[i,x]+d[i,y])>d[x,y]){
        d.new.v[i]<-.5*(d[i,x]+d[i,y]-d[x,y]) #eq7
      }
      else d.new.v[i]<-0
    }
    return(d.new.v)
  }
  
  d.old1<-d[partial[1],]
  d.old2<-d[partial[2],]
  d[new.locus,c(1:n.mrk)]<-d[c(1:n.mrk),new.locus]<-d[new.locus,c(1:n.mrk)]<-din1(partial[1],partial[2],d,n.mrk)
  d[,partial[2]]<-d[partial[2],]<-d[,partial[1]]<- d[partial[1],]<-NA
  H.temp<-c()
  for(i in 1:n.mrk){
    H.temp[i]<-(n.mrk-2)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE) #eq8
  }
  
  H<-H.temp
  next.mark<-pick.pos.min(H)
  side<-c(d.old1[next.mark],d.old2[next.mark])
  if(pick.pos.min(side)==1){
    a<-partial[1]
    partial[1]<-partial[2]
    partial[2]<-a
  }
  partial<-c(partial, next.mark)
  
  ## If there are three markers, do not go to the second step
  if(n.mrk==3)
    return(map(make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(partial)],twopt=input.seq$twopt), 
               tol=10E-5, parallelization.type= parallelization.type))
  
  for (k in 2:(n.mrk-2)){
    ##step 2
    new.locus<-n.mrk+k
    for(i in 1:(new.locus-1)){
      d[i,new.locus]<-d[new.locus,i]<-min(d[(new.locus-1),i], d[partial[k+1],i])#eq9
    }
    d[partial[k+1],]<-NA
    d[,partial[k+1]]<-NA
    d[new.locus-1,]<-NA
    d[,new.locus-1]<-NA
    
    ##step 3
    H.temp<-c()
    for(i in 1:n.mrk){
      H.temp[i]<-(n.mrk-k-1)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE)#eq10
    }
    
    H<-rbind(H,H.temp)
    next.mark<-pick.pos.min(H[k,])
    partial<-c(partial, next.mark)
  }
  complete<-partial
  ## end of UG algorithm
  if(hmm){
    if(verbose) cat("\norder obtained using UG algorithm:\n\n", input.seq$seq.num[avoid_reverse(complete)], "\n\ncalculating multipoint map using tol ", tol, ".\n\n")
    if(phase_cores == 1 | inherits(input.seq$data.name, c("backcross", "riself", "risib"))){
      ug_map <- map(make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(complete)],
                             twopt=input.seq$twopt), tol=tol, rm_unlinked = rm_unlinked, 
                    parallelization.type= parallelization.type)
    } else{
      if(is.null(size) | is.null(overlap)){
        stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you must also define `size` and `overlap` arguments.")
      } else {
        ug_map <- map_overlapping_batches(make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(complete)],
                                                   twopt=input.seq$twopt), 
                                          tol=tol,
                                          size = size, overlap = overlap, 
                                          phase_cores = phase_cores,
                                          rm_unlinked = rm_unlinked, parallelization.type= parallelization.type)
      }
    }
    
    if(!is.list(ug_map)) {
      new.seq <- make_seq(input.seq$twopt, ug_map)
      ug_map <- ug(input.seq = new.seq, 
                   LOD=LOD, 
                   max.rf=max.rf, tol=tol, 
                   rm_unlinked= rm_unlinked,
                   size = size, 
                   overlap = overlap, 
                   phase_cores = phase_cores, parallelization.type= parallelization.type)
    }
    return(ug_map)
  } else {
    ug.seq <- make_seq(input.seq$twopt,input.seq$seq.num[avoid_reverse(complete)],
                       twopt=input.seq$twopt)
    return(ug.seq) 
  }
}
## end of file


















