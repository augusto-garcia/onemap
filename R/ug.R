#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: ug.R                                                          #
# Contains: ug                                                        #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2007-9, Marcelo Mollinari                             #
#                                                                     #
# First version: 11/29/2009                                           #
# Last update: 01/21/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
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
##' recombination fractions between markers. Next is an adapted excerpt from
##' \cite{Mollinari et al (2009)} describing the \emph{UG} algorithm:
##' 
##' Based on the \eqn{R} (recombination fraction) matrix, the distance
##' between all \eqn{m} loci is calculated by \eqn{d_{ij} = }{d_ij = r_ij +
##' (2/n_ij) Sum_k r_ik r_jk}\eqn{ \hat{r}_{ij} + (\frac{2}{n_{ij}}) \sum_{k}
##' \hat{r}_{ik} }{d_ij = r_ij + (2/n_ij) Sum_k r_ik r_jk}\eqn{
##' \hat{r}_{jk}}{d_ij = r_ij + (2/n_ij) Sum_k r_ik r_jk}, for every \eqn{k},
##' with \eqn{\hat{r}_{ij} > \hat{r}_{ik}, \hat{r}_{ij} > }{r_ij > r_ik, r_ij >
##' r_jk}\eqn{ \hat{r}_{jk}}{r_ij > r_ik, r_ij > r_jk}, and \eqn{n_{ij}}{n_ij}
##' individuals. The value \eqn{T_{ij} = 2 d_{ij} - (\sum_{k \neq i} }{T_ij = 2
##' d_ij - (Sum_{k != i} d_ik + Sum_{k != j} d_jk)}\eqn{ d_{ik} + \sum_{k \neq
##' j} d_{jk})}{T_ij = 2 d_ij - (Sum_{k != i} d_ik + Sum_{k != j} d_jk)} is
##' calculated for every \eqn{i < j}. The terminal end of the map is defined by
##' taking the pair of markers \eqn{(f, g)} that presents the smallest value of
##' \eqn{T}. The pair \eqn{(f, g)} is then denoted locus \eqn{m + 1} and its
##' distance to the remaining markers is determined by \eqn{d_{i m+1} =
##' \frac{1}{2}(d_{if} + }{d_im+1 = (1/2)(d_if + d_ig - d_fg)}\eqn{ d_{ig} -
##' d_{fg})}{d_im+1 = (1/2)(d_if + d_ig - d_fg)} if \eqn{(d_{if} + }{(d_if +
##' d_ig) > d_fg}\eqn{ d_{ig}) > d_{fg}}{(d_if + d_ig) > d_fg}, if not,
##' \eqn{d_{i m+1} = }{d_im+1 = 0}\eqn{ 0}{d_im+1 = 0}. The calculation
##' \eqn{W_{i m+1} = (m - 2) d_{i m+1} - }{W_im+1 = (m-2) d_im+1 - Sum_{k != i}
##' d_ik}\eqn{ \sum_{k \neq i} d_{ik}}{W_im+1 = (m-2) d_im+1 - Sum_{k != i}
##' d_ik} is also performed and the locus that minimizes the value \eqn{W_{i
##' }{W_im+1}\eqn{ m+1}}{W_im+1} (called locus \eqn{h}) is placed on the map.
##' The partial resultant map is \emph{f-g-h} if \eqn{d_{fh} > d_{gh}}{d_fh >
##' d_gh} or \emph{h-f-g} otherwise. Considering \eqn{k = 2}, the partial
##' distance of the map with the remaining markers is updated: \eqn{d_{i m+k} =
##' \mbox{min}(d_{i m+k-1}, d_{ij})}{d_im+k = min(d_im+k-1, d_ij)}. The value
##' \eqn{W_{i m+k} = (m - k - 1) d_{i }{W_im+k = (m-k-1) d_im+k - Sum_{k != i}
##' d_ik}\eqn{ m+k} - \sum_{k \neq i} d_{ik}}{W_im+k = (m-k-1) d_im+k - Sum_{k
##' != i} d_ik} is calculated and the locus that minimizes \eqn{W} is added to
##' the map. The last two steps are repeated, taking \eqn{k = 3, }{k = 3, ...,
##' m-1}\eqn{ \ldots, m - 1}{k = 3, ..., m-1} to obtain the complete map.
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
##' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means taht there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{outcross} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.}
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{make.seq}}, \code{\link[onemap]{map}}
##' @references Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
##' A. A. F. (2009) Evaluation of algorithms used to order markers on genetics
##' maps. \emph{Heredity} 103: 494-502.
##' 
##' Tan, Y. and Fu, Y. (2006) A novel method for estimating linkage maps.
##' \emph{Genetics} 173: 2383-2390.
##' @keywords utilities
##' @examples
##' 
##' \dontrun{
##'   #outcross example
##'   data(example.out)
##'   twopt <- rf.2pts(example.out)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.ug <- ug(LG1)
##' 
##'   #F2 example
##'   data(fake.f2.onemap)
##'   twopt <- rf.2pts(fake.f2.onemap)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.ug <- ug(LG1)
##'   LG1.ug
##' }
##' 
ug<-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5){
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
  n.mark<-ncol(r)
  
  ## For two markers
  if(n.mark==2)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))

  ## UG algorithm
  ## The equation numbers below refer to the ones in the article (Tan and Fu, 2006)
  
  ##eq 1 and eq 2
  calc.T<-function(d,i,j) 2*d[i,j]-(sum(d[i,-i], na.rm = TRUE)+sum(d[j,-j], na.rm = TRUE))
  
  ##Function to pick the position that has the minimum value (considering ties)
  pick.pos.min<-function(x){
    if (length(which(min(x,na.rm=TRUE)==x))>1) return(sample((which(min(x,na.rm=TRUE)==x)),1))
    else return(which(min(x,na.rm=TRUE)==x))
  } 
  lambda<-matrix(1,n.mark,n.mark) # see the last paragraph on page 2385
  
  ##Function to calculate distance between loci i and j
  ##eq.13
  r.to.d<-function(r,i,j,lambda){
    N<-0;q<-0
    for(k in 1:n.mark){
      if(r[i,j]>r[i,k] && r[i,j]>r[j,k]){ 
        q<-q+(r[i,k]*r[j,k])
        N<-N+1
      }
    }
    if(N==0) N<-1
    return(r[i,j]+(2*lambda[i,j]/N*q))
  }
  
  ##Calculating d
  d<-matrix(NA,(2*n.mark-1),(2*n.mark-1))
  for(i in 1:n.mark){
    for(j in i:n.mark){
      d[i,j]<-d[j,i]<-r.to.d(r,i,j,lambda)
    }
  }
  
  ##Calculating T 
  T<-matrix(NA,n.mark,n.mark)
  for(i in 1:n.mark){
    for(j in i:n.mark){
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
  new.locus<-n.mark+1
  
  din1<-function(x,y,d,n.mark){
    d.new.v<-rep(NA, n.mark)
    for(i in 1:n.mark){
      if((d[i,x]+d[i,y])>d[x,y]){
        d.new.v[i]<-.5*(d[i,x]+d[i,y]-d[x,y]) #eq7
      }
      else d.new.v[i]<-0
    }
    return(d.new.v)
  }
  
  d.old1<-d[partial[1],]
  d.old2<-d[partial[2],]
  d[new.locus,c(1:n.mark)]<-d[c(1:n.mark),new.locus]<-d[new.locus,c(1:n.mark)]<-din1(partial[1],partial[2],d,n.mark)
  d[,partial[2]]<-d[partial[2],]<-d[,partial[1]]<- d[partial[1],]<-NA
  H.temp<-c()
  for(i in 1:n.mark){
    H.temp[i]<-(n.mark-2)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE) #eq8
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
  if(n.mark==3)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(partial)],twopt=input.seq$twopt), tol=10E-5))
  
  for (k in 2:(n.mark-2)){
    ##step 2
    new.locus<-n.mark+k
    for(i in 1:(new.locus-1)){
      d[i,new.locus]<-d[new.locus,i]<-min(d[(new.locus-1),i], d[partial[k+1],i])#eq9
    }
    d[partial[k+1],]<-NA
    d[,partial[k+1]]<-NA
    d[new.locus-1,]<-NA
    d[,new.locus-1]<-NA
    
    ##step 3
    H.temp<-c()
    for(i in 1:n.mark){
      H.temp[i]<-(n.mark-k-1)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE)#eq10
    }
    
    H<-rbind(H,H.temp)
    next.mark<-pick.pos.min(H[k,])
    partial<-c(partial, next.mark)
  }
  complete<-partial
  ## end of UG algorithm
  cat("\norder obtained using UG algorithm:\n\n", input.seq$seq.num[avoid.reverse(complete)], "\n\ncalculating multipoint map using tol ", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(complete)],twopt=input.seq$twopt), tol=tol)
  
}
## end of file
  
















  
