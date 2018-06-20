#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: seriation.R                                                   ##
## Contains: seriation, ser_ord, Cindex                                ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido                        ##
## copyright (c) 2007-9, Gabriel R A Margarido                         ##
##                                                                     ##
## First version: 11/29/2009                                           ##
## Last update: 12/2015                                                ##
## Description was update by Augusto Garcia on 2015/07/25              ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################



##' Seriation
##'
##' Implements the marker ordering algorithm \emph{Seriation} (\cite{Buetow &
##' Chakravarti, 1987}).
##'
##' \emph{Seriation} is an algorithm for marker ordering in linkage groups. It
##' is not an exhaustive search method and, therefore, is not computationally
##' intensive. However, it does not guarantee that the best order is always
##' found. The only requirement is a matrix with recombination fractions
##' between markers.
##'
##' NOTE: When there are to many pairs of markers with the same value in the
##' recombination fraction matrix, it can result in ties during the ordination
##' process and the \emph{Seriation} algorithm may not work properly. This is
##' particularly relevant for outcrossing populations with mixture of markers
##' of type \code{D1} and \code{D2}. When this occurs, the function shows the
##' following error message: \code{There are too many ties in the ordination
##' process - please, consider using another ordering algorithm}.
##'
##' After determining the order with \emph{Seriation}, the final map is
##' constructed using the multipoint approach (function
##' \code{\link[onemap]{map}}).
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
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{onemap} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
##' 2-point analyses.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{make_seq}}, \code{\link[onemap]{map}}
##' @references Buetow, K. H. and Chakravarti, A. (1987) Multipoint gene
##' mapping using seriation. I. General methods. \emph{American Journal of
##' Human Genetics} 41: 180-188.
##'
##' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A. A. F.
##' (2009) Evaluation of algorithms used to order markers on genetics maps.
##' \emph{Heredity} 103: 494-502.
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'   ##outcross example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG3 <- make_seq(groups,3)
##'   LG3.ser <- seriation(LG3)
##'
##'   ##F2 example
##'   data(onemap_example_f2)
##'   twopt <- rf_2pts(onemap_example_f2)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.ser <- seriation(LG1)
##'   LG1.ser
##' }
##'@export
seriation<-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5)
{
    ## checking for correct object
    if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
    n.mrk <- length(input.seq$seq.num)

    ## create reconmbination fraction matrix

    if(class(get(input.seq$twopt))[2]=="outcross")
        r<-get_mat_rf_out(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
    else
        r<-get_mat_rf_in(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
    r[is.na(r)]<-0.5
    diag(r)<-0

    ## SERIATION algorithm
    n.mrk<-ncol(r)
    orders <- array(0,dim=c(n.mrk,n.mrk))
    CI <- numeric(n.mrk)
    for (i in 1:n.mrk) {
        orders[i,] <- ser_ord(r,i)
        CI[i] <- Cindex(orders[i,],r)
    }
    best <- which(CI==CI[which.min(CI)])
    if (length(best) == 0) stop ("Cannot find any order")
    else {
        if (length(best) != 1) best <- best[sample(length(best),1)]
        complete<-orders[best,]
    }

    ## end of SERIATION algorithm
    cat("\norder obtained using SERIATION algorithm:\n\n", input.seq$seq.num[complete], "\n\ncalculating multipoint map using tol = ", tol, ".\n\n")
    map(make_seq(get(input.seq$twopt),input.seq$seq.num[complete],twopt=input.seq$twopt), tol=tol)
}

##Provides an order given the recombination
##fraction matrix and the starting marker.
ser_ord <- function(r,i) {
    n.mrk <- ncol(r)
    x <- 1:n.mrk
    unres1 <- 0
    unres2 <- 0
    esq <- numeric(0)
    dir <- numeric(0)
    order <- i
    x[i] <- NaN
    j <- which(r[i,][x]==r[i,which.min(r[i,][x])])
    if (length(j) > 1) j <- j[sample(length(j),1)]
    order <- c(order,j)
    x[j] <- NaN
    while (length(order) != n.mrk || unres1 || unres2[1]) {
        if (unres1) {
            if ((length(order) + length(esq) + length(dir) + 1)==n.mrk) {
                duv <- sample(2,1)
                if (duv==1) order <- c(esq,unres1,order,dir) else order <- c(esq,order,unres1,dir)
                unres1 <- 0; dir <- numeric(0); esq <- numeric(0)
            }
      else {
          pri <- if (length(esq)==0) order[1] else esq[1]
          ult <- if (length(dir)==0) order[length(order)] else dir[length(dir)]
          e <- which(r[i,][x]==r[i,which.min(r[i,][x])])
          if (length(e) > 1) e <- e[sample(length(e),1)]
          rand <- 0
          if (r[pri,e] == r[ult,e]) rand <- sample(2,1)
          if (r[pri,e] < r[ult,e] || rand==1) {
              if (r[e,unres1] < r[ult,unres1]) { order <- c(e,esq,unres1,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else if (r[e,unres1] > r[ult,unres1]) { order <- c(e,esq,order,unres1,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else { esq <- c(e,esq); x[e] <- NaN }
          }
        else if (r[pri,e] > r[ult,e] || rand==2) {
            if (r[e,unres1] < r[pri,unres1]) { order <- c(esq,order,unres1,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else if (r[e,unres1] > r[pri,unres1]) { order <- c(esq,unres1,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else { dir <- c(dir,e); x[e] <- NaN }
        }
      }
        }
    else if (unres2[1]) {
        if ((length(order) + length(esq) + length(dir) + 2)==n.mrk) {
            if (unres2[1]==1) {
                duv <- sample(2,1)
                if (duv==1) order <- c(esq,unres2[2],unres2[3],order,dir) else order <- c(esq,unres2[3],unres2[2],order,dir)
                unres2 <- 0; dir <- numeric(0); esq <- numeric(0)
            }
        else if (unres2[1]==2) {
            duv <- sample(2,1)
            if (duv==1) order <- c(esq,order,unres2[2],unres2[3],dir) else order <- c(esq,order,unres2[3],unres2[2],dir)
            unres2 <- 0; dir <- numeric(0); esq <- numeric(0)
        }
        }
      else {
          pri <- if (length(esq)==0) order[1] else esq[1]
          ult <- if (length(dir)==0) order[length(order)] else dir[length(dir)]
          e <- which(r[i,][x]==r[i,which.min(r[i,][x])])
          if (length(e) > 1) e <- e[sample(length(e),1)]
          rand <- 0
          if (r[pri,e] == r[ult,e]) rand <- sample(2,1)
          m1 <- unres2[2]; m2 <- unres2[3]
          if (r[pri,e] < r[ult,e] || rand==1) {
              if (unres2[1]==1) {
                  if (r[e,m1] < r[e,m2]) { order <- c(e,esq,m1,m2,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(e,esq,m2,m1,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { esq <- c(e,esq); x[e] <- NaN }
              }
          else if (unres2[1]==2) {
              if (r[e,m1] < r[e,m2]) { order <- c(e,esq,order,m1,m2,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(e,esq,order,m2,m1,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { esq <- c(e,esq); x[e] <- NaN }
          }
          }
        else if (r[pri,e] > r[ult,e] || rand==2) {
            if (unres2[1]==1) {
                if (r[e,m1] < r[e,m2]) { order <- c(esq,m2,m1,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(esq,m1,m2,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { dir <- c(dir,e); x[e] <- NaN }
            }
          else if (unres2[1]==2) {
              if (r[e,m1] < r[e,m2]) { order <- c(esq,order,m2,m1,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(esq,order,m1,m2,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { dir <- c(dir,e); x[e] <- NaN }
          }
        }
      }
    }
    else {
        if (length(which(r[i,][x]==r[i,which.min(r[i,][x])]))==1) {
            j <- which.min(r[i,][x])
            if (r[order[1],j] < r[order[length(order)],j]) { order <- c(j,order); x[j] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j]) { order <- c(order,j); x[j] <- NaN }
        else {
            light <- 1
            k <- 2
            l <- length(order)-1
            while (k < l && light) {
                if (r[order[k],j] < r[order[l],j]) { order <- c(j,order); x[j] <- NaN; light <- 0 }
            else if (r[order[k],j] > r[order[l],j]) { order <- c(order,j); x[j] <- NaN; light <- 0 }
            else { k <- k+1; l <- l-1 }
            }
            if (light) { unres1 <- j; x[j] <- NaN }
        }
        }
      else if (length(which(r[i,][x]==r[i,which.min(r[i,][x])]))==2) {
          j <- which(r[i,][x]==r[i,which.min(r[i,][x])])[1]
          k <- which(r[i,][x]==r[i,which.min(r[i,][x])])[2]
          if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) { order <- c(j,order,k); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) { order <- c(k,order,j); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) {
            if (r[order[1],j] < r[order[1],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[1],j] > r[order[1],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN }
          else {
              l <- 2
              light <- 1
              while (l <= length(order) && light) {
                  if (r[order[l],j] < r[order[l],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l+1
              }
              if (light) { unres2 <- c(1,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) {
            if (r[order[length(order)],j] < r[order[length(order)],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[length(order)],j] > r[order[length(order)],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN }
          else {
              l <- length(order)-1
              light <- 1
              while (l >= 1 && light) {
                  if (r[order[l],j] < r[order[l],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l-1
              }
              if (light) { unres2 <- c(2,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else stop("There are too many ties in the ordering process - please, consider using another ordering algorithm.")
      }
      else {
          temp <- sample(length(which(r[i,][x]==r[i,which.min(r[i,][x])])),2)
          m1 <- temp[1]; m2 <- temp[2]
          j <- which(r[i,][x]==r[i,which.min(r[i,][x])])[m1]
          k <- which(r[i,][x]==r[i,which.min(r[i,][x])])[m2]
          if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) { order <- c(j,order,k); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) { order <- c(k,order,j); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) {
            if (r[order[1],j] < r[order[1],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[1],j] > r[order[1],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN }
          else {
              l <- 2
              light <- 1
              while (l <= length(order) && light) {
                  if (r[order[l],j] < r[order[l],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l+1
              }
              if (light) { unres2 <- c(1,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) {
            if (r[order[length(order)],j] < r[order[length(order)],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[length(order)],j] > r[order[length(order)],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN }
          else {
              l <- length(order)-1
              light <- 1
              while (l >= 1 && light) {
                  if (r[order[l],j] < r[order[l],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l-1
              }
              if (light) { unres2 <- c(2,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else stop("There are too many ties in the ordination process - please, consider using another ordering algorithm.")
      }
    }
    }
    return(avoid_reverse(order))
}

##Continuity Index
Cindex <- function (order,r) {
    n.mrk <- dim(r)[1]
    CI <- 0
    for (i in 1:(n.mrk-1))
        for (j in (i+1):n.mrk)
            CI <- CI + r[order[i],order[j]]/((j-i)^2)
    return (CI)
}

## end of file
