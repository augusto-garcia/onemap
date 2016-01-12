#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: rcd.R                                                         ##
## Contains: rcd                                                       ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido                        ##
## copyright (c) 2007-9, Gabriel R A Margarido                         ##
##                                                                     ##
## First version: 11/07/2007                                           ##
## Last update: 12/2015                                                ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################



##' Rapid Chain Delineation
##' 
##' Implements the marker ordering algorithm \emph{Rapid Chain Delineation}
##' (\cite{Doerge, 1996}).
##' 
##' \emph{Rapid Chain Delineation} (\emph{RCD}) is an algorithm for marker
##' ordering in linkage groups. It is not an exhaustive search method and,
##' therefore, is not computationally intensive. However, it does not guarantee
##' that the best order is always found. The only requirement is a matrix with
##' recombination fractions between markers.  Next is an excerpt from QTL
##' Cartographer Version 1.17 Manual describing the \emph{RCD} algorithm
##' (\cite{Basten et al., 2005}):
##' 
##' \emph{The linkage group is initiated with the pair of markers having the
##' smallest recombination fraction. The remaining markers are placed in a
##' \dQuote{pool} awaiting placement on the map. The linkage group is extended
##' by adding markers from the pool of unlinked markers. Each terminal marker
##' of the linkage group is a candidate for extension of the chain: The
##' unlinked marker that has the smallest recombination fraction with either is
##' added to the chain subject to the provision that the recombination fraction
##' is statistically significant at a prespecified level. This process is
##' repeated as long as markers can be added to the chain.}
##' 
##' After determining the order with \emph{RCD}, the final map is constructed
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
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{outcross} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{make.seq}}, \code{\link[onemap]{map}}
##' @references Basten, C. J., Weir, B. S. and Zeng, Z.-B. (2005) \emph{QTL
##' Cartographer Version 1.17: A Reference Manual and Tutorial for QTL
##' Mapping}.
##' 
##' Doerge, R. W. (1996) Constructing genetic maps by rapid chain delineation.
##' \emph{Journal of Quantitative Trait Loci} 2: 121-132.
##' 
##' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A. A. F.
##' (2009) Evaluation of algorithms used to order markers on genetics maps.
##' \emph{Heredity} 103: 494-502.
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
##'   LG1.rcd <- rcd(LG1)
##' 
##'   #F2 example
##'   data(fake.f2.onemap)
##'   twopt <- rf.2pts(fake.f2.onemap)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rcd <- rcd(LG1)
##'   LG1.rcd
##' }
##' 
rcd <-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5)
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
    diag(r)<-NA
    
    ## x defines the non-positioned markers
    x <- 1:n.mrk
    
    ## the group starts with the closest two markers
    i <- which(r==r[which.min(r)])
    i <- i[sample(length(i),1)]
    
    ## 'first' and 'last' are the markers on the edges of the ordered group
    first <- ceiling(i/n.mrk)
    last <- i-((first-1)*n.mrk)
    
    ## the two markers are set next to each other  
    order <- c(first,last)
    
    ## markers already ordered are marked as NaN
    x[first] <- NA
    x[last] <- NA
    
    ## extending the chain
    while (length(order) < n.mrk) {
        ## get the markers closest to the left end of the group
        j <- which(r[first,][x]==r[first,which.min(r[first,][x])])
        ## randomly choose one of them
        if (length(j) > 1) j <- j[sample(length(j),1)]
        
        ## get the markers closest to the right end of the group
        k <- which(r[last,][x]==r[last,which.min(r[last,][x])])
        ## randomly choose one of them
        if (length(k) > 1) k <- k[sample(length(k),1)]
        
        if (r[first,j] < r[last,k]) {  ## place a marker at the left side
            order <- c(j,order)
            x[j] <- NA
            first <- j
        }
    else if (r[first,j] > r[last,k]) {  ## or at the right side
        order <- c(order,k)
        x[k] <- NA
        last <- k } 
    else {  ## if the distance is the same, randomly choose one side
        rand <- sample(2,1) 
        if (rand==1) { order <- c(j,order); x[j] <- NA; first <- j }
      else { order <- c(order,k); x[k] <- NA; last <- k }
    }
    }
    ## end of chain
    cat("\norder obtained using RCD algorithm:\n\n", input.seq$seq.num[avoid.reverse(order)], "\n\ncalculating multipoint map using tol = ", tol, ".\n\n")
    map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(order)],twopt=input.seq$twopt), tol=tol)
}
## end of file





