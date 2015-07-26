#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                      #
# File: map.R                                                         #
# Contains: map                                                       #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function constructs the linkage map for a set of markers in a given order


##' Construct the linkage map for a sequence of markers
##' 
##' Estimates the multipoint log-likelihood, linkage phases and recombination
##' frequencies for a sequence of markers in a given order.
##' 
##' Markers are mapped in the order defined in the object \code{input.seq}. If
##' this object also contains a user-defined combination of linkage phases,
##' recombination frequencies and log-likelihood are estimated for that
##' particular case. Otherwise, the best linkage phase combination is also
##' estimated. The multipoint likelihood is calculated according to Wu et al.
##' (2002b)(Eqs. 7a to 11), assuming that the recombination fraction is the
##' same in both parents. Hidden Markov chain codes adapted from Broman et al.
##' (2008) were used.
##' 
##' @param input.seq an object of class \code{sequence}.
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
##' @author Adapted from Karl Broman (package 'qtl') by Gabriel R A Margarido,
##' \email{gramarga@@usp.br} and Marcelo Mollinari, \email{mmollina@@gmail.com}
##' @seealso \code{\link[onemap]{make.seq}}
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##' 
##' Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci with
##' dominant and missing markers in various crosses from two inbred lines.
##' \emph{Genetica} 101: 47-58.
##' 
##' Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J., Lincoln,
##' S. E. and Newburg, L. (1987) MAPMAKER: An interactive computer package for
##' constructing primary genetic linkage maps of experimental and natural
##' populations. \emph{Genomics} 1: 174-181.
##' 
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##' 
##' Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage mapping of
##' sex-specific differences. \emph{Genetical Research} 79: 85-96
##' @keywords utilities
##' @examples
##' 
##'   data(example.out)
##'   twopt <- rf.2pts(example.out)
##' 
##'   markers <- make.seq(twopt,c(30,12,3,14,2)) # correct phases
##'   map(markers)
##' 
##'   markers <- make.seq(twopt,c(30,12,3,14,2),phase=c(4,1,4,3)) # incorrect phases
##'   map(markers)
##' 
map <-
function(input.seq,tol=10E-5) {
  # checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequnece'")
  if(length(input.seq$seq.num) < 2) stop("The sequence must have at least 2 markers")

  ##For F2, BC and rils
  if(class(get(input.seq$data.name, pos=1))=="f2.onemap"){
    ph <- numeric(length(input.seq$seq.num) - 1)
    for(i in 1:(length(input.seq$seq.num)-1)) {
      if(input.seq$seq.num[i] > input.seq$seq.num[i+1])
        ph[i] <- which.max(get(input.seq$twopt)$analysis[acum(input.seq$seq.num[i]-2)+input.seq$seq.num[i+1],,2])
      else
        ph[i] <- which.max(get(input.seq$twopt)$analysis[acum(input.seq$seq.num[i+1]-2)+input.seq$seq.num[i],,2])
    }
    input.seq$seq.phases<-ph
  }
  else if(class(get(input.seq$data.name, pos=1))=="bc.onemap" || class(get(input.seq$data.name, pos=1))=="riself.onemap" || class(get(input.seq$data.name, pos=1))=="risib.onemap"){
    ph <- rep(1,(length(input.seq$seq.num) - 1))
    input.seq$seq.phases<-ph
  }
  
  if((input.seq$seq.phases == -1) && (input.seq$seq.rf == -1) && is.null(input.seq$seq.like)) {
    ## if only the marker order is provided, without predefined linkage phases,
    ## a search for the best combination of phases is performed and recombination
    ## fractions are estimated
    seq.phase <- numeric(length(input.seq$seq.num)-1)
    results <- list(rep(NA,4),rep(-Inf,4))
    
    ## linkage map is started with the first two markers in the sequence
    ## gather two-point information for this pair
    phase.init <- vector("list",1)
	list.init <- phases(make.seq(get(input.seq$twopt),c(input.seq$seq.num[1],input.seq$seq.num[2]),twopt=input.seq$twopt))
    phase.init[[1]] <- list.init$phase.init[[1]]
    Ph.Init <- comb.ger(phase.init)  
    for(j in 1:nrow(Ph.Init)) {
      ## call to 'map' function with predefined linkage phase
      temp <- map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],phase=Ph.Init[j],twopt=input.seq$twopt))
      results[[1]][j] <- temp$seq.phases
      results[[2]][j] <- temp$seq.like
    }
    seq.phase[1] <- results[[1]][which.max(results[[2]])] # best linkage phase is chosen
    
    if(length(input.seq$seq.num) > 2) {
      ## for sequences with three or more markers, these are added sequentially
      for(mrk in 2:(length(input.seq$seq.num)-1)) {
        results <- list(rep(NA,4),rep(-Inf,4))
        ## gather two-point information
        phase.init <- vector("list",mrk)
        list.init <- phases(make.seq(get(input.seq$twopt),c(input.seq$seq.num[mrk],input.seq$seq.num[mrk+1]),twopt=input.seq$twopt))
        phase.init[[mrk]] <- list.init$phase.init[[1]]
        for(j in 1:(mrk-1)) phase.init[[j]] <- seq.phase[j]
        Ph.Init <- comb.ger(phase.init)      
        for(j in 1:nrow(Ph.Init)) {
          ## call to 'map' function with predefined linkage phases
          temp <- map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:(mrk+1)],phase=Ph.Init[j,],twopt=input.seq$twopt))
          results[[1]][j] <- temp$seq.phases[mrk]
          results[[2]][j] <- temp$seq.like
        }
        seq.phase[mrk] <- results[[1]][which.max(results[[2]])] # best combination of phases is chosen
      }
    }
    ## one last call to map function, with the final map
    map(make.seq(get(input.seq$twopt),input.seq$seq.num,phase=seq.phase,twopt=input.seq$twopt))
  }
  
  else {
    ## if the linkage phases are provided but the recombination fractions have
    ## not yet been estimated or need to be reestimated, this is done here
    ## gather two-point information
    rf.init <- numeric(length(input.seq$seq.num)-1)
    for(i in 1:(length(input.seq$seq.num)-1)) {
      if(input.seq$seq.num[i] > input.seq$seq.num[i+1])
        rf.init[i] <- get(input.seq$twopt)$analysis[acum(input.seq$seq.num[i]-2)+input.seq$seq.num[i+1],input.seq$seq.phases[i],1]
      else
        rf.init[i] <- get(input.seq$twopt)$analysis[acum(input.seq$seq.num[i+1]-2)+input.seq$seq.num[i],input.seq$seq.phases[i],1]
    }
    ## estimate parameters
    final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,input.seq$seq.num],
                           type=get(input.seq$data.name, pos=1)$segr.type.num[input.seq$seq.num],
                           phase=input.seq$seq.phases,
                           rec=rf.init,
                           verbose=FALSE,
                           tol=tol)
    
    if(class(get(input.seq$data.name, pos=1))=="riself.onemap" || class(get(input.seq$data.name, pos=1))=="risib.onemap")
      final.map$rf<-adjust.rf.ril(final.map$rf, type=class(get(input.seq$data.name, pos=1)), expand = FALSE)
    
    structure(list(seq.num=input.seq$seq.num, seq.phases=input.seq$seq.phases, seq.rf=final.map$rf,
                   seq.like=final.map$loglike, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "sequence")
  }
}

## end of file
