#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: map.R                                                         ##
## Contains: map                                                       ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  ##
## with minor changes by Cristiane Taniguti
## copyright (c) 2009, Gabriel R A Margarido                           ##
##                                                                     ##
## First version: 02/27/2009                                           ##
## Last update: 01/14/2016                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

## This function constructs the linkage map for a set of markers in a given order


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
##' @param verbose If \code{TRUE}, print tracing information.
##' @param mds.seq When some pair of markers do not follow the linkage criteria, 
##' if \code{TRUE} one of the markers is removed and returns a vector with remaining 
##' marker numbers (useful for mds_onemap function).
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
##' @author Adapted from Karl Broman (package 'qtl') by Gabriel R A Margarido,
##' \email{gramarga@@usp.br} and Marcelo Mollinari, \email{mmollina@@gmail.com}
##' @seealso \code{\link[onemap]{make_seq}}
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
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'
##'   markers <- make_seq(twopt,c(30,12,3,14,2)) # correct phases
##'   map(markers)
##'
##'   markers <- make_seq(twopt,c(30,12,3,14,2),phase=c(4,1,4,3)) # incorrect phases
##'   map(markers)
##'@export
map <- function(input.seq,tol=10E-5, verbose=FALSE, mds.seq=FALSE)
{
  ## checking for correct object
  if(!("sequence" %in% class(input.seq)))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  ##Gathering sequence information
  seq.num<-input.seq$seq.num
  seq.phases<-input.seq$seq.phases
  seq.rf<-input.seq$seq.rf
  seq.like<-input.seq$seq.like
  ##Checking for appropriate number of markers
  if(length(seq.num) < 2) stop("The sequence must have at least 2 markers")
  ##For F2, BC and rils
  if(class(get(input.seq$data.name, pos=1))[2] == "f2")
  {
    final.map<-est_map_hmm_f2(geno=t(get(input.seq$data.name, pos=1)$geno[,seq.num]),
    error=t(get(input.seq$data.name, pos=1)$error[,seq.num]),
                              rf.vec=get_vec_rf_in(input.seq),
                              verbose=verbose,
                              tol=tol)
    return(structure(list(seq.num=seq.num,
                          seq.phases=seq.phases,
                          seq.rf=final.map$rf,
                          seq.like=final.map$loglike,
                          data.name=input.seq$data.name,
                          twopt=input.seq$twopt),
                     class = "sequence"))
  }
  else if(class(get(input.seq$data.name, pos=1))[2] == "backcross" ||
          class(get(input.seq$data.name, pos=1))[2] == "riself" ||
          class(get(input.seq$data.name, pos=1))[2] == "risib")
  {
    final.map<-est_map_hmm_bc(geno=t(get(input.seq$data.name, pos=1)$geno[,seq.num]),
    error=t(get(input.seq$data.name, pos=1)$error[,seq.num]),
                              rf.vec=get_vec_rf_in(input.seq),
                              verbose=verbose,
                              tol=tol)
    if(class(get(input.seq$data.name, pos=1))[2] == "riself" ||
       class(get(input.seq$data.name, pos=1))[2] == "risib")
      final.map$rf<-adjust_rf_ril(final.map$rf,
                                  type=class(get(input.seq$data.name, pos=1))[2],
                                  expand = FALSE)
    return(structure(list(seq.num=seq.num,
                          seq.phases=seq.phases,
                          seq.rf=final.map$rf,
                          seq.like=final.map$loglike,
                          data.name=input.seq$data.name,
                          twopt=input.seq$twopt), class = "sequence"))
  }
  
  if((seq.phases == -1) && (seq.rf == -1) && is.null(seq.like)) {
    ## if only the marker order is provided, without predefined linkage phases,
    ## a search for the best combination of phases is performed and recombination
    ## fractions are estimated
    seq.phase <- numeric(length(seq.num)-1)
    results <- list(rep(NA,4),rep(-Inf,4))
    
    ## linkage map is started with the first two markers in the sequence
    ## gather two-point information for this pair
    phase.init <- vector("list",1)
    list.init <- phases(make_seq(get(input.seq$twopt,pos = 1),seq.num[1:2],twopt=input.seq$twopt))
    phase.init[[1]] <- list.init$phase.init[[1]]
    Ph.Init <- comb_ger(phase.init)
    for(j in 1:nrow(Ph.Init)) {
      ## call to 'map' function with predefined linkage phase
      temp <- map(make_seq(get(input.seq$twopt),seq.num[1:2],phase=Ph.Init[j],twopt=input.seq$twopt))
      results[[1]][j] <- temp$seq.phases
      results[[2]][j] <- temp$seq.like
    }
    if (all(is.na(results[[2]]))) {
      if (mds.seq) {
        warning(cat("The linkage between markers", 
                    seq.num[mrk], "and", seq.num[mrk + 1], 
                    "did not reached the OneMap default criteria. They are probably segregating independently. Marker", 
                    seq.num[mrk + 1], "will be removed.\n"))
        return(seq.num[-(mrk + 1)])
        browser()
      }
      else {
        stop(paste("The linkage between markers", 
                   seq.num[mrk], "and", seq.num[mrk + 1], 
                   "did not reached the OneMap default criteria. They are probably segregating independently.\n"))
      }
    }
    seq.phase[1] <- results[[1]][which.max(results[[2]])] # best linkage phase is chosen
    
    if(length(seq.num) > 2) {
      ## for sequences with three or more markers, these are added sequentially
      for(mrk in 2:(length(seq.num)-1)) {
        results <- list(rep(NA,4),rep(-Inf,4))
        ## gather two-point information
        phase.init <- vector("list",mrk)
        list.init <- phases(make_seq(get(input.seq$twopt),c(seq.num[mrk],seq.num[mrk+1]),twopt=input.seq$twopt))
        phase.init[[mrk]] <- list.init$phase.init[[1]]
        for(j in 1:(mrk-1)) phase.init[[j]] <- seq.phase[j]
        Ph.Init <- comb_ger(phase.init)
        for(j in 1:nrow(Ph.Init)) {
          ## call to 'map' function with predefined linkage phases
          temp <- map(make_seq(get(input.seq$twopt),seq.num[1:(mrk+1)],phase=Ph.Init[j,],twopt=input.seq$twopt))
          results[[1]][j] <- temp$seq.phases[mrk]
          results[[2]][j] <- temp$seq.like
        }
        if(all(is.na(results[[2]]))) {
          if(mds.seq){
            warning(cat("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], "did not reached the OneMap default criteria. They are probably segregating independently. Marker", seq.num[mrk+1], "will be removed.\n"))
            return(seq.num[-(mrk+1)])
            browser()
          } else{
            stop(paste("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], "did not reached the OneMap default criteria. They are probably segregating independently.\n"))
          }
        }
        seq.phase[mrk] <- results[[1]][which.max(results[[2]])] # best combination of phases is chosen
      }
    }
    ## one last call to map function, with the final map
    map(make_seq(get(input.seq$twopt),seq.num,phase=seq.phase,twopt=input.seq$twopt))
  }
  else {
    ## if the linkage phases are provided but the recombination fractions have
    ## not yet been estimated or need to be reestimated, this is done here
    ## gather two-point information
    rf.init <- get_vec_rf_out(input.seq, acum=FALSE)
    ## estimate parameters
    final.map <- est_map_hmm_out(geno=t(get(input.seq$data.name, pos=1)$geno[,seq.num]),
    error = t(get(input.seq$data.name, pos=1)$error[,seq.num]),
                                 type=get(input.seq$data.name, pos=1)$segr.type.num[seq.num],
                                 phase=seq.phases,
                                 rf.vec=rf.init,
                                 verbose=FALSE,
                                 tol=tol)
    return(structure(list(seq.num=seq.num, seq.phases=seq.phases, seq.rf=final.map$rf,
                          seq.like=final.map$loglike, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "sequence"))
  }
}

## end of file
