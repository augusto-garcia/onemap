#######################################################################
##                                                                     ##
## Package: Onemap                                                     ##
##                                                                     ##
## File: seeded_map.R                                                  ##
## Contains: seeded_map                                                ##
##                                                                     ##
## Adaptaded by Cristiane Taniguti from BatchMap                       ##
## copyright (c) 2019 Cristiane Taniguti                               ##
##                                                                     ##
##                                                                     ##
## First version: 22/11/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' Construct the linkage map for a sequence of markers after seeding phases
##'
##' Estimates the multipoint log-likelihood, linkage phases and recombination
##' frequencies for a sequence of markers in a given order using seeded phases.
##'
##' Markers are mapped in the order defined in the object \code{input.seq}. The
##' best combination of linkage phases is also estimated starting from the first
##' position not in the given seeds.The multipoint likelihood is calculated
##' according to Wu et al. (2002b)(Eqs. 7a to 11), assuming that the
##' recombination fraction is the same in both parents. Hidden Markov chain
##' codes adapted from Broman et al. (2008) were used.
##'
##' @param input.seq an object of class \code{sequence}.
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @param verbose A logical, if TRUE it output progress status
##' information.
##' @param seeds A vector given the integer encoding of phases for the first
##' \emph{N} positions of the map
##' @param phase_cores The number of parallel processes to use when estimating
##' the phase of a marker. (Should be no more than 4)
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
##' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
##' if \code{TRUE} one of the markers is removed and map is performed again.
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
##' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
##' 2-point analyses.}
##' @author Adapted from Karl Broman (package 'qtl') by Gabriel R A Margarido,
##' \email{gramarga@@usp.br} and Marcelo Mollinari, \email{mmollina@@gmail.com}.
##' Modified to use seeded phases by Bastian Schiffthaler
##' \email{bastian.schiffthaler@umu.se}
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
##'   markers <- make_seq(twopt,c(30,12,3,14,2))
##'   seeded_map(markers, seeds = c(4,2))
##'
##' @export
seeded_map <- function(input.seq, tol=10E-5, phase_cores = 1,
                       seeds, verbose = F, rm_unlinked=F, parallelization.type = "PSOCK")
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
  
  seq.phase <- numeric(length(seq.num)-1)
  results <- list(rep(NA,4),rep(-Inf,4))
  #Add seeds as known phases
  seq.phase[1:length(seeds)] <- seeds
  
  #Skip all seeds i the phase estimation
  for(mrk in (length(seeds)+1):(length(seq.num)-1)) {
    results <- list(rep(NA,4),rep(-Inf,4))
    if(verbose)
    {
      message("Phasing marker ", input.seq$seq.num[mrk])
    }
    ## gather two-point information
    phase.init <- vector("list",mrk)
    list.init <- phases(make_seq(input.seq$twopt,
                                 c(seq.num[mrk],seq.num[mrk+1]),
                                 twopt=input.seq$twopt))
    phase.init[[mrk]] <- list.init$phase.init[[1]]
    for(j in 1:(mrk-1)) phase.init[[j]] <- seq.phase[j]
    Ph.Init <- comb_ger(phase.init)
    if(phase_cores == 1) {
      phases <-  lapply(1:nrow(Ph.Init), function(j) {
        map(make_seq(input.seq$twopt,
                     seq.num[1:(mrk+1)],
                     phase=Ph.Init[j,],
                     twopt=input.seq$twopt), 
            tol=tol, rm_unlinked = rm_unlinked)
      })
    } else {
      cl <- makeCluster(phase_cores, type = parallelization.type)
      clusterEvalQ(cl, c(library(onemap)))
      clusterExport(cl=cl, varlist=c('map'))
      phases <- parLapply(cl, 1:nrow(Ph.Init),
                          function(j) {
                            ## call to 'map' function with predefined linkage phases
                            map(make_seq(input.seq$twopt,
                                         seq.num[1:(mrk+1)],
                                         phase=Ph.Init[j,],
                                         twopt=input.seq$twopt), 
                                tol=tol, rm_unlinked = rm_unlinked, 
                                parallelization.type = parallelization.type)
                          })
      stopCluster(cl)
      gc(verbose = F)
    }
    if(!all(sapply(phases, function(x) is(x, "sequence")))) {
      if(rm_unlinked){
        warning(cat("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], 
                    "did not reached the OneMap default criteria. They are probably segregating independently. Marker", seq.num[mrk+1], "will be removed.
                    Use function map_avoid_unlinked to remove these markers automatically.\n"))
        return(seq.num[-(mrk+1)])
      } else{
        stop(paste("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], 
                   "did not reached the OneMap default criteria. They are probably segregating independently.
                   Use function map_avoid_unlinked to remove these markers automatically.\n"))
      }
    }
    
    for(j in 1:nrow(Ph.Init))
    {
      results[[1]][j] <- phases[[j]]$seq.phases[mrk]
      results[[2]][j] <- phases[[j]]$seq.like
    }
    if(all(is.na(results[[2]])))
    {
      if(rm_unlinked){
        warning(cat("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], 
                    "did not reached the OneMap default criteria. They are probably segregating independently. Marker", seq.num[mrk+1], "will be removed.
                    Use function map_avoid_unlinked to remove these markers automatically.\n"))
        return(seq.num[-(mrk+1)])
      } else{
        stop(paste("The linkage between markers", seq.num[mrk], "and", seq.num[mrk + 1], 
                   "did not reached the OneMap default criteria. They are probably segregating independently.
                   Use function map_avoid_unlinked to remove these markers automatically.\n"))
      }
    }
    # best combination of phases is chosen
    seq.phase[mrk] <- results[[1]][which.max(results[[2]])]
  }
  ## one last call to map function, with the final map
  map(make_seq(input.seq$twopt,seq.num,phase=seq.phase,
               twopt=input.seq$twopt), rm_unlinked = rm_unlinked, tol=tol, parallelization.type = parallelization.type)
}

## end of file
