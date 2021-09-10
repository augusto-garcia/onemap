#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: map.R                                                         ##
## Contains: map, map_save_ram, map_avoid_unlinked                     ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido, Marcelo Mollinari and ##
## Cristiane Taniguti                                                  ##
## copyright (c) 2009, Gabriel R A Margarido                           ##
##                                                                     ##
## First version: 02/27/2009                                           ##
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
##' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
##' if \code{TRUE} one of the markers is removed and returns a vector with remaining 
##' marker numbers (useful for mds_onemap and map_avoid_unlinked functions).
##' @param phase_cores number of computer cores to be used in analysis
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
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
##' \email{gramarga@@usp.br} and Marcelo Mollinari, \email{mmollina@@gmail.com},
##' with minor changes by Cristiane Taniguti and Bastian Schiffthaler 
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
##' \dontrun{
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'
##'   markers <- make_seq(twopt,c(30,12,3,14,2)) # correct phases
##'   map(markers)
##'
##'   markers <- make_seq(twopt,c(30,12,3,14,2),phase=c(4,1,4,3)) # incorrect phases
##'   map(markers)
##'   }
##'@import parallel
##'
##'@export
map <- function(input.seq,tol=10E-5, verbose=FALSE, 
                rm_unlinked=FALSE, phase_cores = 1, 
                parallelization.type = "PSOCK")
{
  ## checking for correct object
  if(!(is(input.seq, "sequence")))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  ## Checking phase_cores
  if(is(input.seq$data.name, c("riself", "risib", "backcross")) & phase_cores != 1){
    warning("For RILs and backcross populations, we do not need to estimate phase. Therefore, the parallelization is not possible with our approach.")
    phase_cores <- 1
  }
  ## Checking version
  if(!any(names(input.seq$data.name) %in% "error")){
    input.seq$data.name <- create_probs(input.seq$data.name, global_error = 10^(-5))
  }
  
  ##Gathering sequence information
  seq.num<-input.seq$seq.num
  seq.phases<-input.seq$seq.phases
  seq.rf<-input.seq$seq.rf
  seq.like<-input.seq$seq.like
  ##Checking for appropriate number of markers
  if(length(seq.num) < 2) stop("The sequence must have at least 2 markers")
  ##For F2, BC and rils
  if(is(input.seq$data.name, c("backcross", "riself", "risib")))
  {
    final.map<-est_map_hmm_bc(geno=t(input.seq$data.name$geno[,seq.num]),
                              error=input.seq$data.name$error[seq.num + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(seq.num)),],
                              rf.vec=get_vec_rf_in(input.seq),
                              verbose=verbose,
                              tol=tol)
    
    if(is(input.seq$data.name, c("riself", "risib")))
      final.map$rf<-adjust_rf_ril(final.map$rf,
                                  type=class(input.seq$data.name)[2],
                                  expand = FALSE)
    return(structure(list(seq.num=seq.num,
                          seq.phases=seq.phases,
                          seq.rf=final.map$rf,
                          seq.like=final.map$loglike,
                          data.name=input.seq$data.name,
                          probs = final.map$probs,
                          twopt=input.seq$twopt), class = "sequence"))
  }
  
  if(all(seq.phases == -1) && all(seq.rf == -1) && all(is.null(seq.like))) {
    ## if only the marker order is provided, without predefined linkage phases,
    ## a search for the best combination of phases is performed and recombination
    ## fractions are estimated
    seq.phase <- numeric(length(seq.num)-1)
    results <- list(rep(NA,4),rep(-Inf,4))
    
    ## linkage map is started with the first two markers in the sequence
    ## gather two-point information for this pair
    phase.init <- vector("list",1)
    list.init <- phases(make_seq(input.seq$twopt,seq.num[1:2],twopt=input.seq$twopt))
    phase.init[[1]] <- list.init$phase.init[[1]]
    Ph.Init <- comb_ger(phase.init)
    if(phase_cores == 1) {
      phases <-  lapply(1:nrow(Ph.Init), function(j) {
        map(make_seq(input.seq$twopt,
                     seq.num[1:2],
                     phase=Ph.Init[j]), 
            tol=tol, 
            rm_unlinked = rm_unlinked)
      })
    } else {
      cl <- makeCluster(phase_cores, type = parallelization.type)
      clusterEvalQ(cl, c(library(onemap)))
      clusterExport(cl=cl, varlist=c('map'))
      phases <- parLapply(cl, 1:nrow(Ph.Init), 
                          function(j) {
                            ## call to 'map' function with predefined linkage phase
                            map(make_seq(input.seq$twopt,
                                         seq.num[1:2],
                                         phase=Ph.Init[j]), 
                                tol=tol, 
                                rm_unlinked = rm_unlinked,
                                parallelization.type = parallelization.type)
                          })
      stopCluster(cl)
      gc(verbose = F)
    }
    if(!all(sapply(phases, function(x) is(x, "sequence")))){
      if (rm_unlinked) {
        warning(cat("The linkage between markers", 
                    seq.num[1], "and", seq.num[2], 
                    "did not reached the OneMap default criteria. They are probably segregating independently. Marker", 
                    seq.num[2], "will be removed. Use argument rm_unlinked = TRUE if you are ordering markers or function map_avoid_unlinked to remove these markers automatically.\n"))
        return(seq.num[-2])
      }
      else {
        stop(paste("The linkage between markers", 
                   seq.num[1], "and", seq.num[2], 
                   "did not reached the OneMap default criteria. They are probably segregating independently.
                   Use argument rm_unlinked = TRUE if you are ordering markers or function map_avoid_unlinked to remove these markers automatically.\n"))
      }
    }
    for(j in 1:nrow(Ph.Init)) {
      ## call to 'map' function with predefined linkage phase
      #temp <- map(make_seq(input.seq$twopt,seq.num[1:2],phase=Ph.Init[j],twopt=input.seq$twopt))
      results[[1]][j] <- phases[[j]]$seq.phases
      results[[2]][j] <- phases[[j]]$seq.like
    }
    if (all(is.na(results[[2]]))) {
      if (rm_unlinked) {
        warning(cat("The linkage between markers", 
                    seq.num[1], "and", seq.num[2], 
                    "did not reached the OneMap default criteria. They are probably segregating independently. Marker", 
                    seq.num[2], "will be removed. Use function map_avoid_unlinked to remove these markers automatically.\n"))
        return(seq.num[-(2)])
      }
      else {
        stop(paste("The linkage between markers", 
                   seq.num[1], "and", seq.num[2], 
                   "did not reached the OneMap default criteria. They are probably segregating independently. 
                   Use function map_avoid_unlinked to remove these markers automatically.\n"))
      }
    }
    seq.phase[1] <- results[[1]][which.max(results[[2]])] # best linkage phase is chosen
    
    if(length(seq.num) > 2) {
      ## for sequences with three or more markers, these are added sequentially
      for(mrk in 2:(length(seq.num)-1)) {
        results <- list(rep(NA,4),rep(-Inf,4))
        ## gather two-point information
        phase.init <- vector("list",mrk)
        list.init <- phases(make_seq(input.seq$twopt,c(seq.num[mrk],seq.num[mrk+1]),twopt=input.seq$twopt))
        phase.init[[mrk]] <- list.init$phase.init[[1]]
        for(j in 1:(mrk-1)) phase.init[[j]] <- seq.phase[j]
        Ph.Init <- comb_ger(phase.init)
        if(phase_cores == 1) {
          phases <-  lapply(1:nrow(Ph.Init), function(j) {
            map(make_seq(input.seq$twopt,
                         seq.num[1:(mrk+1)],
                         phase=Ph.Init[j,]), 
                tol=tol, 
                rm_unlinked = rm_unlinked)
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
                                             phase=Ph.Init[j,]), 
                                    tol=tol, 
                                    rm_unlinked = rm_unlinked,
                                    parallelization.type = parallelization.type)
                              })
          stopCluster(cl)
          gc(verbose = F)
        }
        if(!all(sapply(phases, function(x) is(x, "sequence")))){
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
        for(j in 1:nrow(Ph.Init)) {
          ## call to 'map' function with predefined linkage phases
          #temp <- map(make_seq(input.seq$twopt,seq.num[1:(mrk+1)],phase=Ph.Init[j,],twopt=input.seq$twopt))
          results[[1]][j] <- phases[[j]]$seq.phases[mrk]
          results[[2]][j] <- phases[[j]]$seq.like
        }
        if(all(is.na(results[[2]]))) {
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
        seq.phase[mrk] <- results[[1]][which.max(results[[2]])] # best combination of phases is chosen
      }
    }
    ## one last call to map function, with the final map
    map(make_seq(input.seq$twopt,
                 seq.num,
                 phase=seq.phase), 
        tol=tol, 
        rm_unlinked = rm_unlinked,  parallelization.type = parallelization.type)
  }
  else {
    ## if the linkage phases are provided but the recombination fractions have
    ## not yet been estimated or need to be reestimated, this is done here
    ## gather two-point information
    rf.init <- get_vec_rf_out(input.seq, acum=FALSE)
    if(any(is.na(rf.init))) {
      stop("Linkage criterias could not be reached")
    }
    ## estimate parameters
    final.map <- est_map_hmm_out(geno=t(input.seq$data.name$geno[,seq.num]),
                                 error = input.seq$data.name$error[seq.num + 
                                                                     rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, 
                                                                         each=length(seq.num)),],
                                 type=input.seq$data.name$segr.type.num[seq.num],
                                 phase=seq.phases,
                                 rf.vec=rf.init,
                                 verbose=FALSE,
                                 tol=tol)
    
    return(structure(list(seq.num=seq.num, seq.phases=seq.phases, seq.rf=final.map$rf,
                          seq.like=final.map$loglike, data.name=input.seq$data.name, 
                          probs = final.map$probs, twopt=input.seq$twopt), class = "sequence"))
  }
}

#' Perform map using background objects with only selected markers. It saves ram memory during the procedure.
#' It is useful if dealing with many markers in total data set.
#' 
#' @param input.seq object of class sequence
#' @param size The center size around which an optimum is to be searched
#' @param overlap The desired overlap between batches
#' @param phase_cores The number of parallel processes to use when estimating
#' the phase of a marker. (Should be no more than 4)
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
#' @param tol tolerance for the C routine, i.e., the value used to evaluate
#' convergence.
##' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
##' if \code{TRUE} one of the markers is removed and returns a vector with remaining 
##' marker numbers (useful for mds_onemap and map_avoid_unlinked functions).
##' @param verbose If \code{TRUE}, print tracing information.
##' 
map_save_ram <- function(input.seq,
                         tol=10E-5, 
                         verbose=FALSE, 
                         rm_unlinked=FALSE, 
                         phase_cores = 1, 
                         size = NULL, 
                         overlap = NULL,
                         parallelization.type = "PSOCK"){
  
  input.seq.tot <- input.seq
  input.seq_ram <- input.seq
  if(length(input.seq$seq.num) < input.seq.tot$data.name$n.mar){
    split.twopts <- split_2pts(twopts.obj = input.seq$twopt, mks = input.seq$seq.num) 
    input.seq_ram <- make_seq(split.twopts, "all")
  }
  if(phase_cores == 1){
    return.map <- map(input.seq = input.seq_ram, tol = tol, 
                      verbose = verbose, 
                      rm_unlinked = rm_unlinked, 
                      phase_cores = phase_cores,
                      parallelization.type = parallelization.type)
  } else {
    if(is.null(size) | is.null(overlap)){
      stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you should also define `size` and `overlap` arguments. See ?map_avoid_unlinked and ?pick_batch_sizes")
    } else {
      return.map <- map_overlapping_batches(input.seq = input.seq_ram,
                                            size = size, overlap = overlap, 
                                            phase_cores = phase_cores, 
                                            tol=tol, rm_unlinked = rm_unlinked,
                                            parallelization.type = parallelization.type)
    }
  }
  if(length(input.seq_ram$seq.num) < input.seq.tot$data.name$n.mar){
    if(!is(return.map, "integer")){ # When rm_unlinked == F
      return.map$seq.num <- input.seq.tot$seq.num
      return.map$data.name <- input.seq.tot$data.name
      return.map$twopt <- input.seq.tot$twopt
    } else {
      remain <- colnames(input.seq_ram$data.name$geno)[return.map]
      old <- colnames(input.seq.tot$data.name$geno)[input.seq.tot$seq.num]
      return.map <- input.seq.tot$seq.num[old %in% remain] 
    }
  }
  return(return.map)
}

#' Repeat HMM if map find unlinked maker
#'
#' @param input.seq object of class sequence
#' @param size The center size around which an optimum is to be searched
#' @param overlap The desired overlap between batches
#' @param phase_cores The number of parallel processes to use when estimating
#' the phase of a marker. (Should be no more than 4)
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
#' @param tol tolerance for the C routine, i.e., the value used to evaluate
#' convergence.
#' 
#' @export
map_avoid_unlinked <- function(input.seq, 
                               size = NULL, 
                               overlap = NULL,
                               phase_cores = 1, 
                               tol = 10E-5,
                               parallelization.type = "PSOCK"){
  #TODO: error checks...
  map_df <- map_save_ram(input.seq, rm_unlinked = T, 
                         size = size, 
                         overlap = overlap, 
                         tol=tol, 
                         phase_cores = phase_cores,
                         parallelization.type = parallelization.type)
  
  while(is(map_df, "integer")){
    seq_true <- make_seq(input.seq$twopt, map_df)
    map_df <- map_save_ram(input.seq = seq_true, 
                           rm_unlinked = T, 
                           tol=tol, 
                           size = size, 
                           overlap = overlap, 
                           phase_cores = phase_cores,
                           parallelization.type = parallelization.type)
  }
  return(map_df)
}
