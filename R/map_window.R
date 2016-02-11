#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: map_window.R                                                  ##
## Contains: map_window                                                ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido                        ##
## copyright (c) 2016, Gabriel R A Margarido                           ##
##                                                                     ##
## First version: 01/22/2016                                           ##
## Last update: 01/22/2016                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

## This function constructs the linkage map for a set of markers in a given
## order, using a sliding window to construct partial maps of adjacent
## subsets of markers


##' Construct the linkage map for a sequence of markers using sliding windows
##'
##' Estimates the multipoint log-likelihood, linkage phases and recombination
##' frequencies for a sequence of markers in a given order, using sliding
##' windows to construct partial maps of adjacent subsets of markers
##'
##' Markers are mapped in the order defined in the object \code{input.seq}.
##' A sliding window of size \code{window} runs along the sequence, in steps
##' of length given by argument \code{step} (which must be between \code{1}
##' and \code{window - 1}). For each window of markers, the best linkage phase
##' combination and recombination fractions are estimated and the multipoint
##' likelihood is calculated (for more details, check the documentation of
##' function \code{\link[onemap]{map}}).
##'
##' @param input.seq an object of class \code{sequence}.
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @param verbose If \code{TRUE}, print tracing information.
##' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions.} \item{seq.rf}{a \code{vector}
##' with the recombination frequencies between markers in the sequence.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{onemap} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.}
##' @author Gabriel R A Margarido \email{gramarga@@usp.br}
##' @seealso \code{\link[onemap]{make.seq}}. If you want to use a predefined
##' set of linkage phases, see \code{\link[onemap]{map}}
##' @keywords utilities
map_window <- function(input.seq, window = 40, step = 20, tol = 10E-5, verbose = FALSE)
{
    ## Check for correct object
    if(!is(input.seq, "sequence"))
        stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
    
    ## Check for appropriate arguments
    n.mar <- length(input.seq$seq.num)
    if(n.mar < 2) stop("The input sequence must have at least 2 markers.")
    if(length(window) != 1 || !is.numeric(window) || window < 2)
        stop("Argument \"window\" must be an integer greater than 1.")
    if(length(step) != 1 || !is.numeric(step) || step < 1)
        stop("Argument \"step\" must be an integer greater than 0.")

    ## Begin constructing map
    if (window >= n.mar) {
        ## Short sequence: no need to use sliding windows
        warning("Window size larger than sequence: estimating map for all markers.")
        map(input.seq = input.seq, tol = tol, verbose = verbose)
    }
    else {
        if (step >= window) {
            ## In this case, windows would not be adjacent or overlapping
            warning("Step size too large: using a step of ", window - 1, " marker(s).")
            step <- window - 1
        }
        
        starts <- seq(1, n.mar - window + 1, step)
        if ((tail(starts, 1) + window - 1) < n.mar) {
            ## We need to make sure the last window covers the last two markers
            ##starts <- c(starts, tail(starts, 1) + step)
            starts <- c(starts, n.mar - window + 1)
        }

        ## Allocate structures to store recombination fractions and linkage phases from partial maps
        counts <- rep(0, n.mar - 1)
        for (i in 1:length(starts)) {
            cur_end <- starts[i] + window - 1
            ## Number of partial maps in which a given marker interval appears
            counts[starts[i]:(cur_end - 1)] <- counts[starts[i]:(cur_end - 1)] + 1
        }
        recomb_fractions <- vector("list", n.mar - 1)
        link_phases <- vector("list", n.mar - 1)
        for (i in 1:(n.mar - 1)) {
            recomb_fractions[[i]] <- rep(NA, counts[i])
            link_phases[[i]] <- rep(NA, counts[i])
        }
        idx <- rep(1, n.mar - 1)

        ## Create partial maps and extract all necessary information
        for (i in 1:length(starts)) {
            ## Create auxiliary sequence with subset of markers
            cur_end <- starts[i] + window - 1
            temp_seq_num <- input.seq$seq.num[starts[i]:cur_end]
            if (length(input.seq$seq.phases) > 1 || input.seq$seq.phases != -1) {
                temp_seq_phase <- input.seq$seq.phases[starts[i]:(cur_end - 1)]
            }
            else {
                temp_seq_phase <- NULL
            }
            temp_seq <- make.seq(get(input.seq$data.name, pos = 1),
                                 arg = temp_seq_num,
                                 phase = temp_seq_phase,
                                 data.name = input.seq$data.name)
            
            #####################################################
            temp_seq$twopt <- "twopt" #### REMOVE THIS LINE
            #####################################################
            
            ## Construct partial linkage map
            temp_map <- map(temp_seq, tol = tol, verbose = verbose)
            print(temp_map)

            ## Extract results
            for (j in 1:length(temp_map$seq.rf)) {
                cur_interval <- starts[i] + j - 1
                recomb_fractions[[cur_interval]][idx[cur_interval]] <- temp_map$seq.rf[j]
                if (length(temp_map$seq.phases) > 1 || temp_map$seq.phases != -1) {
                    ## Linkage phases are only estimated for outcross datasets
                    link_phases[[cur_interval]][idx[cur_interval]] <- temp_map$seq.phases[j]
                }
                idx[cur_interval] <- idx[cur_interval] + 1
            }
        }

        ## Final recombination fractions: we use the mean of all partial maps for each interval
        seq.rf <- sapply(recomb_fractions, mean)
##################################################
        ## Should we check for non-informative markers? ##
##################################################

        print(link_phases)
        
        if (all(is.na(unlist(link_phases)))) {
            ## Linkage phases were not estimated
            seq.phases <- -1
        }
        else {
            if (any(is.na(unlist(link_phases)))) {
                ## If one linkage phase was estimated, there should be no missing data
                ## Should not get here
                stop("Problem estimating linkage phases.")
            }
            seq.phases <- rep(NA, n.mar - 1)
            for (i in 1:(n.mar - 1)) {
                temp_cur_phases <- unique(link_phases[[i]])
                if (length(temp_cur_phases) > 1) {
                    ## Check for redundant linkage phases
                    mrk1_type <- strsplit(get(input.seq$data.name, pos = 1)$segr.type[input.seq$seq.num[i]], "\\.")[[1]][1]
                    mrk2_type <- strsplit(get(input.seq$data.name, pos = 1)$segr.type[input.seq$seq.num[i+1]], "\\.")[[1]][1]
                    if (mrk1_type == "D1" || mrk2_type == "D1") {
                        link_phases[[i]][link_phases[[i]] == 2] <- 1
                        link_phases[[i]][link_phases[[i]] == 4] <- 3
                    }
                    if (mrk1_type == "D2" || mrk2_type == "D2") {
                        link_phases[[i]][link_phases[[i]] == 3] <- 1
                        link_phases[[i]][link_phases[[i]] == 4] <- 2
                    }
                    temp_cur_phases <- unique(link_phases[[i]])
                    if (length(temp_cur_phases) > 1) {
                        stop("Discordant linkage phases estimated. Try using a larger window size.")
                        ##print(i)
                        ##print(link_phases[[i]])
                        ##print(temp_cur_phases)
                        ##cat("\n")
                    }
                }
                seq.phases[i] <- temp_cur_phases
            }
        }
        
        structure(list(seq.num=input.seq$seq.num, seq.phases=seq.phases, seq.rf=seq.rf,
                       seq.like=NULL, data.name=input.seq$data.name, twopt=input.seq$twopt),
                  class = "sequence")
    }
}
## end of file
