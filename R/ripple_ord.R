#######################################################################
##                                                                     ##
## Package: BatchMap                                                     ##
##                                                                     ##
## File: ripple_ord                                                    ##
## Contains: generate_all, generate_rand, generate_one, ripple_window, ##
##           ripple_ord                                                ##
##                                                                     ##
## Written by Bastian Schiffthaler                                     ##
## copyright (c) 2017 Bastian Schiffthaler                             ##
##                                                                     ##
##                                                                     ##
## First version: 07/03/2017                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

generate_all <- function(input.seq, p, ws)
{
  all.ord <- t(apply(perm.tot(input.seq$seq.num[p:(p+ws-1)]),1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
  }))
  return(all.ord)
}

generate_rand <- function(input.seq, p, ws, n, pref)
{
  if(is.null(n))
  {
    n <- prod(ws:1)/2
  }
  
  ref <- input.seq$seq.num[p:(p+ws-1)]
  
  probs <- drop(cor(ref,t(perm.tot(input.seq$seq.num[p:(p+ws-1)])),
                    method = "spearman"))
  
  probs <- (probs - min(probs)) / (max(probs) - min(probs))
  probs[probs == 1] <- 0
  
  all.ord <- t(apply(perm.tot(input.seq$seq.num[p:(p+ws-1)]),1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
  }))
  
  if(pref == "similar")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE,probs),]
  }
  if(pref == "dissimilar")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE,-probs),]
  }
  if(pref == "neutral")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE),]
  }
  return(all.ord)
}

generate_one <- function(input.seq, p, ws, no_reverse)
{
  if(no_reverse)
  {
    all.ord <- matrix(NA,(sum((ws-1):1) + 1), ws)
  }
  else
  {
    all.ord <- matrix(NA,(sum((ws-1):1) + 1) * 2,ws)
  }
  all.ord[1,] <- input.seq$seq.num[p:(p+ws-1)]
  
  r <- 2
  for(i in 1:(ws - 1))
  {
    for(j in ws:(i+1))
    {
      all.ord[r,] <- all.ord[1,]
      tmp <- all.ord[r,i]
      all.ord[r,i] <- all.ord[r,j]
      all.ord[r,j] <- tmp
      r <- r+1
    }
  }
  if(! no_reverse)
  {
    for(r in (sum((ws-1):1) + 2):((sum((ws-1):1) + 1) * 2))
    {
      all.ord[r,] <- rev(all.ord[r - sum((ws-1):1) + 1,])
    }
  }
  all.ord <- all.ord[!duplicated(apply(all.ord,1,paste,collapse="-")),]
  all.ord <- all.ord[-1,]
  all.ord <- t(apply(all.ord,1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
  }))
  return(all.ord)
}


#' Update linkage map with alternative orders at a given position
#'
#' This function carries out re-ordering of one single window according to user
#' defined criteria. The best order is chosen based on the difference in log
#' likelihood. Different heuristics are avaible to select which orders to test.
#' Note that testing all orders has factorial complexity (N!/2) meaning that
#' it's not feasible for window sizes larger than 6.
#'
#' Methods:
#' \emph{all:}{Checks for all possible permutations in the window. Will be very
#' very slow for large window size.}
#' \emph{one:}{Checks for all possible pairwise switches in the window. The
#' complexity scales as \code{sum(ws:1)}, \code{ws} being the window size.}
#' \emph{rand:}{First, generates all possible permutations. Then samples
#' \code{n} sequences from those and tests those. Time complexity is N.}
#'
#' The "rand" method can be further tuned to preferentially select similar,
#' dissimilar sequences or to perform unbiased sampling. In the first two cases
#' sampling probability is adjusted via a spearman correlation of the sequences
#' to all possible sequences.
#'
#' @param input.seq An object of class \code{sequence}.
#' @param ws The window size in which ti consider re-ordering
#' @param tol The tolerance for checking convergence of the EM model
#' @param phase_cores The number of parallel processes to use to estimate phases.
#' Should be no higher than 4.
#' @param ripple_cores The number of parallel processes that should be used when
#' testing different marker orders.
#' @param start The position of the first marker of the window within the
#' input sequence
#' @param verbose A character vector that includes any or all of "batch",
#' "order", "position", "time" and "phase" to output progress status
#' information.
#' @param type One of "one", "all" or "rand".
#' @param n For method "rand": The number of random samples to be tested.
#' @param pref For method "rand": One of "similar", "dissimilar" or "neutral".
#' Controls if sampling probability should be adjusted based on similarity to
#' the input sequence. See description.
#' @param no_reverse For method "one": If \code{FALSE}, the method will also
#' create all possible reverse sequences if the marker swaps.
#' @return An object of class \code{sequence} that is the best order for the
#' re-ordered window within the input sequence.
#' @import parallel
#' @export
ripple_window <- function(input.seq, ws=4, tol=10E-4, phase_cores = 1,
                          ripple_cores = 1, start = 1, verbose = NULL,
                          type = "one", n = NULL, pref = NULL,
                          no_reverse = TRUE) {
  
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) {
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'sequence'")
  }
  if(ws < 2) stop("ws must be greater than or equal to 2")
  len <- length(input.seq$seq.num)
  ## computations unnecessary in this case
  if (len <= ws) stop("Length of sequence ",
                      deparse(substitute(input.seq)),
                      " is smaller than ws. You can use the ",
                      "compare function instead")
  
  ## allocate variables
  rf.init <- rep(NA,len-1)
  phase <- rep(NA,len-1)
  # tot <- prod(1:ws)
  # best.ord.phase <- matrix(NA,tot,len-1)
  # best.ord.like <- best.ord.LOD <- rep(-Inf,tot)
  # all.data <- list()
  
  ## gather two-point information
  list.init <- phases(input.seq)
  
  p <- start
  if("order" %in% verbose)
  {
    message(p-1, "...", input.seq$seq.num[p-1], "|",
            paste(input.seq$seq.num[p:(p+ws-1)], collapse = "-"),"|",
            input.seq$seq.num[p+ws],"...", p+ws+1)
  }
  
  if(type == "all") all.ord <- generate_all(input.seq, p, ws)
  if(type == "rand") all.ord <- generate_rand(input.seq, p, ws, n, pref)
  if(type == "one") all.ord <- generate_one(input.seq, p, ws, no_reverse)
  
  poss <- mclapply(1:nrow(all.ord), mc.allow.recursive = TRUE,
                   mc.cores = ripple_cores, function(i){
                     if("position" %in% verbose)
                     {
                       message("Trying order ",i," of ",nrow(all.ord),
                               " for start position ",p)
                     }
                     mp <- list(seq.like = -Inf, seq.rf = -1)
                     tryCatch({
                       if(start > 1)
                       {
                         seeds <- input.seq$seq.phases[1:(start - 1)]
                         mp <- seeded_map(make_seq(input.seq$twopt,
                                                   all.ord[i,]),
                                          phase_cores = phase_cores,
                                          verbose = F,
                                          seeds = seeds)
                       }
                       else
                       {
                         mp <- map(make_seq(input.seq$twopt, all.ord[i,]),
                                   phase_cores = phase_cores,
                                   verbose = F)
                       }
                     }, error = function(e){},
                     finally = {
                       return(mp)
                     })
                   })
  
  
  
  best <- which.max(sapply(poss,"[[","seq.like"))
  if(poss[[best]]$seq.like > input.seq$seq.like)
  {
    return(poss[[best]])
  } else {
    return(input.seq)
  }
  
  
}


#' Update linkage map with alternative orders at all positions
#'
#' This function carries out re-ordering of all markers using a sliding window
#' according to user defined criteria. The best order is chosen based on the
#' difference in log likelihood. Different heuristics are avaible to select
#' which orders to test. Note that testing all orders has factorial complexity
#' (N!/2) meaning that it's not feasible for window sizes larger than 6.
#'
#' Methods:
#' \emph{all:}{Checks for all possible permutations in the window. Will be very
#' very slow for large window size.}
#' \emph{one:}{Checks for all possible pairwise switches in the window. The
#' complexity scales as \code{sum(ws:1)}, \code{ws} being the window size.}
#' \emph{rand:}{First, generates all possible permutations. Then samples
#' \code{n} sequences from those and tests those. Time complexity is N.}
#' The "rand" method can be further tuned to preferentially select similar,
#' dissimilar sequences or to perform unbiased sampling. In the first two cases
#' sampling probability is adjusted via a spearman correlation of the sequences
#' to all possible sequences.
#'
#' @param input.seq An object of class \code{sequence}.
#' @param ws The window size in which ti consider re-ordering
#' @param tol The tolerance for checking convergence of the EM model
#' @param phase_cores The number of parallel processes to use to estimate phases.
#' Should be no higher than 4.
#' @param ripple_cores The number of parallel processes that should be used when
#' testing different marker orders.
#' @param start The position of the first marker of the window within the
#' input sequence
#' @param verbose A character vector that includes any or all of "batch",
#' "order", "position", "time" and "phase" to output progress status
#' information.
#' @param method One of "one", "all" or "rand". Which algorithm to use to
#' generate alternative orders to test (see Details)
#' @param n For method "rand": The number of random samples to be tested.
#' @param pref For method "rand": One of "similar", "dissimilar" or "neutral".
#' Controls if sampling probability should be adjusted based on similarity to
#' the input sequence. See description.
#' @param batches List with batches that are being processed. Only used when
#' estimating finish time.
#' @param no_reverse For method "one": If \code{FALSE}, the method will also
#' create all possible reverse sequences if the marker swaps.
#' @return An object of class \code{sequence} that is the best order for the
#' re-ordered input sequence.
#' @import parallel
#' @export
ripple_ord <- function(input.seq, ws = 4, tol = 10E-4, phase_cores = 1,
                       ripple_cores = 1, method = "one", n = NULL,
                       pref = "neutral", start = 1, verbose = NULL,
                       batches = NULL, no_reverse = TRUE)
{
  LG <- input.seq
  if(start + ws > length(input.seq$seq.num)) return(LG)
  
  if("time" %in% verbose){
    message("Progress:\n")
    pb <- txtProgressBar(min = start, max = (length(input.seq$seq.num) - ws), style = 3)
  }
  for(i in start:(length(input.seq$seq.num) - ws))
  {
    LG <- ripple_window(input.seq = LG, ws = ws, tol = tol,
                        phase_cores = phase_cores, ripple_cores = ripple_cores,
                        start = i, verbose = verbose,
                        no_reverse = no_reverse, type = method, n = n,
                        pref = pref)
    if("time" %in% verbose)
      setTxtProgressBar(pb, i)
  }
  close(pb)
  return(LG)
}

predict_time <- function(batches, ws, timings)
{
  if(is.null(batches)) return(NULL)
  m <- median(timings)
  s <- sd(timings)
  overlap <- length(intersect(batches[[1]], batches[[2]]))
  tot <- length(batches[[1]]) - ws
  for(f in 2:length(batches))
  {
    tot <- tot + length(batches[[f]] - overlap - ws)
  }
  m <- round(m*tot)
  s <- round(s*tot)
  message("The current best estimate for the total time is ", m, "s +/- ",
          s, "s.")
  return(c(m,s))
}