#######################################################################
##                                                                     ##
## Package: OneMap                                                     ##
##                                                                     ##
## File: overlapping_batches.R                                         ##
## Contains: map_overlapping_batches, generate_overlapping_batches,    ##
##          pick_batch_sizes                                           ##
##                                                                     ##
## Adaptated by Cristiane Taniguti from BatchMap package               ##
## copyright (c) 2019 Cristiane Taniguti                               ##
##                                                                     ##
##                                                                     ##
## First version: 13/11/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' @export
generate_overlapping_batches <- function(input.seq, size = 50, overlap = 15,
                                         silent = FALSE)
{
  start <- 1
  end <- size
  current <- 1
  res <- list()
  while(end <= length(input.seq$seq.num))
  {
    res[[current]] <- input.seq$seq.num[start:end]
    current <- current + 1
    start <- end - overlap
    if(end == length(input.seq$seq.num)) break
    end <- end + size - overlap - 1
    if(end > length(input.seq$seq.num)) end <- length(input.seq$seq.num)
  }
  sizes <- unlist(lapply(res, length))
  if(length(sizes) < 2 & ! silent)
  {
    warning("You should at least have two overlapping batches.",
            " Reconsider the size parameter.")
  }
  if(any(sizes/size > 1.25) & ! silent)
  {
    warning("One group is 25% bigger than the group size. ",
            "Consider adjusting parameters.")
  }
  return(res)
}


##' Picking optimal batch size values
##'
##' Suggest an optimal batch size value for use in
##' \code{\link[BatchMap]{map_overlapping_batches}}
##'
##' @param input.seq an object of class \code{sequence}.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param around The range around the center which is maximally allowed
##' to be searched.
##' @return An integer value for the size which most evenly divides batches. In
##' case of ties, bigger batch sizes are preferred.
##' @author Bastian Schiffthaler, \email{bastian.schiffthaler@umu.se}
##' @seealso \code{\link[BatchMap]{map_overlapping_batches}}
##'
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'   LG <- structure(list(seq.num = seq(1,800)), class = "sequence")
##'   batchsize <- pick_batch_sizes(LG, 50, 19)
##' }
##' @export
pick_batch_sizes <- function(input.seq, size = 50, overlap = 15, around = 5)
{
  test.sizes <- c(size, (size - 1):(size - around), (size + 1):(size + around))
  all.batches <- lapply(test.sizes, function(s){
    generate_overlapping_batches(input.seq, s, overlap, silent = TRUE)
  })
  x <- unlist(lapply(all.batches, function(f){
    ran <- range(unlist(lapply(f,length)))
    ran[2] - ran[1]
  }))
  x <- which(x == min(x))
  test.sizes[x[length(x)]] #prefer larger maps
}

##' Mapping overlapping batches
##'
##' Apply the batch mapping algorithm using overlapping windows.
##'
##' This algorithm implements the overlapping batch maps for high density
##' marker sets. The mapping problem is reduced to a number of subsets (batches)
##' which carry information forward in order to more accurately estimate
##' recombination fractions and phasing. Further the user has the option of
##' setting \code{fun.order} to a function that tries different orders and
##' iteratively reorders markers to improve the map. See
##' \code{\link[BatchMap]{ripple.ord}} for such an implementation. The ordering
##' function is triggered at least \code{min.tries} times per batch, or as
##' long as a batch has two markers with a distance greater than
##' \code{max.dist}.
##'
##' @param input.seq an object of class \code{sequence}.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param fun.order A function that is applied to each batch to improve
##' marker order. See \code{\link[BatchMap]{ripple.ord}}
##' @param phase_cores The number of parallel processes to use when estimating
##' the phase of a marker. (Should be no more than 4)
##' @param ripple.cores The number of parallel processes to use when calculating
##' alternative order.
##' @param verbosity A character vector that includes any or all of "batch",
##' "order", "position", "time" and "phase" to output progress status
##' information.
##' @param max.dist The maximum distance (in cM) two markers can have in a batch
##' before automatic reordering is triggered (given that \code{fun.order} was
##' set).
##' @param ws The window size that the reordering function should use
##' @param increase.every Increase the window size by one every n-th round when
##' re-ordering.
##' @param max.tries The maximum number of re-ordering tries. Failing to order
##' after max.tries outputs a warning.
##' @param min.tries The minimum number of re-ordering tries.
##' @param seeds A vector of phase information used as seeds for the first
##' batch
##' @param optimize Either "likelihood" or "count". Passed to \code{ripple.ord}
##' in order to optimize the map's likelihood or the RECORD COUNT criterion.
##' Unless you are absolutely sure why, you should use "likelihood".
##' @param ... Other arguments passed to \code{fun.order}
##' @return A list with the first element \code{Map} being an object of class
##' \code{sequence}, which is a list containing the
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
##' 2-point analyses.} Secondly \code{batches}, a list of \code{Map}s for each
##' of the batches.
##' @author Bastian Schiffthaler, \email{bastian.schiffthaler@umu.se}
##' @seealso \code{\link[BatchMap]{pick_batch_sizes}}, \code{\link[BatchMap]{map}}
##'
##' @keywords utilities
##' @export
map_overlapping_batches <- function(input.seq, size = 50, overlap = 15,
                                    phase_cores = 1, verbosity = NULL, 
                                    ws = 4, seeds = NULL, ...)
{
  #TODO: error checks...
  #Create initial set of batches
  batches <- generate_overlapping_batches(input.seq, size, overlap)
  if("batch" %in% verbosity)
  {
    message("Have ", length(batches), " batches.")
    message("The number of markers in the final batch is: ",
            length(batches[[length(batches)]]))
    message("Processing batch 1...")
  }
  LGs <- list()
  #The first batch is run in full again to get all necessary data (phases etc.)
  if(is.null(seeds))
   {
    LG <- map(input.seq = make_seq(input.seq$twopt, batches[[1]]), phase_cores = phase_cores)
  } else {
     LG <- seeded_map(input.seq = make_seq(input.seq$twopt, batches[[1]],
                               twopt = input.seq$twopt), phase_cores = phase_cores,
                      verbosity = verbosity, seeds = seeds)
   }
  
  round <- 1
  increment <- 0

  LGs[[1]] <- LG
  #Start processing all following batches
  for(i in 2:length(batches))
  {
    if("batch" %in% verbosity) #Print previous batch-map segment
    {
      print(LGs[[i - 1]])
      message("Processing batch ",i,"...")
    }
    #Need to use a seeded map in order to not mess with the overlapping area
    #which we trust more from the previous batch (as that had more information)
    seeds <- tail(LGs[[i - 1]]$seq.phases, overlap)
    batches[[i]][1:(overlap+1)] <- tail(LGs[[i - 1]]$seq.num, overlap + 1)
    LG <- seeded_map(make_seq(input.seq$twopt,
                              batches[[i]],
                              twopt = input.seq$twopt),
                     verbosity = verbosity,
                     seeds = seeds)

    LGs[[i]] <- LG
  }
  #Initialize final order, phases and rfs with the first batch
  final.seq <- LGs[[1]]$seq.num
  final.phase <- LGs[[1]]$seq.phases
  final.rf <- LGs[[1]]$seq.rf
  #Iteratively add data from other batches to the first sequence
  for(i in 2:length(batches))
  {
    start <- length(final.seq) - overlap #start position of overlap with next
    #Add marker order from the next batch to the sequence starting from the
    #start of the overlap
    final.seq[start:length(final.seq)] <- head(LGs[[i]]$seq.num, overlap + 1)
    final.seq <- c(final.seq,
                   LGs[[i]]$seq.num[(overlap + 2):length(LGs[[i]]$seq.num)])
    #Add phases and RFs. We need to shift the indices left by 1
    start <- length(final.phase) - overlap + 1
    final.phase[start:length(final.phase)] <- head(LGs[[i]]$seq.phases, overlap)
    final.phase <- c(final.phase,
                   LGs[[i]]$seq.phases[(overlap + 1):length(LGs[[i]]$seq.phases)])
    final.rf[start:length(final.rf)] <- head(LGs[[i]]$seq.rf, overlap)
    final.rf <- c(final.rf,
                     LGs[[i]]$seq.rf[(overlap + 1):length(LGs[[i]]$seq.rf)])

  }
  if("batch" %in% verbosity)
  {
    message("Final call to map...")
  }
  #Create final sequence and run
  #final.rf is currently only used for debugging purposes
  s <- make_seq(input.seq$twopt, final.seq, final.phase, input.seq$twopt)
  mp <- map(s)
  return(list(Map = mp, batches = LGs))
}
