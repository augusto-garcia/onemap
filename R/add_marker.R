#######################################################################
##                                                                     #
## Package: onemap                                                     #
##                                                                     #
## File: add_marker.R                                                  #
## Contains: add_marker                                                #
##                                                                     #
## Written by Marcelo Mollinari                                        #
## copyright (c) 2009, Marcelo Mollinari                               #
##                                                                     #
## First version: 02/27/2009                                           #
## Last update: 01/02/2016                                             #
## License: GNU General Public License version 2 (June, 1991) or later #
##                                                                     #
#######################################################################

##'   Creates a new sequence by adding markers.
##'
##'   Creates a new sequence by adding markers from a predetermined
##'   one. The markers are added in the end of the sequence.
##'
##' @param input.seq an object of class \code{sequence}.
##'
##' @param mrks a vector containing the markers to be added from the
##'     \code{sequence}.
##'
##' @return An object of class \code{sequence}, which is a list
##'     containing the following components:
##'
##' \item{seq.num}{a \code{vector} containing the (ordered) indices of
##'     markers in the sequence, according to the input file.}
##'
##' \item{seq.phases}{a \code{vector} with the linkage phases between
##'     markers in the sequence, in corresponding positions. \code{-1}
##'     means that there are no defined linkage phases.}
##'
##' \item{seq.rf}{a \code{vector} with the recombination fractions
##'     between markers in the sequence. \code{-1} means that there
##'     are no estimated recombination fractions.}
##'
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##'     \item{data.name}{name of the object of class \code{onemap}
##'     with the raw data.}
##'
##' \item{twopt}{name of the object of class \code{rf_2pts} with the
##'     2-point analyses.}
##'
##'  @author Marcelo Mollinari, \email{mmollina@@usp.br}
##'
##' @seealso \code{\link[onemap]{drop_marker}}
##'
##' @examples
##' data(example_out)
##' twopt <- rf_2pts(example_out)
##' all_mark <- make_seq(twopt,"all")
##' groups <- group(all_mark)
##' (LG1 <- make_seq(groups,1))
##' (LG.aug<-add_marker(LG1, c(4,7)))
##'
##' @export
add_marker<-function(input.seq, mrks)
  {
    if (!any(class(input.seq) == "sequence"))
      stop(sQuote(deparse(substitute(input.seq))), " is not an object of class 'sequence'")
    seq.num<-c(input.seq$seq.num,mrks)
    return(make_seq(get(input.seq$twopt),seq.num, twopt=input.seq$twopt))
  }
