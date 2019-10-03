#######################################################################
##                                                                     #
## Package: onemap                                                     #
##                                                                     #
## File: drop_marker.R                                                 #
## Contains: drop_marker                                               #
##                                                                     #
## Written by Marcelo Mollinari                                        #
## copyright (c) 2009, Marcelo Mollinari                               #
##                                                                     #
## First version: 02/27/2009                                           #
## Last update: 01/02/2009                                             #
## License: GNU General Public License version 2 (June, 1991) or later #
##                                                                     #
#######################################################################

##'   Creates a new sequence by dropping markers.
##'
##'   Creates a new sequence by dropping markers from a predetermined
##'   one.
##'
##' @param input.seq an object of class \code{sequence}.
##'
##' @param mrks a vector containing the markers to be removed
##'     from the \code{sequence}.
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
##' @seealso \code{\link[onemap]{add_marker}}
##'
##' @examples
##' data(onemap_example_out)
##' twopt <- rf_2pts(onemap_example_out)
##' all_mark <- make_seq(twopt,"all")
##' groups <- group(all_mark)
##' (LG1 <- make_seq(groups,1))
##' (LG.aug<-drop_marker(LG1, c(10,14)))
##'
##' @export

drop_marker<-function(input.seq, mrks)
  {
    if (!any(class(input.seq) == "sequence"))
      stop(deparse(substitute(input.seq)), " is not an object of class 'sequence'")
    seq.num<-input.seq$seq.num;count<-numeric()
    for(i in 1:length(mrks)){
      if(all(seq.num!=mrks[i])) count<-c(count, i)
      seq.num<-seq.num[seq.num!=mrks[i]]
    }
    if(length(count)> 0){
      mrklist<-paste(mrks[count], collapse = ", ")
      msg <- sprintf(ngettext(length(count),
                              "marker %s was not in the sequence",
                              "markers %s were not in the sequence"), mrklist)
      warning(msg, domain=NA)
      mrks<-mrks[-count]
    }
    return(make_seq(input.seq$twopt,seq.num,twopt=input.seq$twopt))
  }
