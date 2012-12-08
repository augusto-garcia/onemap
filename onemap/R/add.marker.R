#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: add.marker.R                                                  #
# Contains: add.marker                                                #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function add markers to a sequence
add.marker<-function(input.seq, mrks)
  {
    if (!any(class(input.seq) == "sequence")) 
      stop(sQuote(deparse(substitute(input.seq))), " is not an object of class 'sequence'")
    seq.num<-c(input.seq$seq.num,mrks)
    return(make.seq(get(input.seq$twopt),seq.num, twopt=input.seq$twopt))
  }
