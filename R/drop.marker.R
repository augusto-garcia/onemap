#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: drop.marker.R                                                 #
# Contains: drop.marker                                               #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function remove markers from a sequence
drop.marker<-function(input.seq, mrks)
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
    return(make.seq(get(input.seq$twopt),seq.num,twopt=input.seq$twopt))
  }
