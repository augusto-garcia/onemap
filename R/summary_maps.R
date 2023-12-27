#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: summary_maps.R                                                #
# Contains: summary_maps_onemap                                       #
#                                                                     #
# Written by Jeekin Lau  with minor modifications                     #
#   from Cristiane Taniguti                                           #
# copyright (c) 2023, Jeekin Lau                                      #
#                                                                     #
# First version: 04/24/2023                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Create table with summary information about the linkage map
##'
##' @param map.list a map, i.e. an object of class \code{sequence} with a
##' predefined order, linkage phases, recombination fraction and likelihood;
##' also it could be a list of maps.
##' 
##' @param mapping_function either "kosambi" or "haldane"
##' 
##' @return data.frame with basic summary statistics 
##' 
##' @author Jeekin Lau, \email{jeekinlau@@gmail.com}
##' 
##' @import tidyr
##' 
##' @export summary_maps_onemap 
##' 
summary_maps_onemap  = function(map.list, mapping_function="kosambi"){
  
  if(!(inherits(map.list,c("list", "sequence")))) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequence'")
  
  ## if map.list is just a single chromosome, convert it  into a list
  if(inherits(map.list,"sequence")) map.list<-list(map.list)
  
  mk_types <- lapply(map.list, function(x) as.data.frame(table(x$data.name$segr.type[x$seq.num])))
  for(i in 1:length(mk_types)) mk_types[[i]] <- cbind(Var2=i, mk_types[[i]])
  
  mk_types <- do.call(rbind, mk_types)
  mk_types <- pivot_wider(mk_types, names_from =  "Var1", values_from = "Freq")
  mk_types[is.na(mk_types)] <- 0
  
  summary=data.frame(LG = 1:length(map.list),
                     n_mks = unlist(lapply(map.list, function(x) length(x$seq.num))),
                     mk_types[,-1],
                     map_length = unlist(lapply(map.list, function(x) sum(c(0,kosambi(x$seq.rf))))),
                     max_gap = {if(mapping_function=="kosambi") {
                       unlist(lapply(map.list, function(x) kosambi(max(x$seq.rf))))
                     } else  unlist(lapply(map.list, function(x) haldane(max(x$seq.rf))))
                     })
  
  last_line= apply(summary, 2, sum) 
  
  stats=rbind(summary,last_line)
  stats[dim(stats)[1],1] <- "Total"
  
  return(stats)
}