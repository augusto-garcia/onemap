#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: summary_maps.R                                                #
# Contains: summary_maps                                              #
#                                                                     #
# Written by Jeekin Lau                                               #
# copyright (c) 2023, Jeekin Lau                                      #
#                                                                     #
# First version: 04/24/2023                                           #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
##' Draws a physical vs cM map
##'
##' Provides simple genetic to physical ggplot.
##' @param map.list a map, i.e. an object of class \code{sequence} with a
##' predefined order, linkage phases, recombination fraction and likelihood;
##' also it could be a list of maps.
##' 
##' @param mapping_function either "kosambi" or "haldane"
##' 
##' @return dataframe with basic summary statistics 
##' 
##' @author Jeekin Lau, \email{jeekinlau@@gmail.com}
##' 
##' 
##' @export summary_maps
##' 
##' 
##' 
##' 
summary_maps = function(map.list,mapping_function="kosambi"){
  
  
  if(!(inherits(map.list,c("list", "sequence")))) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequnece'")
  ## if map.list is just a single chormosome, convert it  into a list
  if(inherits(map.list,"sequence")) map.list<-list(map.list)
  
  
  if(mapping_function=="kosambi"){
    summary=data.frame(LG = 1:length(map.list),
                       nMrks = unlist(lapply(map.list, function(x) length(x$seq.num))),
                       map_length = unlist(lapply(map.list, function(x) sum(c(0,kosambi(x$seq.rf))))),
                       max_gap = unlist(lapply(map.list, function(x) kosambi(max(x$seq.rf))))) 
    
    
    last_line=data.frame(LG="All", nMrks=sum(summary$nMrks), map_length=sum(summary$map_length),max_gap=max(summary$max_gap))
    
    stats=rbind(summary,last_line)
  }
  if(mapping_function=="haldane"){
    summary=data.frame(LG = 1:length(map.list),
                       nMrks = unlist(lapply(map.list, function(x) length(x$seq.num))),
                       map_length = unlist(lapply(map.list, function(x) sum(c(0,haldane(x$seq.rf))))),
                       max_gap = unlist(lapply(map.list, function(x) haldane(max(x$seq.rf))))) 
    
    
    last_line=data.frame(LG="All", nMrks=sum(summary$nMrks), map_length=sum(summary$map_length),max_gap=max(summary$max_gap))
    
    stats=rbind(summary,last_line)
  }
  return(stats)
}