#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: plot_genome_vs_cm.R                                           #
# Contains: plot_genome_vs_cm                                         #
#                                                                     #
# Written by Jeekin Lau with minor                                    #
#  modifications by Cristiane Taniguti                                #
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
##' @param group.names vector with group name for each sequence object in the map.list
##' 
##' @return ggplot with cM on x-axis and physical position on y-axis
##' 
##' @author Jeekin Lau, \email{jeekinlau@@gmail.com}
##' 
##' @export plot_genome_vs_cm
##' 
##' @import ggplot2 ggpubr
plot_genome_vs_cm = function(map.list,mapping_function="kosambi", group.names=NULL){
  if(!(inherits(map.list,c("list", "sequence")))) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequence'")
  
  ## if map.list is just a single chormosome, convert it  into a list
  if(inherits(map.list,"sequence")) map.list<-list(map.list)
  
  if(is.null(map.list[[1]]$data.name$POS) | is.null(map.list[[1]]$data.name$CHROM)) stop("Reference genome chromosome and position informaion is not available in the dataset.")
  
  
  number_chromomes = length(map.list)
  plot=list()
  for (i in 1:number_chromomes){
    data_for_plot = data.frame(Marker = map.list[[i]]$seq.num, 
                               Chrom=map.list[[i]]$data.name$CHROM[map.list[[i]]$seq.num], 
                               Position = map.list[[i]]$data.name$POS[map.list[[i]]$seq.num],
                               cM = {if(mapping_function=="kosambi") {
                                 cM =cumsum(c(0,kosambi(map.list[[i]]$seq.rf)))  
                               } else {
                                 cM=cumsum(c(0,haldane(map.list[[i]]$seq.rf)))
                               }})
    plot[[i]]=ggplot(data_for_plot,mapping=aes(cM,Position))+geom_point()+ ylab("Genomic position") +
      {if(!is.null(group.names)) ggtitle(paste0(group.names[i]))} + theme_bw() 
  }
  
  a=ggarrange(plotlist=plot)
  return(a)
}