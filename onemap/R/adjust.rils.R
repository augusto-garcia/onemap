#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: adjust.rils.R                                                 #
# Contains: adjust.rf.ril                                             #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# Adapted from read.cross.mm (found in the R package qtl)             #
# copyright (c) 2001-10, Karl W Broman                                #
#                                                                     #
# First version: Feb 2001                                             #
# Last update:   Oct 2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later # 
#                                                                     #
#######################################################################

adjust.rf.ril <-
  function(r, type=c("riself.onemap","risib.onemap"), expand=TRUE)
{
  ## type of RI lines
  type <- match.arg(type)
  if(type=="riself.onemap") {
    if(expand) return(r*2/(1+2*r))
    else return(r/2/(1-r))
  }
  else {
    if(expand) return(r*4/(1+6*r))
    else return(r/(4-6*r))
  }
}
