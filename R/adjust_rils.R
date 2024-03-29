#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: adjust_rils.R                                                 #
# Contains: adjust_rf_ril                                             #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# Adapted from qtl package                                            #
# copyright (c) 2001-10, Karl W Broman                                #
#                                                                     #
# First version: Feb 2001                                             #
# License: GNU General Public License version 2 (June, 1991) or later # 
#                                                                     #
#######################################################################

adjust_rf_ril <-
  function(r, type=c("riself", "risib"), expand = TRUE)
{
  ## type of RI lines
  type <- match.arg(type)
  if(type=="riself") {
    if(expand) return(r*2/(1+2*r))
    else return(r/2/(1-r))
  }
  else {
    if(expand) return(r*4/(1+6*r))
    else return(r/(4-6*r))
  }
}
