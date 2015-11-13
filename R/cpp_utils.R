#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: cpp_utils.R                                                   #
# Contains: get.bins                                                  #
# These functions are for internal use only                           #
#                                                                     #
# Written Marcelo Mollinari                                           #
#                                                                     #
# First version: 09/2015                                              #
# Last update: 11/2015                                                #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

get.bins <- function(geno, exact=TRUE)
{
  bins<-.Call("get_bins",
              geno,
              as.numeric(exact),
              options()$width-6,
              PACKAGE = "onemap" )
  return(bins)
}
