#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: acum.R                                                        #
# Contains: acum                                                      #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

acum <- function(w) {
  if (w<0) stop("'w' should be equal to or higher than zero")

  # the famous gaussian sum from 1 to w
  w*(w+1)/2
}

# end of file
