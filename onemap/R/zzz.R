#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# zzz.R                                                               #
# Contains: .First.lib                                                #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 03/12/2012                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

.First.lib <- function(lib, pkg) library.dynam("onemap", pkg, lib)
.onemapEnv <- new.env()
assign(".map.fun",  "kosambi", envir = .onemapEnv)
# end of file
