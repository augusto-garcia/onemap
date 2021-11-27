#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: suggest_lod.R                                                 #
# Contains: suggest_lod                                               #
#                                                                     #
# Written by Antonio Augusto Franco Garcia and updated by Cristiane   #
# Taniguti                                                            #
# copyright (c) 2015 Antonio Augusto Franco Garcia                    #
#                                                                     #
# First version: 2015/04/21                                           #
# License: GNU General Public License version 3 or later              #
#                                                                     #
#######################################################################

##' Suggests a LOD Score for two point tests
##'
##' It suggests a LOD Score for declaring statistical significance for two-point tests
##' for linkage between all pairs of markers, considering that multiple tests are being performed.
##'
##' In a somehow naive approach, the function calculates the number of two-point tests that
##' will be performed for all markers in the data set, and then using this to calculate
##' the global alpha required to control type I error using Bonferroni's correction.
##'
##' From this global alpha, the corresponding quantile from the chi-square distribution is taken
##' and then converted to LOD Score.
##'
##' This can be seen as just an initial approximation to help users to select a LOD Score for two
##' point tests.
##'
##' @param x an object of class \code{sequence} or \code{onemap}
##'
##' @return the suggested LOD to be used for testing linkage
##'
##' @examples
##' 
##' data(onemap_example_bc) # Loads a fake backcross dataset installed with onemap
##' suggest_lod(onemap_example_bc) # An value that should be used to start the analysis
##' 
##' @export
suggest_lod <- function(x) {
    if (is(x,c("sequence", "onemap"))) { # Keep onemap class just to be compatible with older versions
        if(is(x, "onemap"))
            num.tests <- choose(x$n.mar, 2) #Number of pairwise tests
        if(is(x, "sequence"))
            num.tests <- choose(length(x$seq.num), 2)
        LOD <- 0.2172 * qchisq(1-0.05/num.tests, 1) #Corresponding LOD
        return(LOD)
    }
    else stop("This is not a onemap object with raw data")
}
##'
