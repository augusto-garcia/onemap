#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: est.map.c.R                                                   #
# Contains: est.map.c                                                 #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# Adapted from est.map.c (package: R/qtl)                             #
# copyright (c) 2000-6, Karl W Broman                                 #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to estimate the parameters of a linkage map (recombination
# frequencies and log-likelihood)
est.map.c <-
function(geno,type,phase,rec,verbose,tol=1E-2) {
  error.prob <- 1E-50
  maxit <- 1000 # maximum number of iterations
  z1 <- .C("est_map_outbred",
           as.integer(nrow(geno)),    # number of individuals
           as.integer(ncol(geno)),    # number of markers
           as.integer(type),
           as.integer(phase),
           as.integer(geno),          # genotype data
           rf=as.double(rec),         # initial recombination fractions
           as.double(error.prob),
           loglike=as.double(0),      # log-likelihood
           as.integer(maxit),
           as.double(tol),
           as.integer(verbose),
           PACKAGE="onemap")
  list(rf=z1$rf, loglike=z1$loglike)
}

# end of file
