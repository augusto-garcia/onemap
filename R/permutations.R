#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: permutations.R                                                #
# Contains: permt.tot, perm.pars                                      #
#                                                                     #
# Written Marcelo Mollinari                                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#N! combination
perm.tot <- function(v){
  n <- length(v)
  result <- v 
  if (n>1){
    M <- perm.tot(v[2:n]) 
    result <- cbind(v[1],M) 
    if (n>2){
      for (i in 2:(n-1)){
        N <- cbind(M[,1:(i-1)],v[1],M[,i:(n-1)])  
        result <- rbind(result,N) 
      }
    }
    N <- cbind(M,v[1]) 
    result <- rbind(result,N) 
  }
  return(result) 
}

#N!/2 combination
perm.pars <- function(v){
  n <- length(v)
  result <- v 
  if (n>2){
    Mt <- perm.tot(v[2:n]) 
    result <- cbind(v[1],Mt) 
    f <- floor(n/2)
    c <- ceiling(n/2)
    if (n>3){
      for (i in 2:f){
        N <- cbind(Mt[,1:(i-1)],v[1],Mt[,i:(n-1)])
        result <- rbind(result,N) 
      }
    }
    if (c>f){ 
      Ms <- perm.pars(v[2:n]) 
      if (n>3) {N <- cbind(Ms[,1:f],v[1],Ms[,c:(n-1)])} 
      else {N <- cbind(Ms[1:f],v[1],Ms[c:(n-1)])}   
      result <- rbind(result,N) 
    }
  }
  return(result) 
}
