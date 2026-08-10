#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: cpp_utils.R                                                   #
# Contains: get_bins est_rf_out est_rf_bc est_map_hmm_f2              #
#  est_map_hmm_bc est_map_hmm_out                                     #
# These functions are for internal use only                           #
#                                                                     #
# Written Marcelo Mollinari                                           #
#                                                                     #
# First version: 09/2015                                              #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function calls C++ routine to find markers with redundant information
##' @useDynLib onemap
##' @import Rcpp
get_bins <- function(geno, exact=TRUE)
{
  bins<-.Call("get_bins",
              geno,
              as.numeric(exact),
              PACKAGE = "onemap" )
  return(bins)
}

# This function calls C++ routine for two-point analysis (outcross)
##' @useDynLib onemap
##' @import Rcpp
est_rf_out<-function(geno, mrk=0, seg_type=NULL, nind, verbose=TRUE)
{
  r<-.Call("est_rf_out_wrap",
           geno,
           mrk=mrk-1,
           as.numeric(seg_type),
           as.numeric(nind),
           as.numeric(verbose),
           PACKAGE = "onemap" )
  
  if(mrk <= 0)
  {
    names(r)<-c("CC", "CR", "RC", "RR")
    for(i in 1:4) dimnames(r[[i]])<-list(colnames(geno), colnames(geno))
    
    # Bug fix with D1D2 - It can be estimated than receive 0 for LOD and 0.25 for rf
    for(i in 1:length(seg_type))
      for(j in 1:(length(seg_type)-1))
        if((seg_type[i] == 7 & seg_type[j] == 6) | (seg_type[i] == 6 & seg_type[j] == 7)){
          r[[1]][i,j] <- r[[2]][i,j] <- r[[3]][i,j] <- r[[4]][i,j] <- 0.25
          r[[1]][j,i] <- r[[2]][j,i] <- r[[3]][j,i] <- r[[4]][j,i] <- 0
        }
    
    # Bug fix: If rf is very close to 0.5, LOD can be very close to zero, but with negative value
    # This is causing numerical problems in further analysis
    r[[1]][which(r[[1]] < 0)] <- 10^(-5)
    r[[2]][which(r[[2]] < 0)] <- 10^(-5)
    r[[3]][which(r[[3]] < 0)] <- 10^(-5)
    r[[4]][which(r[[4]] < 0)] <- 10^(-5)
    
    return(r)
  }
  else
  {
    rownames(r[[1]])<-c("rCC", "rCR", "rRC", "rRR")
    colnames(r[[1]])<-colnames(geno)
    rownames(r[[2]])<-c("LODCC", "LODCR", "LODRC", "LODRR")
    colnames(r[[2]])<-colnames(geno)
    return(r)
  }
}

# This function calls C++ routine for two-point analysis (bc)
##' @useDynLib onemap
##' @import Rcpp
est_rf_bc<-function(geno, mrk=0,  nind, type=0, verbose=TRUE)
{
  r<-.Call("est_rf_bc_wrap",
           geno,
           mrk-1,
           as.numeric(nind),
           as.numeric(type), #0=bc; 1=riself; 2=risib
           as.numeric(verbose),
           PACKAGE = "onemap" )
  if(mrk <= 0)
    dimnames(r)<-list(colnames(geno), colnames(geno))
  else
    dimnames(r)<-list(c("rf", "LOD"), colnames(geno))
  return(r)
}

# This function calls C++ routine for multipoint analysis (f2)
##' @useDynLib onemap
##' @import Rcpp
est_map_hmm_f2<-function(geno, error, rf.vec=NULL, verbose=TRUE, tol=1e-6)
{
  if(length(rf.vec) != (nrow(geno)-1))
    rf.vec = rep(0.1, (nrow(geno)-1))
  r<-.Call("est_hmm_f2",
           geno,
           error,
           as.numeric(rf.vec),
           as.numeric(verbose),
           as.numeric(tol),
           PACKAGE = "onemap" )
  names(r)<-c("rf", "loglike", "probs")
  return(r)
}

# This function calls C++ routine for multipoint analysis (bc)
##' @useDynLib onemap
##' @import Rcpp
est_map_hmm_bc<-function(geno, error, rf.vec=NULL, verbose=TRUE, tol=1e-6)
{
  if(length(rf.vec) != (nrow(geno)-1))
    rf.vec = rep(0.1, (nrow(geno)-1))
  r<-.Call("est_hmm_bc",
           geno,
           error,
           as.numeric(rf.vec),
           as.numeric(verbose),
           as.numeric(tol),
           PACKAGE = "onemap" )
  names(r)<-c("rf", "loglike", "probs")
  return(r)
}

#


##' C++ routine for multipoint analysis in outcrossing populations
##'
##' It calls C++ routine that implements the methodology of Hidden
##' Markov Models (HMM) to construct multipoint linkage maps in
##' outcrossing species
##'
##' @param geno matrix of genotypes. Rows represent marker and columns
##'     represent individuals.
##'
##' @param type a vector indicating the type of marker. For more
##'     information see \code{\link[onemap]{read_onemap}}
##'
##' @param phase a vector indicating the linkage phases between
##'     markers. For more information see
##'     \code{\link[onemap]{make_seq}}
##'
##' @param rf.vec a vector containing the recombination fraction
##'     initial values
##'
##' @param verbose If \code{TRUE}, print tracing information.
##'
##' @param tol tolerance for the C routine, i.e., the value used to
##'     evaluate convergence.
##'
##' @return a list containing the re-estimated vector of recombination
##'     fractions and the logarithm of the likelihood
##'
##' @keywords internal
##'
##' @useDynLib onemap
##' @import Rcpp
##' 
##' @export
##' 
est_map_hmm_out<-function(geno, error, type,  phase, rf.vec=NULL, verbose=TRUE, tol=1e-6)
{
  if(is.null(rf.vec))
    rf.vec<-rep(0.1, (nrow(geno)-1))
  r<-.Call("est_hmm_out",
           geno,
           error,
           as.numeric(type),
           as.numeric(phase),
           as.numeric(rf.vec),
           as.numeric(verbose),
           as.numeric(tol), PACKAGE = "onemap")
  names(r)<-c("rf", "loglike", "probs")
  return(r)
}
#end of the file