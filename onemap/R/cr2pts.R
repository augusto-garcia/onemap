#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: cr2pts.R                                                      #
# Contains: cr2pts                                                    #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function allocates matrices and calls C routine for two-point analysis
cr2pts <- 
function (mrk1, mrk2, segr.type1, segr.type2) {
  # checking for correct types of segregation
  if (!any(!is.na(match(1:7,segr.type1)))) 
    stop("unknown segregation type for 'mrk1'")
  
  if (!any(!is.na(match(1:7,segr.type2)))) 
    stop("unknown segregation type for 'mrk2'")

  # 'I' matrices are defined for both markers, according to the segregation type
  switch(EXPR=segr.type1,
         {I1 <- diag(4)},
         {I1 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE)},
         {I1 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE)},
         {I1 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE)},
         {I1 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE)},
         {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE)},
         {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE)}
		)
  
  switch(EXPR=segr.type2,
         {I2 <- diag(4)},
         {I2 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE)},
         {I2 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE)},
         {I2 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE)},
         {I2 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE)},
         {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE)},
         {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE)}
		)
  
  # procedure to count the number of individuals in each class
  n <- numeric(ncol(I1)*ncol(I2))
  k <- 1
  for(p2 in 1:ncol(I2))
    for(p1 in 1:ncol(I1)) {
      n[k] <- length(which(mrk1==p1 & mrk2==p2))
      k <- k+1
    }
  ntot <- sum(n) # total number of individuals

  # calling C routine
  rcmb <- .C("r2pts",
             as.double(t(I1)),
             as.integer(ncol(I1)),
             as.double(I2),
             as.integer(ncol(I2)),
             as.integer(n),
             as.integer(ntot),
             r=as.double(numeric(4)),
             log_like=as.double(numeric(4)),
             posterior=as.double(numeric(4)),
             LOD=as.double(numeric(4)),
             PACKAGE="onemap")
  
  if ((segr.type1 == 6 && segr.type2 == 7) || (segr.type1 == 7 && segr.type2 == 6))
    rcmb$r <- rep(0.5, 4)
  
  if (all(n == 0)) {
    rcmb$r <- rep(0.5, 4)
    rcmb$log_like <- rcmb$LOD <- rep(0.0, 4)
    rcmb$posterior <- rep(0.25, 4)
  }
  
  # results
  final <- cbind(rcmb$r,rcmb$log_like,rcmb$posterior,rcmb$LOD)
  dimnames(final) <- list(c("1","2","3","4"),list("Theta","log-Like","Posterior","LODs"))
  final
}

# end of file
