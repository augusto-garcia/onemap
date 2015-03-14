#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: defunct.R                                                     #
# Contains: cr3pts, def.rf.3pts, print.rf.3pts                        #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function allocates matrices and calls C routine for three-point analysis
cr3pts <-
function (mrk1, mrk2, mrk3, segr.type1, segr.type2, segr.type3) {
  # checking for correct types of segregation
  if (!any(!is.na(match(c("A.1","A.2","A.3","A.4","B1.5","B2.6","B3.7",
                          "C.8","D1.9","D1.10","D1.11","D1.12","D1.13",
                          "D2.14","D2.15","D2.16","D2.17","D2.18"),
                        segr.type1)))) 
    stop("unknown segregation type for 'mrk1'")

  if (!any(!is.na(match(c("A.1","A.2","A.3","A.4","B1.5","B2.6","B3.7",
                          "C.8","D1.9","D1.10","D1.11","D1.12","D1.13",
                          "D2.14","D2.15","D2.16","D2.17","D2.18"),
                        segr.type2)))) 
    stop("unknown segregation type for 'mrk2'")

  if (!any(!is.na(match(c("A.1","A.2","A.3","A.4","B1.5","B2.6","B3.7",
                          "C.8","D1.9","D1.10","D1.11","D1.12","D1.13",
                          "D2.14","D2.15","D2.16","D2.17","D2.18"),
                        segr.type3)))) 
    stop("unknown segregation type for 'mrk3'")

  # 'I' matrixes and 'm' vectors are defined for the three markers,
  # according to the segregation type
  switch(EXPR=segr.type1,
         'A.1' =, 'A.2' =, 'A.3' =, 'A.4' = {I1 <- m1 <- diag(4)},
         'B1.5' = {I1 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3)},
         'B2.6' = {I1 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3)},
         'B3.7' = {I1 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m1 <- diag(3)},
         'C.8' = {I1 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);
                  m1 <- diag(2)},
         'D1.9' =, 'D1.10' =, 'D1.11' =, 'D1.12' =, 'D1.13' = {I1 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                   m1 <- diag(2)},
		 'D2.14' =, 'D2.15' =, 'D2.16' =, 'D2.17' =, 'D2.18' = {I1 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m1 <- diag(2)})
  
  switch(EXPR=segr.type2,
         'A.1' =, 'A.2' =, 'A.3' =, 'A.4' = {I2 <- m2 <- diag(4)},
         'B1.5' = {I2 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3)},
         'B2.6' = {I2 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3)},
         'B3.7' = {I2 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m2 <- diag(3)},
         'C.8' = {I2 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);
                  m2 <- diag(2)},
         'D1.9' =, 'D1.10' =, 'D1.11' =, 'D1.12' =, 'D1.13' = {I2 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                   m2 <- diag(2)},
		 'D2.14' =, 'D2.15' =, 'D2.16' =, 'D2.17' =, 'D2.18' = {I2 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m2 <- diag(2)})
  
  switch(EXPR=segr.type3,
         'A.1' =, 'A.2' =, 'A.3' =, 'A.4' = {I3 <- m3 <- diag(4)},
         'B1.5' = {I3 <- matrix(c(1,0,0,1,0,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m3 <- diag(3)},
         'B2.6' = {I3 <- matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),4,3,byrow=TRUE);
                   m3 <- diag(3)},
         'B3.7' = {I3 <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),4,3,byrow=TRUE);
                   m3 <- diag(3)},
         'C.8' = {I3 <- matrix(c(1,0,1,0,1,0,0,1),4,2,byrow=TRUE);
                  m3 <- diag(2)},
         'D1.9' =, 'D1.10' =, 'D1.11' =, 'D1.12' =, 'D1.13' = {I3 <- matrix(c(1,0,1,0,0,1,0,1),4,2,byrow=TRUE);
                   m3 <- diag(2)},
		 'D2.14' =, 'D2.15' =, 'D2.16' =, 'D2.17' =, 'D2.18' = {I3 <- matrix(c(1,0,0,1,1,0,0,1),4,2,byrow=TRUE);
                    m3 <- diag(2)})
					
  # procedure to count the number of individuals in each class
  n <- numeric(nrow(m1)*nrow(m2)*nrow(m3))
  k <- 1
  for(p3 in 1:nrow(m3))
    for(p1 in 1:nrow(m1))
      for(p2 in 1:nrow(m2)) {
        n[k] <- length(which(mrk1==p1 & mrk2==p2 & mrk3==p3))
        k <- k+1
      }
  ntot <- sum(n) # total number of individuals

  # calling C routine
  rcmb <- .C("r3pts",
             as.double(t(I1)),
             as.integer(nrow(m1)),
             as.double(t(I2)),
             as.integer(nrow(m2)),
             as.double(I3),
             as.integer(nrow(m3)),
             as.integer(n),
             as.integer(ntot),
             theta12=as.double(numeric(16)),
             theta23=as.double(numeric(16)),
             theta13=as.double(numeric(16)),
             log_like=as.double(numeric(16)),
             posterior=as.double(numeric(16)),
             LOD=as.double(numeric(16)),
             PACKAGE="onemap")
  
  # results
  final <- cbind(rcmb$theta12,rcmb$theta23,rcmb$theta13,rcmb$log_like,
                 rcmb$posterior,rcmb$LOD)
  dimnames(final) <- list(c("A11","A12","A13","A14","A21","A22","A23",
                            "A24","A31","A32","A33","A34","A41","A42",
                            "A43","A44"),list("Theta12","Theta23",
                                              "Theta13","log-Like",
                                              "Posterior","LODs"))
  final
}



# This function is directly called to perform three-point analysis
def.rf.3pts <-
function(w, mrk1name=NULL, mrk2name=NULL, mrk3name=NULL, LOD=5, max.rf=0.35,
         max.nolink=0.55) {
  # checking for correct objects
  if(!any(class(w)=="outcross")) stop(deparse(substitute(w))," is not an object of class 'outcross'")
  if (w$n.mar<3)
    stop("there must be at least three markers to proceed with analysis")
  if (is.null(mrk1name) || is.null(mrk2name) || is.null(mrk3name))
    stop("you must specify 3 markers")
  
  marnames <- colnames(w$geno)

  # checking if markers really exist
  if (is.character(mrk1name) && is.character(mrk2name) && is.character(mrk3name)) {
    mrk1 <- which(marnames==mrk1name)
    mrk2 <- which(marnames==mrk2name)
    mrk3 <- which(marnames==mrk3name)
    if (length(mrk1)==0) stop("marker ", mrk1name, " not found")
    if (length(mrk2)==0) stop("marker ", mrk2name, " not found")
    if (length(mrk3)==0) stop("marker ", mrk3name, " not found")
    if (mrk1==mrk2 || mrk1==mrk3 || mrk2==mrk3)
      stop("the 3 markers must be different")
  }
  else stop("'mrk1name', 'mrk2name' and 'mrk3name' must be of type \"character\"")

  # 'max.nolink' must be greater than 'max.rf'
  if (max.nolink < max.rf)
    stop("'max.rf' must be less than or equal to 'max.nolink'")
  
  goodness <- character(16)
  flag <- NA
  recomb <- c(0.5, 0.5)

  # indirect call to C routine
  test.3pt <- cr3pts(w$geno[,mrk1],w$geno[,mrk2],w$geno[,mrk3],
                     w$segr.type[mrk1],w$segr.type[mrk2],w$segr.type[mrk3])

  # choosing most probable assignments
  probab <- which(test.3pt[,6]>(max(test.3pt[,6]-0.005)) &
                  test.3pt[,6]<=max(test.3pt[,6]))
  for (a in probab) {
    if (test.3pt[a,1] < max.nolink &&
        test.3pt[a,2] < max.nolink &&
        test.3pt[a,3] < max.nolink) {
      if (test.3pt[a,1] < max.rf &&
          test.3pt[a,2] < max.rf &&
          test.3pt[a,3] >= test.3pt[a,1] &&
          test.3pt[a,3] >= test.3pt[a,2])
        goodness[a] <- "****"
      else goodness[a] <- "*"
    }
    else goodness[a] <- "-"
  }
  goodness[-probab] <- "-"
  phase <- which(goodness=="****")

  # no assignment meets 'max.rf' and 'max.nolink' criteria
  if (length(phase)==0) {
    phase <- "none"
  }
  else { # one or more assignments meet both criteria

    # if more than one probable phase, choose the first one
    if (length(phase)>1) {
      phase <- phase[1]
      flag <- 1
    }
    else flag <- 0

    # checking if phase meets 'LOD' criterion
    if (test.3pt[phase,6] >= LOD) {
      recomb <- c(test.3pt[phase,1],test.3pt[phase,2])
      phase <- c("A11","A12","A13","A14","A21","A22","A23","A24","A31",
                 "A32","A33","A34","A41","A42","A43","A44")[phase]
    }
    else { # 'LOD' criterion not met
      phase <- "ns"
      flag <- NA
    }
  }

  # results
  marnames <- marnames[c(mrk1,mrk2,mrk3)]
  dimnames(test.3pt) <- list(c("A11","A12","A13","A14","A21","A22","A23",
                               "A24","A31","A32","A33","A34","A41","A42",
                               "A43","A44"),list("Theta12","Theta23",
                                                 "Theta13","log-Like",
                                                 "Posterior","LODs"))
  rf.3pts <- list(LOD=LOD, max.rf=max.rf, max.nolink=max.nolink,
                  marnames=marnames, recomb=recomb, phase=phase,
                  analysis=test.3pt, goodness=goodness, flag=flag)
  class(rf.3pts) <- "rf.3pts"
  rf.3pts
}



# print method for object class 'rf.3pts'
print.rf.3pts <-
function(x,...) {
  # checking for correct objects
  if(!any(class(x)=="rf.3pts")) stop(deparse(substitute(x))," is not an object of class 'rf.3pts'")

  # printing basic information (criteria)
  cat("  This is an object of class 'rf.3pts'\n")
  cat(paste("  Results of the 3-point analysis for markers: ",
            x$marnames[1], ", ", x$marnames[2], " and ",
            x$marnames[3], "\n",sep=""))
  cat("  Criteria: LOD = ", x$LOD)
  cat("\n             Maximum recombination fraction between adjacent markers = ",
      x$max.rf)
  cat("\n             Maximum recombination fraction between markers on the two ends = ",
      x$max.nolink, "\n\n")

  # results of the three-point analysis
  print(x$analysis)

  if (x$phase=="none") {  # order possibly has problems
    cat("\n  The results do not meet the criteria under any assignment.")
    if (any(x$goodness=="*"))
      cat("\n  The order of the three markers may be wrong.\n")
  }
  else {
    if (x$phase=="ns") {  # LOD smaller than threshold
      cat("\n  The test was not significant (LOD score smaller than specified).\n")
    }
    else {  # results meet all criteria for at least one assignment
      dist <- c(x$recomb,0)
      phases <- character(3)
      if (substr(x$phase,2,2)=="1")
        phases[1] <- "coupling/coupling"
      else if (substr(x$phase,2,2)=="2")
        phases[1] <- "coupling/repulsion"
      else if (substr(x$phase,2,2)=="3")
        phases[1] <- "repulsion/coupling"
      else if (substr(x$phase,2,2)=="4")
        phases[1] <- "repulsion/repulsion"
      
      if (substr(x$phase,3,3)=="1")
        phases[2] <- "coupling/coupling"
      else if (substr(x$phase,3,3)=="2")
        phases[2] <- "coupling/repulsion"
      else if (substr(x$phase,3,3)=="3")
        phases[2] <- "repulsion/coupling"
      else if (substr(x$phase,3,3)=="4")
        phases[2] <- "repulsion/repulsion"

      # final results
      results <- cbind(x$marnames,round(kosambi(dist),3),
                       round(haldane(dist),3),phases)
      results[length(dist),2:3] <- ""
      results <- as.data.frame(results)
      colnames(results) <- c("Markers", "cM Kosambi", "cM Haldane",
                             "Linkage Phases")
      cat("\n  Results under the most probable linkage phases: assignment",
          x$phase, "\n")
      print(results)
      if (!is.na(x$flag) && x$flag==1)
        cat("\n  Warning: more than 1 probable assignment!\n")
    }
  }
}

# end of file
