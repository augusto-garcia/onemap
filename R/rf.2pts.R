#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: rf.2pts.R                                                     #
# Contains: rf.2pts, print.rf.2pts                                    #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## Function to perform two-point analyses for all markers in a data set
rf.2pts <- 
function(input.obj, LOD=3, max.rf=0.50, verbose = TRUE) {
  ## checking for correct object
  if(!any(class(input.obj)=="outcross"||class(input.obj)=="f2.onemap" || class(input.obj)=="bc.onemap" || class(input.obj)=="riself.onemap" || class(input.obj)=="risib.onemap")) stop(deparse(substitute(input.obj))," is not an object of class 'outcross', 'bc.onemap', 'f2.onemap', 'riself.onemap' or 'risib.onemap'")
  if (input.obj$n.mar<2) stop("there must be at least two markers to proceed with analysis")

  ## creating variables (result storage and progress output)
  count <- 0
  prtd <- 0
  progr.past <- 0
  tot <- choose(input.obj$n.mar,2)
  analysis <- array(dim=c(tot,4,2))
  
  if (verbose==TRUE) {
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)
##    cat("--Progress: 0%")
    
    for (i in 2:input.obj$n.mar) {
      for (j in 1:(i-1)) {
        count <- count+1

        ## indirect call to C routine
        if ((input.obj$segr.type.num[i] == 6 && input.obj$segr.type.num[j] == 7) ||
            (input.obj$segr.type.num[i] == 7 && input.obj$segr.type.num[j] == 6)) {
          analysis[acum(i-2)+j,,] <- matrix(rep(c(0.25,0), each = 4), nrow = 4)
        }
        
        else{
          current <- cr2pts(input.obj$geno[,i],input.obj$geno[,j],input.obj$segr.type.num[i],input.obj$segr.type.num[j])
          analysis[acum(i-2)+j,,] <- current[,c(1,4)]
        }

        setTxtProgressBar(pb, count/tot)
        
        ## current progress output
##        if (count==tot) cat("\n--Finished\n")
##        else {
##          progr <- round(100*(count/tot),0)
##          if ((progr - progr.past) >= 5) {
##            if (prtd==10) {
##              cat("\n            ")
##              prtd <- 0
##            }
##            cat(paste("....",progr,"%",sep=""))
##            flush.console()
##            progr.past <- progr
##            prtd <- prtd+1
##          }
##        }
      }
    }
    close(pb)
  }
  else {
    ## two-point analysis for each pair of markers
    for (i in 2:input.obj$n.mar) {
      for (j in 1:(i-1)) {
        count <- count+1

        ## indirect call to C routine
        if ((input.obj$segr.type.num[i] == 6 && input.obj$segr.type.num[j] == 7) ||
            (input.obj$segr.type.num[i] == 7 && input.obj$segr.type.num[j] == 6)) {          
          analysis[acum(i-2)+j,,] <- matrix(rep(c(0.25,0), each = 4), nrow = 4)
        } else{
          current <- cr2pts(input.obj$geno[,i],input.obj$geno[,j],input.obj$segr.type.num[i],input.obj$segr.type.num[j])
          analysis[acum(i-2)+j,,] <- current[,c(1,4)]
        }
      }
    }
  }
  ## results
  ## If type of cross = f2, bc riself or risib, discard inmpossible phases
  if(class(input.obj)=="f2.onemap" || class(input.obj)=="bc.onemap"){
    for (i in 2:input.obj$n.mar) {
      for (j in 1:(i - 1)) {
        if(abs(sum(input.obj$phase[c(i,j)]))==2)
          analysis[acum(i - 2) + j, , ]<-c(analysis[acum(i - 2) + j, , ][1],.5,.5,.5,analysis[acum(i - 2) + j, , ][5],0,0,0)
        else if(abs(sum(input.obj$phase[c(i,j)]))==0)
          analysis[acum(i - 2) + j, , ]<-c(.5,.5,.5,analysis[acum(i - 2) + j, , ][4],0,0,0,analysis[acum(i - 2) + j, , ][8])
        else ("Shoud not get here. Linkage phase error.")
      }
    }
    analysis[,,1][analysis[,,1] > 0.5] <- 0.5
  }
  else if(class(input.obj)=="riself.onemap" || class(input.obj)=="risib.onemap"){
    for (i in 2:input.obj$n.mar) {
      for (j in 1:(i - 1)) {
        if(abs(sum(input.obj$phase[c(i,j)]))==2)
          analysis[acum(i - 2) + j, , ]<-c(adjust.rf.ril(analysis[acum(i - 2) + j, , ][1], type=class(input.obj),
                                                        expand = FALSE),.5,.5,.5,analysis[acum(i - 2) + j, , ][5],0,0,0)
        else if(abs(sum(input.obj$phase[c(i,j)]))==0)
          analysis[acum(i - 2) + j, , ]<-c(.5,.5,.5,adjust.rf.ril(analysis[acum(i - 2) + j, , ][4], type=class(input.obj),
                                                        expand = FALSE),0,0,0,analysis[acum(i - 2) + j, , ][8])
        else ("Shoud not get here. Linkage phase error.")
      }
    }
    analysis[,,1][analysis[,,1] > 0.5] <- 0.5
  }
  dimnames(analysis) <- list(NULL,c("1","2","3","4"),
                             c("Theta","LODs"))
  structure(list(data.name=as.character(sys.call())[2], n.mar=input.obj$n.mar,marnames=colnames(input.obj$geno),
                 LOD=LOD, max.rf=max.rf, input=input.obj$input, analysis=analysis), class = "rf.2pts")
}



## print method for object class 'rf.2pts'
print.rf.2pts <-
function(x, mrk1=NULL, mrk2=NULL,...) {
  ## checking for correct object
  if(!any(class(x)=="rf.2pts")) stop(deparse(substitute(x))," is not an object of class 'rf.2pts'")

  if (is.null(mrk1) || is.null(mrk2)) {
    ## printing a brief summary
    cat("  This is an object of class 'rf.2pts'\n")
    cat("\n  Criteria: LOD =", x$LOD, ", Maximum recombination fraction =",
        x$max.rf, "\n")
    cat("\n  This object is too complex to print\n")
    cat("  Type 'print(object,mrk1=marker,mrk2=marker)' to see the analysis for two markers\n")
    cat("    mrk1 and mrk2 can be the names or numbers of both markers\n")
  }
  else {
    ## printing detailed results for two markers
	
    ## checking if markers exist and converting character to numeric
    if (is.character(mrk1) && is.character(mrk2)) {
      mrk1name <- mrk1
      mrk2name <- mrk2
      mrk1 <- which(x$marnames==mrk1)
      mrk2 <- which(x$marnames==mrk2)
      if (length(mrk1)==0) stop("marker ", mrk1name, " not found")
      if (length(mrk2)==0) stop("marker ", mrk2name, " not found")
    }
    if (is.numeric(mrk1) && is.numeric(mrk2)) {
      cat("  Results of the 2-point analysis for markers:", x$marnames[mrk1],
          "and", x$marnames[mrk2], "\n")

        ## results found
        cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
            x$max.rf, "\n\n")
        ## do not print anything for the same marker 
        if(mrk1 == mrk2) stop("mrk1 and mrk2 are the same")
        ## results of the two-point analysis
      if (mrk1 > mrk2){
        output<-x$analysis[acum(mrk1-2)+mrk2,,]
        ## checking data type
        if(class(get(x$data.name, pos=1))=="outcross") print(output)
        else{
          colnames(output)<-c("rf","LOD")
          print(output[output[,2]!=0,])
        }
      }
      else{
        output<-x$analysis[acum(mrk2-2)+mrk1,,]
        ## checking data type
        if(class(get(x$data.name, pos=1))=="outcross") print(output)
        else{
          colnames(output)<-c("rf","LOD")
          print(output[output[,2]!=0,])
        }
      }
    }
    else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
  }
}

## end of file


