#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: make.seq.R                                                    #
# Contains: make.seq, print.sequence                                  #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to create sequences based on other object types
make.seq <- 
function(input.obj, arg=NULL, phase=NULL, twopt=NULL) {
  # checking for correct object
  if(all(is.na(match(class(input.obj),c("rf.2pts","group","compare","try","order")))))
    stop(deparse(substitute(input.obj))," is not an object of classes 'rf.2pts', 'group', 'compare', 'try' or 'order'")

  switch(EXPR=class(input.obj),
         'rf.2pts' = {
           if (length(arg) == 1 && arg == "all") seq.num <- 1:input.obj$n.mar # generally used for grouping markers
		   else if(is.vector(arg) && is.numeric(arg)) seq.num <- arg
		   else stop("for an object of class 'rf.2pts', \"arg\" must be a vector of integers or the string 'all'")
		   ### CHECK IF MARKERS REALLY EXIST
           if (is.null(phase)) seq.phases <- -1 # no predefined linkage phases
           else if(length(phase) == (length(seq.num)-1)) seq.phases <- phase
           else stop("the length of 'phase' must be equal to the length of the sequence minus 1")
           seq.rf <- -1
           seq.like <- NULL
           if(is.null(twopt)) twopt <- deparse(substitute(input.obj))
         },
         'group' = {
           if(length(arg) == 1 && is.numeric(arg) && arg <= input.obj$n.groups) seq.num <- input.obj$seq.num[which(input.obj$groups == arg)]
		   else stop("for this object of class 'group', \"arg\" must be an integer less than or equal to ",input.obj$n.groups)
           seq.phases <- -1
           seq.rf <- -1
           seq.like <- NULL
           twopt <- input.obj$twopt
         },
		 'compare' = {
           n.ord <- max(which(head(input.obj$best.ord.LOD,-1) != -Inf))
           unique.orders <- unique(input.obj$best.ord[1:n.ord,])
		   if(is.null(arg)) seq.num <- unique.orders[1,] # NULL = 1 is the best order
           else if(length(arg) == 1 && is.numeric(arg) && arg <= nrow(unique.orders)) seq.num <- unique.orders[arg,]
		   else stop("for this object of class 'compare', \"arg\" must be an integer less than or equal to ",nrow(unique.orders))
           if (is.null(phase)) phase <- 1 # NULL = 1 is the best combination of phases
           chosen <- which(apply(input.obj$best.ord[1:n.ord,],1,function(x) all(x==seq.num)))[phase]
           seq.phases <- input.obj$best.ord.phase[chosen,]
           seq.rf <- input.obj$best.ord.rf[chosen,]
           seq.like <- input.obj$best.ord.like[chosen]
           twopt <- input.obj$twopt
         },
		 'try' = {
		   if(length(arg) != 1 || !is.numeric(arg) || arg > length(input.obj$ord))
		     stop("for this object of class 'try', \"arg\" must be an integer less than or equal to ",length(input.obj$ord))
           if (is.null(phase)) phase <- 1 # NULL = 1 is the best combination of phases
           seq.num <- input.obj$try.ord[arg,]
           seq.phases <- input.obj$ord[[arg]]$phase[phase,]
           seq.rf <- input.obj$ord[[arg]]$rf[phase,]
           seq.like <- input.obj$ord[[arg]]$like[phase]
           twopt <- input.obj$twopt
         },
		 'order' = {
           arg <- match.arg(arg,c("safe","force"))
           if (arg == "safe") {
		     # order with safely mapped markers
             seq.num <- input.obj$ord$seq.num
             seq.phases <- input.obj$ord$seq.phases
             seq.rf <- input.obj$ord$seq.rf
             seq.like <- input.obj$ord$seq.like
           }
           else {
		     # order with all markers
             seq.num <- input.obj$ord.all$seq.num
             seq.phases <- input.obj$ord.all$seq.phases
             seq.rf <- input.obj$ord.all$seq.rf
             seq.like <- input.obj$ord.all$seq.like
           }
           twopt <- input.obj$twopt
         }
		 )

  # check if any marker appears more than once in the sequence
  if(length(seq.num) != length(unique(seq.num))) stop("there are duplicated markers in the sequence")
  
  structure(list(seq.num=seq.num, seq.phases=seq.phases, seq.rf=seq.rf, seq.like=seq.like,
                 data.name=input.obj$data.name, twopt=twopt), class = "sequence")
}

# print method for object class 'sequence'
print.sequence <- function(x,...) {
  marnames <- colnames(get(x$data.name, pos=1)$geno)[x$seq.num]
  if(is.null(x$seq.like)) {
    # no information available for the order
    cat("\nNumber of markers:",length(marnames))
    cat("\nMarkers in the sequence:\n")
    cat(marnames,fill=TRUE)
    cat("\nParameters not estimated.\n\n")
  }
  else {
    # convert numerical linkage phases to strings
    link.phases <- matrix(NA,length(x$seq.num),2)
    link.phases[1,] <- rep(1,2)
    for (i in 1:length(x$seq.phases)) {
      switch(EXPR=x$seq.phases[i],
             link.phases[i+1,] <- link.phases[i,]*c(1,1),
             link.phases[i+1,] <- link.phases[i,]*c(1,-1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
             )
    }

    ## display results
    longest.name <- max(nchar(marnames))
    marnames <- formatC(marnames,flag="-")
    longest.number <- max(nchar(x$seq.num))
    marnumbers <- formatC(x$seq.num, format="d", width=longest.number)   
    distances <- formatC(c(0,cumsum(get(get(".map.fun", envir=.onemapEnv))(x$seq.rf))),format="f",digits=2,width=7)
    ## whith diplotypes for class 'outcross'
    if(class(get(x$data.name, pos=1))=="outcross"){
      ## create diplotypes from segregation types and linkage phases
      link.phases <- apply(link.phases,1,function(x) paste(as.character(x),collapse="."))
      parents <- matrix("",length(x$seq.num),4)
      for (i in 1:length(x$seq.num))
        parents[i,] <- return.geno(get(x$data.name, pos=1)$segr.type[x$seq.num[i]],link.phases[i])
      cat("\nPrinting map:\n\n")
      cat("Markers",rep("",max(longest.number+longest.name-7,0)+10),"Position",rep("",10),"Parent 1","     ","Parent 2\n\n")
      for (i in 1:length(x$seq.num)) {
        cat(marnumbers[i],marnames[i],rep("",max(7-longest.name-longest.number,0)+10),distances[i],rep("",10),parents[i,1],"|  |",parents[i,2],"     ",parents[i,3],"|  |",parents[i,4],"\n")
      }
      cat("\n")
      cat(length(marnames),"markers            log-likelihood:",x$seq.like,"\n\n")
    }
    ## whithout diplotypes for another classes
    else if(class(get(x$data.name, pos=1))=="bc.onemap" || class(get(x$data.name, pos=1))=="f2.onemap" || class(get(x$data.name, pos=1))=="riself.onemap" || class(get(x$data.name, pos=1))=="risib.onemap"){
      cat("\nPrinting map:\n\n")
      cat("Markers",rep("",max(longest.number+longest.name-7,0)+10),"Position",rep("",10),"\n\n")
      for (i in 1:length(x$seq.num)) {
        cat(marnumbers[i],marnames[i],rep("",max(7-longest.name-longest.number,0)+10),distances[i],rep("",10),"\n")
      }
      cat("\n",length(marnames),"markers            log-likelihood:",x$seq.like,"\n\n")
    }
    else warning("invalid type of cross") 
  }
}
##end of file
