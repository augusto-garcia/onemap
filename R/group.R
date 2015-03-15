#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: group.R                                                       #
# Contains: group, print.group                                        #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to assign markers to linkage groups


##' Assign markers to linkage groups
##' 
##' Identifies linkage groups of markers, using results from two-point
##' (pairwise) analysis and the \emph{transitive} property of linkage.
##' 
##' If the arguments specifying thresholds used to group markers, i.e., minimum
##' LOD Score and maximum recombination fraction, are \code{NULL} (default),
##' the values used are those contained in object \code{input.seq}. If not
##' using \code{NULL}, the new values overridden the ones in object
##' \code{input.seq}.
##' 
##' @aliases group print.group
##' @param input.seq an object of class \code{sequence}.
##' @param LOD a (positive) real number used as minimum LOD score (threshold)
##' to declare linkage.
##' @param max.rf a real number (usually smaller than 0.5) used as maximum
##' recombination fraction to declare linkage.
##' @param x an object of class \code{group}.
##' @param detailed logical. If \code{FALSE}, only a small summary of the
##' linkage groups is printed. If \code{TRUE} (default), the names of markers
##' in each linkage group are also displayed.
##' @param \dots further arguments, passed to other methods. Currently ignored.
##' @return Returns an object of class \code{group}, which is a list containing
##' the following components: \item{data.name}{name of the object of class
##' \code{outcross} that contains the raw data.} \item{twopt}{name of the
##' object of class \code{rf.2ts} used as input, i.e., containing information
##' used to assign markers to linkage groups.} \item{marnames}{marker names,
##' according to the input file.} \item{n.mar}{total number of markers.}
##' \item{LOD}{minimum LOD Score to declare linkage.} \item{max.rf}{maximum
##' recombination fraction to declare linkage.} \item{n.groups}{number of
##' linkage groups found.} \item{groups}{number of the linkage group to which
##' each marker is assigned.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{rf.2pts}} and \code{\link[onemap]{make.seq}}
##' @references Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993)
##' Constructing genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial
##' and reference manual. \emph{A Whitehead Institute for Biomedical Research
##' Technical Report}.
##' @keywords misc
##' @examples
##' 
##'   data(example.out)
##'   twopts <- rf.2pts(example.out)
##'   
##'   all.data <- make.seq(twopts,"all")
##'   link_gr <- group(all.data)
##'   link_gr
##' 
group <-
function(input.seq, LOD=NULL, max.rf=NULL) {
  # checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  
  n.mar <- length(input.seq$seq.num)
  # 'groups' indicates the linkage group to which each marker is associated
  groups <- rep(NA,n.mar)
  
  # determining thresholds
  if (is.null(LOD)) 
    LOD <- get(input.seq$twopt, pos=1)$LOD
  if (is.null(max.rf)) 
    max.rf <- get(input.seq$twopt, pos=1)$max.rf
  
  # 'group.count' is the current linkage group
  group.count <- 1
  significant <- rep(NA,n.mar)
  temp <- matrix(NA,4,2)
  
  for (m in 1:n.mar) {
    if (is.na(groups[m])) {  # marker with index 'm' is not associated to any group
      
      for (k in (1:n.mar)[-m]) {
        # recover values from two-point analyses
        big <- pmax.int(input.seq$seq.num[m],input.seq$seq.num[k])
        small <- pmin.int(input.seq$seq.num[m],input.seq$seq.num[k])
        temp <- get(input.seq$twopt, pos=1)$analysis[acum(big-2)+small,,]

        # check if any assignment meets the criteria
        if (any(temp[,1] <= max.rf & temp[,2] >= LOD)) significant[k] <- "*"
        else significant[k] <- "ns"
      }
      significant[m] <- NA
      
      # check which markers are linked with 'm'
      grouping <- c(m,which(significant=="*"))
      
      if (length(grouping) > 1) {  # 'm' is linked with at least one marker
        begin <- 1 # 'begin' is the index of the last marker added
        flag <- 1 # 'flag' indicates if a new marker has been added to the group
        while (flag) {
          flag <- 0
          next.begin <- length(grouping) # 'next.begin' holds the future value of 'begin'
          
          for (i in tail(grouping,-begin)) {
            # detect all markers linked to those already in group
            significant <- rep(NA,n.mar)
            
			for (k in (1:n.mar)[-i]) {
			  # recover values from two-point analyses
              big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[k])
              small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[k])
              temp <- get(input.seq$twopt, pos=1)$analysis[acum(big-2)+small,,]
              
              # check if any assignment meets the criteria
              if (any(temp[,1] <= max.rf & temp[,2] >= LOD)) significant[k] <- "*"
              else significant[k] <- "ns"
            }
            significant[i] <- NA
            
            group_parc <- which(significant=="*")

            # check if markers in 'group_parc' are already in the current group
            if (any(new.mrk <- is.na(match(group_parc,grouping)))) {
              grouping <- c(grouping,group_parc[new.mrk])
              flag <- 1
            }
          }
          begin <- next.begin
        }
        
        # finishing the current linkage group
        groups[grouping] <- group.count
        group.count <- group.count + 1
      }
    }
  }
  
  ifelse(all(is.na(groups)), n.groups <- 0, n.groups <- max(groups,na.rm=TRUE))
  
  # results
  structure(list(data.name=get(input.seq$twopt, pos=1)$data.name, input.name=deparse(substitute(input.seq)),
                 twopt=input.seq$twopt, marnames=get(input.seq$twopt, pos=1)$marnames,
                 n.mar=n.mar, seq.num=input.seq$seq.num, LOD=LOD, max.rf=max.rf,
                 n.groups=n.groups, groups=groups), class = "group")
}



# print method for object class 'group'
print.group <-
function(x, detailed=TRUE,...) {
  # checking for correct object
  if(!any(class(x)=="group")) stop(deparse(substitute(x))," is not an object of class 'group'")
  
  cat("  This is an object of class 'group'\n")
  cat(paste("  It was generated from the object \"", x$input.name,
            "\"\n\n",sep=""))
  
  # criteria
  cat("  Criteria used to assign markers to groups:\n")
  cat("    LOD =", x$LOD, ", Maximum recombination fraction =",
      x$max.rf, "\n")

  # printing summary
  cat("\n  No. markers:           ", x$n.mar, "\n")
  cat("  No. groups:            ", x$n.groups, "\n")
  cat("  No. linked markers:    ", sum(!is.na(x$groups)), "\n")
  cat("  No. unlinked markers:  ", sum(is.na(x$groups)), "\n")
  
  if (detailed) {
    # printing detailed results (markers in each linkage group)
    cat("\n  Printing groups:")
    for (i in 1:x$n.groups) {
      cat("\n  Group", i, ":", length(which(x$groups==i)) , "markers\n    ")
      cat(x$marnames[x$seq.num[which(x$groups==i)]], "\n")
    }
    if (any(is.na(x$groups))) {
      cat("\n  Unlinked markers:", length(which(is.na(x$groups))) ," markers\n    ")
      cat(x$marnames[x$seq.num[which(is.na(x$groups))]], "\n")
    }
  }
}
##end of file

