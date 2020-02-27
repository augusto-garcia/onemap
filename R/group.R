#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: group.R                                                       ##
## Contains: check_linkage, group, print.group                         ##
##                                                                     ##
## Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  ##
## copyright (c) 2007-9, Gabriel R A Margarido and Marcelo Mollinari   ##
##                                                                     ##
## First version: 11/07/2007                                           ##
## Last update: 21/06/2016                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

## Function to assign markers to linkage groups

##' Assign markers to linkage groups
##'
##' Identifies linkage groups of markers, using results from two-point
##' (pairwise) analysis and the \emph{transitive} property of linkage.
##'
##' If the arguments specifying thresholds used to group markers, i.e., minimum
##' LOD Score and maximum recombination fraction, are \code{NULL} (default),
##' the values used are those contained in object \code{input.seq}. If not
##' using \code{NULL}, the new values override the ones in object
##' \code{input.seq}.
##'
##' @aliases group
##' @param input.seq an object of class \code{sequence}.
##' @param LOD a (positive) real number used as minimum LOD score
##'     (threshold) to declare linkage.
##' @param max.rf a real number (usually smaller than 0.5) used as
##'     maximum recombination fraction to declare linkage.
##' @param verbose logical. If \code{TRUE}, current progress is shown;
##'     if \code{FALSE}, no output is produced.
##' @return Returns an object of class \code{group}, which is a list
##'     containing the following components: \item{data.name}{name of
##'     the object of class \code{onemap} that contains the raw
##'     data.} \item{twopt}{name of the object of class \code{rf.2ts}
##'     used as input, i.e., containing information used to assign
##'     markers to linkage groups.} \item{marnames}{marker names,
##'     according to the input file.} \item{n.mar}{total number of
##'     markers.}  \item{LOD}{minimum LOD Score to declare linkage.}
##'     \item{max.rf}{maximum recombination fraction to declare
##'     linkage.} \item{n.groups}{number of linkage groups found.}
##'     \item{groups}{number of the linkage group to which each marker
##'     is assigned.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com} and
##'     Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{rf_2pts}} and
##'     \code{\link[onemap]{make_seq}}
##' @references Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993)
##'     Constructing genetic linkage maps with MAPMAKER/EXP Version
##'     3.0: a tutorial and reference manual. \emph{A Whitehead
##'     Institute for Biomedical Research Technical Report}.
##' @keywords misc
##' @examples
##'
##'   data(onemap_example_out)
##'   twopts <- rf_2pts(onemap_example_out)
##'
##'   all.data <- make_seq(twopts,"all")
##'   link_gr <- group(all.data)
##'   link_gr
##'   print(link_gr, details=FALSE) #omit the names of the markers
##'@export
group <- function(input.seq, LOD=NULL, max.rf=NULL, verbose=TRUE)
{
    ## checking for correct object
    if(!is(input.seq,"sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
    ## determining thresholds
    if (is.null(LOD))
        LOD <- get(input.seq$twopt, pos=1)$LOD
    if (is.null(max.rf))
        max.rf <- get(input.seq$twopt, pos=1)$max.rf
    cl<-class(get(input.seq$data.name))[2]
    geno<-get(input.seq$data.name)$geno[,input.seq$seq.num]
    st<-get(input.seq$data.name)$segr.type.num[input.seq$seq.num]
    groups<-rep(0, length(input.seq$seq.num))
    tp<-list(unlk=1:length(input.seq$seq.num))
    i<-1
    if(verbose) cat("   Selecting markers: \n")
    while(length(tp$unlk) > 0)
    {
        g<-tp$unlk[1]
        s<-tp$unlk
        j<-1
        tp<-check_linkage(i=g[j], s=s, cl=cl, geno=geno, st=st, max.rf=max.rf, LOD=LOD)
        gt<-tp$lk
        if(length(gt) > 0)
        {
            if(verbose) cat("\t  group   ", i,"\n\t   ")
            g<-c(g,gt)
            while(!is.na(g[j+1])){
                if(verbose)
                    {
                        cat(".")
                        if(j %% 60 == 0) cat("\n\t   ")
                }
                
                tp<-check_linkage(i=g[j+1], s=tp$unlk, cl=cl, geno=geno, st=st, max.rf=max.rf, LOD=LOD)
                gt<-tp$lk
                g<-c(g,gt)
                j<-j+1
            }
            if(verbose) cat("\n")
            groups[g]<-i
            i<-i+1
        }
    }
    if(all(groups==0)) cat("\t No group found.\n")
    ## results
    structure(list(data.name=input.seq$data.name, input.name=deparse(substitute(input.seq)),
                   twopt=input.seq$twopt, marnames=colnames(geno),
                   n.mar=length(input.seq$seq.num), seq.num=input.seq$seq.num, LOD=LOD, max.rf=max.rf,
                   n.groups=i-1, groups=groups), class = "group")
}

##' Show the results of grouping procedure
##'
##' It shows the linkage groups as well as the unlinked markers.
##'
##' @aliases print.group
##' @param x an object of class onemap_segreg_test
##'
##' @param detailed logical. If \code{TRUE} the markers in each
##'     linkage group are printed.
##'
##' @param ... currently ignored
##' @return \code{NULL}
##' @keywords internal
##' @method print group
##' @export
print.group <-
    function(x, detailed=TRUE,...) {
        ## checking for correct object
        if(!is(x,"group")) stop(deparse(substitute(x))," is not an object of class 'group'")

        cat("  This is an object of class 'group'\n")
        cat(paste("  It was generated from the object \"", x$input.name,
                  "\"\n\n",sep=""))

        ## criteria
        cat("  Criteria used to assign markers to groups:\n")
        cat("    LOD =", x$LOD, ", Maximum recombination fraction =",
            x$max.rf, "\n")

        ## printing summary
        cat("\n  No. markers:           ", x$n.mar, "\n")
        cat("  No. groups:            ", x$n.groups, "\n")
        cat("  No. linked markers:    ", sum(x$groups > 0), "\n")
        cat("  No. unlinked markers:  ", sum(x$groups == 0), "\n")

        if (detailed) {
            ## printing detailed results (markers in each linkage group)
            cat("\n  Printing groups:")
            for (i in 1:x$n.groups) {
                cat("\n  Group", i, ":", length(which(x$groups==i)) , "markers\n    ")
                cat(x$marnames[which(x$groups==i)], "\n")
            }
            if (any(x$groups==0)) {
                cat("\n  Unlinked markers:", length(which(x$groups==0)) ," markers\n    ")
                cat(x$marnames[which(x$groups==0)], "\n") 
            }
        }
    }


##Checks if a marker i is linked with markers in a vector s
check_linkage<-function(i, s, cl, geno, st=NULL, max.rf, LOD)
{
    s<-s[is.na(match(s,i))]
    if(cl=="outcross")
    {
        r<-est_rf_out(geno = geno[,c(i,s)], mrk = 1, seg_type = st[c(i,s)], nind = nrow(geno))
        sig<-apply(r[[1]], 2, function(x,y) min(x) <= y, y=max.rf) &
            apply(r[[2]], 2, function(x,y) max(x) >= y, y=LOD)
    }
    else if(cl=="f2")
    {
        r<-est_rf_f2(geno = geno[,c(i,s)], mrk = 1, seg_type = st[c(i,s)], nind = nrow(geno))
        sig<-r[1,] <= max.rf & r[2,] >=LOD
    }
    else if(cl=="backcross")
    {
        r<-est_rf_bc(geno = geno[,c(i,s)], mrk = 1, type = 0, nind = nrow(geno))
        sig<-r[1,] <= max.rf & r[2,] >=LOD
    }
    else if(cl=="riself")
    {
        r<-est_rf_bc(geno = geno[,c(i,s)], mrk = 1, type = 1, nind = nrow(geno))
        sig<-r[1,] <= max.rf & r[2,] >=LOD
    }

    else if(cl=="risib")
    {
        r<-est_rf_bc(geno = geno[,c(i,s)], mrk = 1, type = 1, nind = nrow(geno))
        sig<-r[1,] <= max.rf & r[2,] >=LOD
    }
    return(list(lk=s[sig[-1]], unlk=s[!(sig[-1])]))
}
###end of file
