#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: marker.type.R                                                 #
# Contains: marker.type                                               #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################



##' Informs the segregation patterns of markers
##'
##' Informs the type of segregation of all markers from an object of class
##' \code{sequence}. For outcross populations it uses the notation by \cite{Wu
##' et al., 2002}. For backcrosses, F2s and RILs, it uses the
##' traditional notation from MAPMAKER i.e. AA, AB, BB, not AA and not BB.
##'
##' The segregation types are (\cite{Wu et al., 2002}): \tabular{lcc}{ Type
##' \tab Cross \tab Segregation \cr A.1 \tab ab x cd \tab 1:1:1:1 \cr A.2 \tab
##' ab x ac \tab 1:1:1:1 \cr A.3 \tab ab x co \tab 1:1:1:1 \cr A.4 \tab ao x bo
##' \tab 1:1:1:1 \cr B1.5 \tab ab x ao \tab 1:2:1 \cr B2.6 \tab ao x ab \tab
##' 1:2:1 \cr B3.7 \tab ab x ab \tab 1:2:1 \cr C8 \tab ao x ao \tab 3:1 \cr
##' D1.9 \tab ab x cc \tab 1:1 \cr D1.10 \tab ab x aa \tab 1:1 \cr D1.11 \tab
##' ab x oo \tab 1:1 \cr D1.12 \tab bo x aa \tab 1:1 \cr D1.13 \tab ao x oo
##' \tab 1:1 \cr D2.14 \tab cc x ab \tab 1:1 \cr D2.15 \tab aa x ab \tab 1:1
##' \cr D2.16 \tab oo x ab \tab 1:1 \cr D2.17 \tab aa x bo \tab 1:1 \cr D2.18
##' \tab oo x ao \tab 1:1 }
##'
##' @param input.seq an object of class \code{sequence}.
##' @return Nothing is returned. Segregation types of all markers in the
##' sequence are displayed on the screen.
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{make.seq}}
##' @references Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002)
##' Simultaneous maximum likelihood estimation of linkage and linkage phases in
##' outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords manip utilities
##' @examples
##'
##'   data(example.out)
##'   twopts <- rf.2pts(example.out)
##'   markers.ex <- make.seq(twopts,c(3,6,8,12,16,25))
##'   marker.type(markers.ex) # segregation type for some markers
##'
##'   data(fake.f2.onemap)
##'   twopts <- rf.2pts(fake.f2.onemap)
##'   all.mrk<-make.seq(twopts, "all")
##'   lgs<-group(all.mrk)
##'   lg1<-make.seq(lgs,1)
##'   marker.type(lg1) # segregation type for linkage group 1
##'
##'
marker.type <-
function(input.seq) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")

  ## printing marker type
  if(class(get(input.seq$data.name, pos=1))=="outcross"){
    for(i in 1:length(input.seq$seq.num))
      cat("  Marker", input.seq$seq.num[i], "(", colnames(get(input.seq$twopt)$analysis[[1]])[input.seq$seq.num[i]], ") has type", get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num[i]], "\n")
  }
  else{
    for(i in 1:length(input.seq$seq.num)){
      mrk.type<-rep("NA",length(input.seq$seq.num))
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="C.A"]<-"Not  AA : AA (3:1) "
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="D.B"]<-"Not  BB : BB (3:1) "
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="A.H.B"]<-"AA : AB : BB (1:2:1) "
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="M.X"]<-"Mixed: Dominant & Co-dominant"
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="A.H"]<-"AA : AB (1:1)"
      mrk.type[get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="A.B"]<-"AA : BB (1:1)"

      cat("  Marker", input.seq$seq.num[i], "(", colnames(get(input.seq$twopt)$analysis)[input.seq$seq.num[i]], ") -->", mrk.type[i], "\n")
    }
  }
}

## end of file
