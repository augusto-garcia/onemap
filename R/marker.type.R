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

marker.type <-
function(input.seq) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")

  ## printing marker type
  if(class(get(input.seq$data.name, pos=1))=="outcross"){
    for(i in 1:length(input.seq$seq.num))
      cat("  Marker", input.seq$seq.num[i], "(", get(input.seq$twopt)$marnames[input.seq$seq.num[i]], ") has type", get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num[i]], "\n")
  }
  else if(class(get(input.seq$data.name, pos=1))=="f2.onemap" || class(get(input.seq$data.name, pos=1))=="bc.onemap"){
    for(i in 1:length(input.seq$seq.num)){
      mrk.type<-rep("NA",length(input.seq$seq.num))
      mrk.type[which(get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="C.8" & get(input.seq$data.name, pos=1)$phase[input.seq$seq.num] == -1)]<-"Not  AA : AA (3:1) "
      mrk.type[which(get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="C.8" & get(input.seq$data.name, pos=1)$phase[input.seq$seq.num] == 1)]<-"Not  BB : BB (3:1) "
      mrk.type[which(get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="B3.7" & get(input.seq$data.name, pos=1)$phase[input.seq$seq.num] == 1)]<-"AA : AB : BB (1:2:1) "      
      mrk.type[which(get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="D1.10")]<-"AA : AB (1:1)"
      
      cat("  Marker", input.seq$seq.num[i], "(", get(input.seq$twopt)$marnames[input.seq$seq.num[i]], ") -->", mrk.type[i], "\n")
    }  
  }
  else if(class(get(input.seq$data.name, pos=1))=="riself.onemap" || class(get(input.seq$data.name, pos=1))=="risib.onemap"){
    for(i in 1:length(input.seq$seq.num)){
      mrk.type<-rep("NA",length(input.seq$seq.num))
      mrk.type[which(get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]=="D1.10")]<-"AA : BB (1:1)"

      cat("  Marker", input.seq$seq.num[i], "(", get(input.seq$twopt)$marnames[input.seq$seq.num[i]], ") -->", mrk.type[i], "\n")
    }  
  }
  else stop("There are invalid class")
}

## end of file
