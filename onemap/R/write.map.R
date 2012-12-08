#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: write.map.R                                                   #
# Contains: write.map                                                 #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2010, Marcelo Mollinari                               #
#                                                                     #
# First version: 10/10/2010                                           #
# Last update: 11/30/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function creates Rqtl input files.
write.map<-function(map.list,file.out){
   # checking for correct object
  if(!any(class(map.list)=="list" | class(map.list)=="sequence")) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequnece'")
  if(!any(class(file.out)=="character")) stop(deparse(substitute(file.out))," is an invalid output file name")
  
  # if map.list is just a single chormosome, convert it  into a list
  if(class(map.list)=="sequence") map.list<-list(map.list)
  
  write(x="",file=file.out)
  for(i in 1:length(map.list)){
    if(!any(class(map.list[[i]])=="sequence")) stop("Object ", i , " in map.list is not an object of class 'sequnece'")
    if(is.null(map.list[[i]]$seq.like))  stop("Parameters are not estimated for object ", i, " in map.list")
    map<-cumsum(c(0,get(get(".map.fun", envir=.onemapEnv))(map.list[[i]]$seq.rf)))
    marnames<-colnames(get(map.list[[i]]$data.name, pos=1)$geno)[map.list[[i]]$seq.num]
    write(t(cbind(i,marnames,map)), file=file.out, ncolumns =3, append=TRUE)
  }
}

