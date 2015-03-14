#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: codif.data.R                                                  #
# Contains: codif.data                                                #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function gets the input in strings and converts into numbers
codif.data <-
function(geno.in, segr.type.in) {
  geno.out <- matrix(NA,nrow(geno.in),ncol(geno.in))
  segr.type.out <- rep(NA,length(segr.type.in))

  # missing data are represented by '0'
  geno.out[is.na(geno.in)] <- 0

  for(i in 1:length(segr.type.in)) {
    # based on the marker type, convert strings to numbers
    switch(EXPR=segr.type.in[i],
           A.1={
             geno.out[which(geno.in[,i]=="ac"),i] <- 1
             geno.out[which(geno.in[,i]=="ad"),i] <- 2
             geno.out[which(geno.in[,i]=="bc"),i] <- 3
             geno.out[which(geno.in[,i]=="bd"),i] <- 4
             segr.type.out[i] <- 1
           },
           A.2={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ac"),i] <- 2
             geno.out[which(geno.in[,i]=="ba"),i] <- 3
             geno.out[which(geno.in[,i]=="bc"),i] <- 4
             segr.type.out[i] <- 1
           },         
           A.3={
             geno.out[which(geno.in[,i]=="ac"),i] <- 1
             geno.out[which(geno.in[,i]=="a"),i] <- 2
             geno.out[which(geno.in[,i]=="bc"),i] <- 3
             geno.out[which(geno.in[,i]=="b"),i] <- 4
             segr.type.out[i] <- 1
           },
           A.4={
             geno.out[which(geno.in[,i]=="ab"),i] <- 1
             geno.out[which(geno.in[,i]=="a"),i] <- 2
             geno.out[which(geno.in[,i]=="b"),i] <- 3
             geno.out[which(geno.in[,i]=="o"),i] <- 4
             segr.type.out[i] <- 1
           },
           B1.5={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ab"),i] <- 2
             geno.out[which(geno.in[,i]=="b"),i] <- 3
             segr.type.out[i] <- 2
           },
           B2.6={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ab"),i] <- 2
             geno.out[which(geno.in[,i]=="b"),i] <- 3
             segr.type.out[i] <- 3
           },
           B3.7={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ab"),i] <- 2
             geno.out[which(geno.in[,i]=="b"),i] <- 3
             segr.type.out[i] <- 4
           },
           C.8={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="o"),i] <- 2
             segr.type.out[i] <- 5
           },
           D1.9={
             geno.out[which(geno.in[,i]=="ac"),i] <- 1
             geno.out[which(geno.in[,i]=="bc"),i] <- 2
             segr.type.out[i] <- 6
           },
           D1.10={      
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ab"),i] <- 2
             segr.type.out[i] <- 6
           },
           D1.11={      
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="b"),i] <- 2
             segr.type.out[i] <- 6
           },
           D1.12={      
             geno.out[which(geno.in[,i]=="ab"),i] <- 1
             geno.out[which(geno.in[,i]=="a"),i] <- 2
             segr.type.out[i] <- 6
           },
           D1.13={      
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="o"),i] <- 2
             segr.type.out[i] <- 6
           },
           D2.14={      
             geno.out[which(geno.in[,i]=="ac"),i] <- 1
             geno.out[which(geno.in[,i]=="bc"),i] <- 2
             segr.type.out[i] <- 7
           },
           D2.15={      
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="ab"),i] <- 2
             segr.type.out[i] <- 7
           },
           D2.16={      
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="b"),i] <- 2
             segr.type.out[i] <- 7
           },
           D2.17={      
             geno.out[which(geno.in[,i]=="ab"),i] <- 1
             geno.out[which(geno.in[,i]=="a"),i] <- 2
             segr.type.out[i] <- 7
           },
           D2.18={
             geno.out[which(geno.in[,i]=="a"),i] <- 1
             geno.out[which(geno.in[,i]=="o"),i] <- 2
             segr.type.out[i] <- 7
           }
           )
    if(any(is.na(geno.out[,i])))
      stop(paste("Invalid marker codification. Please check data for marker", colnames(geno.in)[i]), ".", sep="")
  }
  dimnames(geno.out) <- dimnames(geno.in)
  return(list(geno.out,segr.type.out))
}

# end of file
