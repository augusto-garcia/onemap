#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: comb.ger.R                                                    #
# Contains: comb, comb.ger, diplo, rem.amb.ph                         #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 04/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## This function combines two linkage phase vectors 
comb <-
function(x,y) {
  count <- 0
  M <- matrix(NA, nrow(x)*nrow(y), 1+ncol(y))
  for(i in 1:nrow(x)) {
    for(j in 1:nrow(y)) {
      count <- count+1
      M[count,] <- c(x[i,],y[j,])
    }
  }
  return(M)
}

##This function creates diplotypes from segregation types and linkage phases
diplo <-
function(w, seq.num, seq.phases) {
    # convert numerical linkage phases to strings
    link.phases <- matrix(NA,length(seq.num),2)
    link.phases[1,] <- rep(1,2)
    for (i in 1:length(seq.phases)) {
      switch(EXPR=seq.phases[i],
             link.phases[i+1,] <- link.phases[i,]*c(1,1),
             link.phases[i+1,] <- link.phases[i,]*c(1,-1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
             )
    }
    ## create diplotypes from segregation types and linkage phases
    link.phases <- apply(link.phases,1,function(x) paste(as.character(x),collapse="."))
    parents <- matrix("",length(seq.num),4)
    for (i in 1:length(seq.num))
      parents[i,] <- return.geno(get(w$data.name, pos=1)$segr.type[seq.num[i]],link.phases[i])
      return(parents)
}


##This function removes ambiguous phases based on identical diplotypes
rem.amb.ph <-
function(M,w,seq.num) {
  M.new<-matrix(NA,nrow(M),length(seq.num)*4)
  for(j in 1:nrow(M)){
    M.new[j,]<-as.vector(diplo(w=w, seq.num=seq.num, M[j,]))
  }
  v<-which(duplicated(M.new)==FALSE)
  if(length(v)>1){
    k<-numeric()    
    for(i in 1:(length(v)-1)){
      for(j in (i+1):length(v)){
        if(all(diplo(w=w, seq.num=seq.num, M[v[i],])==diplo(w=w, seq.num=seq.num, M[v[j],])[,c(2,1,3,4)])){
         k<-c(k,j)         
        }
        if(all(diplo(w=w, seq.num=seq.num, M[v[i],])==diplo(w=w, seq.num=seq.num, M[v[j],])[,c(1,2,4,3)])){
          k<-c(k,j)         
        }
        if(all(diplo(w=w, seq.num=seq.num, M[v[i],])==diplo(w=w, seq.num=seq.num, M[v[j],])[,c(2,1,4,3)])){
          k<-c(k,j)
        }
      }
    }
    if(length(k)==0) return(v)
    else return(v[-k])
  }
  else return(v)
}


# This function makes all possible combinations of a list
# containing the linkage phase vectors for each interval
comb.ger <- function(f){
  M <- as.matrix(f[[length(f)]])
  if (length(f)==1) return(M)
  else{
    for(i in (length(f)-1):1){
      M <- comb(as.matrix(f[[i]]),M)
    }
    return(M)
  }
}




# end of file


