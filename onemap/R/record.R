#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: record.R                                                      #
# Contains: record                                                    #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2007-9, Marcelo Mollinari                             #
#                                                                     #
# First version: 11/29/2009                                           #
# Last update: 01/21/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

record<-function(input.seq, times=10, LOD=0, max.rf=0.5, tol=10E-5){
  ## checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
  markers <- length(input.seq$seq.num)

  ## create recombination fraction matrix 
  r <- matrix(NA,markers,markers)
  for(i in 1:(markers-1)) {
    for(j in (i+1):markers) {
      big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[j])
      small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[j])
      temp <- get(input.seq$twopt)$analysis[acum(big-2)+small,,]
      ## check if any assignment meets the criteria
      relevant <- which(temp[,2] > (max(temp[,2])-0.005)) # maximum LOD scores
      phases <- relevant[which((temp[relevant,1] <= max.rf & temp[relevant,2] >= LOD))]
      if(length(phases) == 0) r[i,j] <- r[j,i] <- 0.5
      else r[i,j] <- r[j,i] <- temp[phases[1],1]
    }
  }
  diag(r)<-0
  
  ##RECORD algorithm
  n.mark<-ncol(r)
  X<-r*get(input.seq$data.name, pos=1)$n.ind ## Obtaining X multiplying the MLE of the recombination
                                      ## fraction by the number of individuals
  
  COUNT<-function(X, sequence){ ## See eq. 1 on the paper (Van Os et al., 2005)
    return(sum(diag(X[sequence[-length(sequence)],sequence[-1]]),na.rm=TRUE))
  }

  ## For two markers
  if(n.mark==2)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))
  
  ## For three markers (calculation of 3 possible orders: 3!/2)
  else if(n.mark==3) {
    all.perm<-perm.pars(1:3)
    m.old<-Inf
    for(k in 1:nrow(all.perm)){
      m.new<-COUNT(X, all.perm[k,])
      if(m.new < m.old)
        result.new <- all.perm[k,]; m.old<-m.new
    }
  }
  
  ## For more than three markers (RECORD algorithm itself)
  else{
    marks<-c(1:n.mark)
    result.new<-sample(marks)## randomize markers
    for(k in 1:times){ ## loop for replicates the RECORD procedure
      result<-sample(marks, 2)
      for(i in 2:(n.mark-2)){
        next.mark<-sample(c(1:n.mark)[-result],1)
        partial<-rep(NA,2*length(result)+1)
        partial[seq(2,(length(partial)-1), by = 2)]<-result
        temp<-10e1000
        for(j in seq(1,(length(partial)), by = 2)){
          partial.temp<-partial
          partial.temp[j]<-next.mark
          temp.new<-COUNT(X, partial.temp[is.na(partial.temp)==FALSE])  
          if(temp.new<temp) {
            result<-partial.temp[is.na(partial.temp)==FALSE]
            temp<-temp.new
          }
        }
      }
      next.mark<-c(1:n.mark)[-result]
      partial<-rep(NA,2*length(result)+1)
      partial[seq(2,(length(partial)-1), by = 2)]<-result
      temp<-10e1000
      for(j in seq(1,(length(partial)), by = 2)){
        partial.temp<-partial
        partial.temp[j]<-next.mark
        temp.new<-COUNT(X, partial.temp[ is.na(partial.temp)==FALSE])  
        if(temp.new<temp) {
          result<-partial.temp[ is.na(partial.temp)==FALSE]
          temp<-temp.new
        }
      }
      for(i in 1:(n.mark-1)){
        for(j in 1:(n.mark-i)){
          y<-result[j:(j+i)]
          if(j==1){
            if(i!=n.mark-1){
              perm.order<-c(rev(result[1:(1+i)]),result[(i+2):length(result)])
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp) {
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
          else{
            if(j!=(n.mark-i)){
              perm.order<-c(result[1:(j-1)],rev(result[j:(j+i)]),result[(j+i+1):length(result)])
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }
            else{
              perm.order<-c(result[1:(j-1)],rev(result[j:length(result)]))
              COUNT.temp<-COUNT(X,perm.order)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }      
          }
        }
      }
      if(COUNT(X,result.new) > COUNT(X,result)){
        result.new<-result
        ##print(COUNT(X,result.new))
      }
    }
  }
  ## end of RECORD algorithm
 cat("\norder obtained using RECORD algorithm:\n\n", input.seq$seq.num[avoid.reverse(result.new)], "\n\ncalculating multipoint map using tol", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(result.new)],twopt=input.seq$twopt), tol=tol)
}

##end of file
  
















  
