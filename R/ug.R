#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: ug.R                                                          #
# Contains: ug                                                        #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2007-9, Marcelo Mollinari                             #
#                                                                     #
# First version: 11/29/2009                                           #
# Last update: 01/21/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

ug<-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5){
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
  n.mark<-ncol(r)
  
  ## For two markers
  if(n.mark==2)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))

  ## UG algorithm
  ## The equation numbers below refer to the ones in the article (Tan and Fu, 2006)
  
  ##eq 1 and eq 2
  calc.T<-function(d,i,j) 2*d[i,j]-(sum(d[i,-i], na.rm = TRUE)+sum(d[j,-j], na.rm = TRUE))
  
  ##Function to pick the position that has the minimum value (considering ties)
  pick.pos.min<-function(x){
    if (length(which(min(x,na.rm=TRUE)==x))>1) return(sample((which(min(x,na.rm=TRUE)==x)),1))
    else return(which(min(x,na.rm=TRUE)==x))
  } 
  lambda<-matrix(1,n.mark,n.mark) # see the last paragraph on page 2385
  
  ##Function to calculate distance between loci i and j
  ##eq.13
  r.to.d<-function(r,i,j,lambda){
    N<-0;q<-0
    for(k in 1:n.mark){
      if(r[i,j]>r[i,k] && r[i,j]>r[j,k]){ 
        q<-q+(r[i,k]*r[j,k])
        N<-N+1
      }
    }
    if(N==0) N<-1
    return(r[i,j]+(2*lambda[i,j]/N*q))
  }
  
  ##Calculating d
  d<-matrix(NA,(2*n.mark-1),(2*n.mark-1))
  for(i in 1:n.mark){
    for(j in i:n.mark){
      d[i,j]<-d[j,i]<-r.to.d(r,i,j,lambda)
    }
  }
  
  ##Calculating T 
  T<-matrix(NA,n.mark,n.mark)
  for(i in 1:n.mark){
    for(j in i:n.mark){
      T[i,j]<-calc.T(d,j,i)
    }
    diag(T)<-NA
  }
  
  ##verification of eq 3
  ##linkage<-c(T[1,2],T[19,20],T[20,21])
  ##min(linkage)== min(T, na.rm=T)
  
  
  ##UG algorithm itself
  
  ##step 1
  ##Finding the smallest T-value
  partial<-which(min(T,na.rm=TRUE)==T, arr.ind = TRUE)
  if(length(partial)>2) partial<-partial[sample(c(1:nrow(partial)),1),]
  new.locus<-n.mark+1
  
  din1<-function(x,y,d,n.mark){
    d.new.v<-rep(NA, n.mark)
    for(i in 1:n.mark){
      if((d[i,x]+d[i,y])>d[x,y]){
        d.new.v[i]<-.5*(d[i,x]+d[i,y]-d[x,y]) #eq7
      }
      else d.new.v[i]<-0
    }
    return(d.new.v)
  }
  
  d.old1<-d[partial[1],]
  d.old2<-d[partial[2],]
  d[new.locus,c(1:n.mark)]<-d[c(1:n.mark),new.locus]<-d[new.locus,c(1:n.mark)]<-din1(partial[1],partial[2],d,n.mark)
  d[,partial[2]]<-d[partial[2],]<-d[,partial[1]]<- d[partial[1],]<-NA
  H.temp<-c()
  for(i in 1:n.mark){
    H.temp[i]<-(n.mark-2)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE) #eq8
  }
  
  H<-H.temp
  next.mark<-pick.pos.min(H)
  side<-c(d.old1[next.mark],d.old2[next.mark]) 
  if(pick.pos.min(side)==1){
    a<-partial[1]
    partial[1]<-partial[2]
    partial[2]<-a
  }
  partial<-c(partial, next.mark)
  
  ## If there are three markers, do not go to the second step
  if(n.mark==3)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(partial)],twopt=input.seq$twopt), tol=10E-5))
  
  for (k in 2:(n.mark-2)){
    ##step 2
    new.locus<-n.mark+k
    for(i in 1:(new.locus-1)){
      d[i,new.locus]<-d[new.locus,i]<-min(d[(new.locus-1),i], d[partial[k+1],i])#eq9
    }
    d[partial[k+1],]<-NA
    d[,partial[k+1]]<-NA
    d[new.locus-1,]<-NA
    d[,new.locus-1]<-NA
    
    ##step 3
    H.temp<-c()
    for(i in 1:n.mark){
      H.temp[i]<-(n.mark-k-1)*d[new.locus,i]-sum(d[i,-i], na.rm = TRUE)#eq10
    }
    
    H<-rbind(H,H.temp)
    next.mark<-pick.pos.min(H[k,])
    partial<-c(partial, next.mark)
  }
  complete<-partial
  ## end of UG algorithm
  cat("\norder obtained using UG algorithm:\n\n", input.seq$seq.num[avoid.reverse(complete)], "\n\ncalculating multipoint map using tol ", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(complete)],twopt=input.seq$twopt), tol=tol)
  
}
## end of file
  
















  
