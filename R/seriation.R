#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: seriation.R                                                   #
# Contains: seriation, ser.ord, Cindex                                #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/29/2009                                           #
# Last update: 11/29/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

seriation<-function(input.seq, LOD=0, max.rf=0.5, tol=10E-5){
    ## checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
  markers <- length(input.seq$seq.num)
  r <- matrix(NA,markers,markers)
  for(i in 1:(markers-1)) {
    for(j in (i+1):markers) {
      big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[j])
      small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[j])
      temp <- get(input.seq$twopt)$analysis[acum(big-2)+small,,]

      ## check if any assignment meets the criteria
	  relevant <- which(temp[,2] > (max(temp[,2])-0.005)) ## maximum LOD scores
      phases <- relevant[which((temp[relevant,1] <= max.rf & temp[relevant,2] >= LOD))]
	  if(length(phases) == 0) r[i,j] <- r[j,i] <- 0.5
	  else r[i,j] <- r[j,i] <- temp[phases[1],1]
    }
  }
  diag(r)<-0
  
  ## SERIATION algorithm
  n.mark<-ncol(r)
  orders <- array(0,dim=c(n.mark,n.mark))
  CI <- numeric(n.mark)
  for (i in 1:n.mark) {
    orders[i,] <- ser.ord(r,i)
    CI[i] <- Cindex(orders[i,],r)
  }
  best <- which(CI==CI[which.min(CI)])
  if (length(best) == 0) stop ("Cannot find any order")
  else {
    if (length(best) != 1) best <- best[sample(length(best),1)]
    complete<-orders[best,]
  }

  ## end of SERIATION algorithm
  cat("\norder obtained using SERIATION algorithm:\n\n", input.seq$seq.num[complete], "\n\ncalculating multipoint map using tol = ", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[complete],twopt=input.seq$twopt), tol=tol)
}

##Provides an order given the recombination
##fraction matrix and the starting marker. 
ser.ord <- function(r,i) {
  markers <- ncol(r)
  x <- 1:markers
  unres1 <- 0
  unres2 <- 0
  esq <- numeric(0)
  dir <- numeric(0)
  order <- i
  x[i] <- NaN
  j <- which(r[i,][x]==r[i,which.min(r[i,][x])])
  if (length(j) > 1) j <- j[sample(length(j),1)]
  order <- c(order,j)
  x[j] <- NaN
  while (length(order) != markers || unres1 || unres2[1]) {
    if (unres1) {
      if ((length(order) + length(esq) + length(dir) + 1)==markers) {
        duv <- sample(2,1)
        if (duv==1) order <- c(esq,unres1,order,dir) else order <- c(esq,order,unres1,dir)
        unres1 <- 0; dir <- numeric(0); esq <- numeric(0)
      }
      else {
        pri <- if (length(esq)==0) order[1] else esq[1]
        ult <- if (length(dir)==0) order[length(order)] else dir[length(dir)]
        e <- which(r[i,][x]==r[i,which.min(r[i,][x])])
        if (length(e) > 1) e <- e[sample(length(e),1)]
        rand <- 0
        if (r[pri,e] == r[ult,e]) rand <- sample(2,1)
        if (r[pri,e] < r[ult,e] || rand==1) {
          if (r[e,unres1] < r[ult,unres1]) { order <- c(e,esq,unres1,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else if (r[e,unres1] > r[ult,unres1]) { order <- c(e,esq,order,unres1,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else { esq <- c(e,esq); x[e] <- NaN }
        }
        else if (r[pri,e] > r[ult,e] || rand==2) {
          if (r[e,unres1] < r[pri,unres1]) { order <- c(esq,order,unres1,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else if (r[e,unres1] > r[pri,unres1]) { order <- c(esq,unres1,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres1 <- 0 }
          else { dir <- c(dir,e); x[e] <- NaN }
        }
      }
    }
    else if (unres2[1]) {
      if ((length(order) + length(esq) + length(dir) + 2)==markers) {
        if (unres2[1]==1) {
          duv <- sample(2,1)
          if (duv==1) order <- c(esq,unres2[2],unres2[3],order,dir) else order <- c(esq,unres2[3],unres2[2],order,dir)
          unres2 <- 0; dir <- numeric(0); esq <- numeric(0)
        }
        else if (unres2[1]==2) {
          duv <- sample(2,1)
          if (duv==1) order <- c(esq,order,unres2[2],unres2[3],dir) else order <- c(esq,order,unres2[3],unres2[2],dir)
          unres2 <- 0; dir <- numeric(0); esq <- numeric(0)
        }
      }
      else {
        pri <- if (length(esq)==0) order[1] else esq[1]
        ult <- if (length(dir)==0) order[length(order)] else dir[length(dir)]
        e <- which(r[i,][x]==r[i,which.min(r[i,][x])])
        if (length(e) > 1) e <- e[sample(length(e),1)]
        rand <- 0
        if (r[pri,e] == r[ult,e]) rand <- sample(2,1)
        m1 <- unres2[2]; m2 <- unres2[3]
        if (r[pri,e] < r[ult,e] || rand==1) {
          if (unres2[1]==1) {
            if (r[e,m1] < r[e,m2]) { order <- c(e,esq,m1,m2,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(e,esq,m2,m1,order,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { esq <- c(e,esq); x[e] <- NaN }
          }
          else if (unres2[1]==2) {
            if (r[e,m1] < r[e,m2]) { order <- c(e,esq,order,m1,m2,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(e,esq,order,m2,m1,dir); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { esq <- c(e,esq); x[e] <- NaN }	
          }
        }
        else if (r[pri,e] > r[ult,e] || rand==2) {
          if (unres2[1]==1) {
            if (r[e,m1] < r[e,m2]) { order <- c(esq,m2,m1,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(esq,m1,m2,order,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { dir <- c(dir,e); x[e] <- NaN }
          }
          else if (unres2[1]==2) {
            if (r[e,m1] < r[e,m2]) { order <- c(esq,order,m2,m1,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else if (r[e,m1] > r[e,m2]) { order <- c(esq,order,m1,m2,dir,e); x[e] <- NaN; esq <- numeric(0); dir <- numeric(0); unres2 <- 0 }
            else { dir <- c(dir,e); x[e] <- NaN }	
          }
        }
      }
    }
    else {
      if (length(which(r[i,][x]==r[i,which.min(r[i,][x])]))==1) {
        j <- which.min(r[i,][x])
        if (r[order[1],j] < r[order[length(order)],j]) { order <- c(j,order); x[j] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j]) { order <- c(order,j); x[j] <- NaN }
        else {
          light <- 1
          k <- 2
          l <- length(order)-1
          while (k < l && light) {
            if (r[order[k],j] < r[order[l],j]) { order <- c(j,order); x[j] <- NaN; light <- 0 }
            else if (r[order[k],j] > r[order[l],j]) { order <- c(order,j); x[j] <- NaN; light <- 0 }
            else { k <- k+1; l <- l-1 }
          }
          if (light) { unres1 <- j; x[j] <- NaN }
        }
      }
      else if (length(which(r[i,][x]==r[i,which.min(r[i,][x])]))==2) {
        j <- which(r[i,][x]==r[i,which.min(r[i,][x])])[1]
        k <- which(r[i,][x]==r[i,which.min(r[i,][x])])[2]
        if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) { order <- c(j,order,k); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) { order <- c(k,order,j); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) {
          if (r[order[1],j] < r[order[1],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[1],j] > r[order[1],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN }
          else {
            l <- 2
            light <- 1
            while (l <= length(order) && light) {
              if (r[order[l],j] < r[order[l],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l+1
            }
            if (light) { unres2 <- c(1,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) {
          if (r[order[length(order)],j] < r[order[length(order)],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[length(order)],j] > r[order[length(order)],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN }
          else {
            l <- length(order)-1
            light <- 1
            while (l >= 1 && light) {
              if (r[order[l],j] < r[order[l],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l-1
            }
            if (light) { unres2 <- c(2,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else stop("There are too many ties in the ordering process - please, consider using another ordering algorithm.")
      }
      else {
        temp <- sample(length(which(r[i,][x]==r[i,which.min(r[i,][x])])),2)
        m1 <- temp[1]; m2 <- temp[2]
        j <- which(r[i,][x]==r[i,which.min(r[i,][x])])[m1]
        k <- which(r[i,][x]==r[i,which.min(r[i,][x])])[m2]
        if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) { order <- c(j,order,k); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) { order <- c(k,order,j); x[j] <- NaN; x[k] <- NaN }
        else if (r[order[1],j] < r[order[length(order)],j] && r[order[1],k] < r[order[length(order)],k]) {
          if (r[order[1],j] < r[order[1],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[1],j] > r[order[1],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN }
          else {
            l <- 2
            light <- 1
            while (l <= length(order) && light) {
              if (r[order[l],j] < r[order[l],k]) { order <- c(k,j,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(j,k,order); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l+1
            }
            if (light) { unres2 <- c(1,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else if (r[order[1],j] > r[order[length(order)],j] && r[order[1],k] > r[order[length(order)],k]) {
          if (r[order[length(order)],j] < r[order[length(order)],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN }
          else if (r[order[length(order)],j] > r[order[length(order)],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN }
          else {
            l <- length(order)-1
            light <- 1
            while (l >= 1 && light) {
              if (r[order[l],j] < r[order[l],k]) { order <- c(order,j,k); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else if (r[order[l],j] > r[order[l],k]) { order <- c(order,k,j); x[j] <- NaN; x[k] <- NaN; light <- 0 }
              else l <- l-1
            }
            if (light) { unres2 <- c(2,j,k); x[j] <- NaN; x[k] <- NaN }
          }
        }
        else stop("There are too many ties in the ordination process - please, consider using another ordering algorithm.")
      }
    }
  }
  return(avoid.reverse(order))
}

##Continuity Index
Cindex <- function (order,r) {
  markers <- dim(r)[1]
  CI <- 0
  for (i in 1:(markers-1))
      for (j in (i+1):markers)
        CI <- CI + r[order[i],order[j]]/((j-i)^2)
  return (CI)
}

## end of file
