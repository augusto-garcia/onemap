#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: rcd.R                                                         #
# Contains: rcd                                                       #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/29/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

rcd <-function(input.seq, LOD=0,max.rf=0.5,tol=10E-5) {
  # checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
    not an object of class 'sequence'")
  
  markers <- length(input.seq$seq.num)

  ## create reconmbination fraction matrix 
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

  # x defines the non-positioned markers
  x <- 1:markers

  # the group starts with the closest two markers
  i <- which(r==r[which.min(r)])
  i <- i[sample(length(i),1)]

  # 'first' and 'last' are the markers on the edges of the ordered group
  first <- ceiling(i/markers)
  last <- i-((first-1)*markers)
  
  # the two markers are set next to each other  
  order <- c(first,last)

  # markers already ordered are marked as NaN
  x[first] <- NaN
  x[last] <- NaN

  # extending the chain
  while (length(order) < markers) {
    # get the markers closest to the left end of the group
    j <- which(r[first,][x]==r[first,which.min(r[first,][x])])
    # randomly choose one of them
    if (length(j) > 1) j <- j[sample(length(j),1)]

    # get the markers closest to the right end of the group
    k <- which(r[last,][x]==r[last,which.min(r[last,][x])])
    # randomly choose one of them
    if (length(k) > 1) k <- k[sample(length(k),1)]

    if (r[first,j] < r[last,k]) {  # place a marker at the left side
      order <- c(j,order)
      x[j] <- NaN
      first <- j
    }
    else if (r[first,j] > r[last,k]) {  # or at the right side
      order <- c(order,k)
      x[k] <- NaN
      last <- k } 
    else {  # if the distance is the same, randomly choose one side
      rand <- sample(2,1) 
      if (rand==1) { order <- c(j,order); x[j] <- NaN; first <- j }
      else { order <- c(order,k); x[k] <- NaN; last <- k }
    }
  }
  # end of chain
  cat("\norder obtained using RCD algorithm:\n\n", input.seq$seq.num[avoid.reverse(order)], "\n\ncalculating multipoint map using tol = ", tol, ".\n\n")
  map(make.seq(get(input.seq$twopt),input.seq$seq.num[avoid.reverse(order)],twopt=input.seq$twopt), tol=tol)
}
# end of file





