#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: try.seq.R                                                     #
# Contains: try.seq, print.try, draw.try                              #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 11/29/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function tries to position a marker in a given linkage map
## CHECK IF THE ord$rf MATRIX HAS SOME NaN VALUE AND ISSUE WARNING (CHANGE FOR 'BAD' VALUE)
try.seq <-
function(input.seq,mrk,tol=10E-2,draw.try=FALSE,pos= NULL,verbose=FALSE) {
  # checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(input.seq$seq.phases[1] == -1 || input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like)) stop("You must run 'compare' or 'map' before the 'try.seq' function")
  if(mrk > get(input.seq$data.name, pos=1)$n.mar) stop(deparse(substitute(mrk))," exceeds the number of markers in object ", input.seq$data.name)

  # allocate variables
  rf.init <- vector("list",length(input.seq$seq.num))
  phase.init <- vector("list",length(input.seq$seq.num))
  ord <- list(list(matrix(NA,16,length(input.seq$seq.num)),
                 matrix(NA,16,length(input.seq$seq.num)),
                 rep(-Inf,16)))
  names(ord[[1]]) <- c("rf","phase","like")
  ord <- rep(ord,length(input.seq$seq.num)+1)
  best.seq <- rep(-Inf,length(input.seq$seq.num)+1)
      
  ### positioning before the given sequence
  # get two-point information
  list.init <- phases(make.seq(get(input.seq$twopt),c(mrk,input.seq$seq.num[1]),twopt=input.seq$twopt))
  rf.init[[1]] <- list.init$rf.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) rf.init[[j+1]] <- input.seq$seq.rf[j]
  phase.init[[1]] <- list.init$phase.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) phase.init[[j+1]] <- input.seq$seq.phases[j]
  Ph.Init <- comb.ger(phase.init)
  Rf.Init <- comb.ger(rf.init)
  mark.max<-max(nchar(colnames(get(input.seq$data.name, pos=1)$geno)))
  num.max<-nchar(ncol(get(input.seq$data.name, pos=1)$geno))
  
  # create first order
  try.ord <- c(mrk,input.seq$seq.num)
  if(verbose) cat("TRY", 1,": ", c(mrk,input.seq$seq.num),"\n")
  else cat(format(mrk,width=num.max) , "-->", format(colnames(get(input.seq$data.name, pos=1)$geno)[mrk], width=mark.max), ": .")
  flush.console()
  
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem.amb.ph(M=Ph.Init, w=input.seq, seq.num=c(mrk,input.seq$seq.num))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(class(Ph.Init) == "numeric" || class(Ph.Init) == "integer"){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,c(mrk,input.seq$seq.num)],
                           type=get(input.seq$data.name, pos=1)$segr.type.num[c(mrk,input.seq$seq.num)],
                           phase=Ph.Init[j,],
                           rec=Rf.Init[j,],
                           verbose=FALSE,
                           tol=tol)
    ord[[1]]$rf[j,] <- final.map$rf
    ord[[1]]$phase[j,] <- Ph.Init[j,]
    ord[[1]]$like[j] <- final.map$loglike
    best.seq[1] <- max(best.seq[1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[1]]$like, decreasing=TRUE)
  ord[[1]]$rf <- ord[[1]]$rf[ord.ind,]
  ord[[1]]$phase <- ord[[1]]$phase[ord.ind,]
  ord[[1]]$like <- ord[[1]]$like[ord.ind]
  
  ### positioning between markers of the given sequence
  for(i in 1:(length(input.seq$seq.num)-1)) {
    # get two-point information
    list.init <- phases(make.seq(get(input.seq$twopt),c(input.seq$seq.num[i],mrk,input.seq$seq.num[i+1]),twopt=input.seq$twopt))
    if(i!=1) {
      for(k in 1:(i-1)) {
        rf.init[[k]] <- input.seq$seq.rf[k]
        phase.init[[k]] <- input.seq$seq.phases[k]
      }
    }
    rf.init[[i]] <- list.init$rf.init[[1]]
    phase.init[[i]] <- list.init$phase.init[[1]]
    rf.init[[i+1]] <- list.init$rf.init[[3]]
    phase.init[[i+1]] <- list.init$phase.init[[3]]
    if(i!=(length(input.seq$seq.num)-1)) {
      for(k in (i+2):length(input.seq$seq.num)) {
        rf.init[[k]] <- input.seq$seq.rf[k-1]
        phase.init[[k]] <- input.seq$seq.phases[k-1]
      }
    }   
    Ph.Init <- comb.ger(phase.init)
    Rf.Init <- comb.ger(rf.init)
 
	# create intermediate orders
    try.ord <- rbind(try.ord,c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
    if(verbose) cat("TRY", i+1,": ",c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]) ,"\n")
    else cat(".")
    flush.console()
    
    if(nrow(Ph.Init)>1){
      ##Removing ambigous phases
      rm.ab<-rem.amb.ph(M=Ph.Init, w=input.seq, seq.num=c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
      Ph.Init <- Ph.Init[rm.ab,]
      Rf.Init <- Rf.Init[rm.ab,]
      if(class(Ph.Init) == "numeric" || class(Ph.Init) == "integer"){
        Ph.Init<-matrix(Ph.Init,nrow=1)
        Rf.Init<-matrix(Rf.Init,nrow=1)
      }
    }
    ## estimate parameters for all possible linkage phases for the current order
    for(j in 1:nrow(Ph.Init)) {
      final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)])],
							 type=get(input.seq$data.name, pos=1)$segr.type.num[c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)])],
                             phase=Ph.Init[j,],
                             rec=Rf.Init[j,],
                             verbose=FALSE,
                             tol=tol)
      ord[[i+1]]$rf[j,] <- final.map$rf
      ord[[i+1]]$phase[j,] <- Ph.Init[j,]
      ord[[i+1]]$like[j] <- final.map$loglike
      best.seq[i+1] <- max(best.seq[i+1],final.map$loglike)
    }
	# sort linkage phases by log-likelihood
    ord.ind <- order(ord[[i+1]]$like, decreasing=TRUE)
    ord[[i+1]]$rf <- ord[[i+1]]$rf[ord.ind,]
    ord[[i+1]]$phase <- ord[[i+1]]$phase[ord.ind,]
    ord[[i+1]]$like <- ord[[i+1]]$like[ord.ind] 
  }
  
  ### positioning after the given sequence
  # get two-point information
  list.init <- phases(make.seq(get(input.seq$twopt),c(input.seq$seq.num[length(input.seq$seq.num)],mrk),twopt=input.seq$twopt))
  rf.init[[(length(input.seq$seq.num))]] <- list.init$rf.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) rf.init[[j]] <- input.seq$seq.rf[j]
  phase.init[[(length(input.seq$seq.num))]] <- list.init$phase.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) phase.init[[j]] <- input.seq$seq.phases[j]
  Ph.Init <- comb.ger(phase.init)
  Rf.Init <- comb.ger(rf.init)
  
  # create last order
  try.ord <- rbind(try.ord,c(input.seq$seq.num,mrk))
  if(verbose) cat("TRY",length(input.seq$seq.num)+1,": ", c(input.seq$seq.num,mrk) ,"\n")
  else cat(".\n")
  flush.console()
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem.amb.ph(M=Ph.Init, w=input.seq, seq.num=c(input.seq$seq.num,mrk))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(class(Ph.Init) == "numeric" || class(Ph.Init)=="integer"){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,c(input.seq$seq.num,mrk)],
                           type=get(input.seq$data.name, pos=1)$segr.type.num[c(input.seq$seq.num,mrk)],
                           phase=Ph.Init[j,],
                           rec=Rf.Init[j,],
                           verbose=FALSE,
                           tol=tol)
    ord[[length(input.seq$seq.num)+1]]$rf[j,] <- final.map$rf
    ord[[length(input.seq$seq.num)+1]]$phase[j,] <- Ph.Init[j,]
    ord[[length(input.seq$seq.num)+1]]$like[j] <- final.map$loglike
    best.seq[length(input.seq$seq.num)+1] <- max(best.seq[length(input.seq$seq.num)+1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[length(input.seq$seq.num)+1]]$like, decreasing=TRUE)
  ord[[length(input.seq$seq.num)+1]]$rf <- ord[[length(input.seq$seq.num)+1]]$rf[ord.ind,]
  ord[[length(input.seq$seq.num)+1]]$phase <- ord[[length(input.seq$seq.num)+1]]$phase[ord.ind,]
  ord[[length(input.seq$seq.num)+1]]$like <- ord[[length(input.seq$seq.num)+1]]$like[ord.ind]
  
  # calculate LOD-Scores (best linkage phase combination for each position)
  LOD <- (best.seq-max(best.seq))/log(10)
   if(draw.try==TRUE){
    draw.try(input.seq, structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "try"), pos=pos)
  }
  structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "try")
}



# print method for object class 'try'
print.try <- function(x,j=NULL,...) {
  phases.char <- c("CC","CR","RC","RR")
  marker <- x$try.ord[1,1]
  
  if(is.null(j)) {
    # general summary
    seq <- format(x$try.ord[1,-1], scientific = FALSE)
    size1 <- max(nchar(seq))
    seq.pr <- format(seq,width=size1)
    size2 <- max(nchar(formatC(x$LOD,format="f",digits=2)))
    LOD.pr <- formatC(round(x$LOD,2),format="f",digits=2,width=size2)
    LOD.pr[which(LOD.pr=="  -0.0")] <- "   0.0"
    
    cat("\nLOD scores correspond to the best linkage phase combination\nfor each position\n")
    cat("\nThe symbol \"*\" outside the box indicates that more than one\nlinkage phase is possible for the corresponding position\n")
    cat(paste("\n\n\t\t  Marker tested: ",marker,"\n\n",sep=""))
    cat("\t\t  Markers",rep("",size1+size2-4),"LOD\n")
    cat(paste("\t\t",paste(rep("=",size1+size2+13),collapse=""),"\n",sep=""))
    cat("\t\t|",rep("",size1+size2+10),"|\n")
    cat("\t\t|",rep("",size1+8),LOD.pr[1]," |")
    ifelse(max(which(x$ord[[1]]$like != -Inf)) != 1,pr <- paste("  ",1,"  *\n",sep=""),pr <- paste("  ",1,"  \n",sep=""))
    cat(pr)
    for(i in 1:length(seq.pr)) {
      cat("\t\t| ",seq.pr[i],rep("",size2+8),"|","\n")
      cat("\t\t|",rep("",size1+8),LOD.pr[i+1]," |")
      ifelse(max(which(x$ord[[i+1]]$like != -Inf)) != 1,pr <- paste("  ",i+1,"  *\n",sep=""),pr <- paste("  ",i+1,"  \n",sep=""))
      cat(pr)
    }
    cat("\t\t|",rep("",size1+size2+10),"|\n")
    cat(paste("\t\t",paste(rep("=",size1+size2+13),collapse=""),"\n",sep=""))
  }
  else {
    # detailed output for a given position
    seq <- format(x$try.ord[j,], scientific = FALSE)
    size1 <- max(nchar(seq))
    seq.pr <- format(seq,width=size1)
    n.phase <- max(which(x$ord[[j]]$like != -Inf))
    max.like <- -Inf
    for(i in 1:length(x$ord)) max.like <- max(max.like,max(x$ord[[i]]$like[1:n.phase]))
    LOD <- round((x$ord[[j]]$like[1:n.phase]-max.like)/log(10),1)
    LOD.pr <- formatC(LOD,format="f",digits=1,width=6)
    LOD.pr[which(LOD.pr=="  -0.0")] <- "   0.0"
    nest.LOD <- round((x$ord[[j]]$like[1:n.phase]-max(x$ord[[j]]$like[1:n.phase]))/log(10),1)
    nest.LOD.pr <- formatC(nest.LOD,format="f",digits=1,width=6)
    nest.LOD.pr[which(nest.LOD.pr=="  -0.0")] <- "   0.0"
    
    cat("\nLOD is the overall LOD score (among all orders)\n")
    cat("\nNEST.LOD is the LOD score within the order\n")
    cat(paste("\nMarker tested: ",marker,"\n",sep=""))
    
    cat(paste(rep("-",max(2,size1)+5+7*(n.phase)),collapse=""),"\n")
    cat("|",rep("",max(2,size1)+2),rep("|     ",n.phase),"|\n")
    cat("|",seq.pr[1],rep("",max(3-size1,0)),rep("|     ",n.phase),"|\n|")
    for(i in 2:length(seq.pr)) {
      cat(rep("",max(2,size1)+2))
      for(k in 1:n.phase) {
        cat("  | ",phases.char[x$ord[[j]]$phase[k,i-1]])
      }
      cat("  |",rep("",max(3-size1,0)),"\n")
      cat("|",seq.pr[i],rep("",max(3-size1,0)),rep("|     ",n.phase),"|\n|")
    }
    cat("",rep("",max(2,size1)+2),rep("|     ",n.phase),"|\n")
    cat("|",rep("-",max(2,size1)+3+7*(n.phase)),"|","\n",sep="")
    if (size1 > 2) cat("| LOD ",rep(" ",size1-2), sep="")
    else cat("| LOD ")
    for(k in 1:n.phase) {
      cat(paste("|",LOD.pr[k],sep=""))
    }
    cat("|\n")
    cat("|",rep("-",max(2,size1)+3+7*(n.phase)),"|","\n",sep="")
    if (size1 > 2) cat("|NEST.",rep(" ",size1-2), sep="")
    else cat("|NEST.")
    cat(rep("|     ",n.phase),"|\n")
    if (size1 > 2) cat("| LOD ",rep(" ",size1-2), sep="")
    else cat("| LOD ")
    for(k in 1:n.phase) {
      cat(paste("|",nest.LOD.pr[k],sep=""))
    }
    cat("|\n")
    cat(paste(rep("-",max(2,size1)+5+7*(n.phase)),collapse=""),"\n")
  }
}



draw.try<-function(base.input, try.input, pos=NULL){
  layout(matrix(c(2,1,2,3),2,2), heights = c(1,2.5))
  base.dist<-cumsum(c(0, kosambi(base.input$seq.rf)))
  base.input.len<-length(base.dist)
  try.dist<-c(-1, base.dist[-base.input.len]+kosambi(base.input$seq.rf)/2, base.dist[base.input.len]+1)
  op<-par(mar=c(5,7,7,2), cex=.75)
  plot(x=try.dist, try.input$LOD, typ="l", xlab="Frame", ylab="LOD", axes=FALSE)
  abline(v=try.dist, lty=2, lwd=.5)
  op<-par(mar=c(5,7,7,2), cex=.75, xpd=TRUE)
  text(x=try.dist, y=rep(max(abs(try.input$LOD))/5,length(try.dist)), labels=1:length(try.dist), cex=.7)
  text(x=try.dist[1]-(max(try.dist)/40), y=max(abs(try.input$LOD))/5 ,"Position",  adj=c(1,0.5))
  text(x=try.dist[1]-(max(try.dist)/40), y=max(abs(try.input$LOD))/10 ,"Distance",  adj=c(1,0.05))
  axis(2)
  axis(3, at=round(base.dist,1), lwd.ticks = .5, cex.axis=.75, las=2)
  par(op)
  if(is.null(pos)){
    op<-par(xpd=TRUE)
    points(try.dist[which.max(try.input$LOD)],0, pch=17, col=2, cex=1.5)
    new.map<-make.seq(try.input,which.max(try.input$LOD))
    new.dist<-cumsum(c(0, kosambi(new.map$seq.rf)))
    new.dist.len<-length(new.dist)
    plot(x=new.dist, rep(1,new.dist.len), pch="|", xlab="New Genetic Map", ylab="", axes=FALSE, type="n", main=paste("Adding marker ",try.input$try.ord[1,1]," (", colnames(get(try.input$data.name,pos=1)$geno)[try.input$try.ord[1,1]],")", sep=""))
    axis(1, at=round(new.dist,1), lwd.ticks = .75, cex.axis=.75, las=2)
    text(new.dist, y=rep(1,length(new.dist)), labels=new.map$seq.num, cex=.7)
    points(new.dist[which.max(try.input$LOD)],1.5, col=2, cex=1.5, pch=25, bg = 2)
    text(x=new.dist[1]-(max(new.dist)/40), y=1 ,"Markers",  adj=c(1,0.5))
    text(x=new.dist[1]-(max(new.dist)/40), y=0 ,"Distance",  adj=c(1,0.2))
    par(op)
  }
  else{
    op<-par(xpd=TRUE)
    points(try.dist[pos],0, pch=17, col=2, cex=1.5)
    new.map<-make.seq(try.input,pos)
    new.dist<-cumsum(c(0, kosambi(new.map$seq.rf)))
    new.dist.len<-length(new.dist)
    plot(x=new.dist, rep(1,new.dist.len), xlab="New Genetic Map", ylab="", axes=FALSE, type="n", main=paste("Adding marker ",try.input$try.ord[1,1]," (", colnames(get(try.input$data.name,pos=1)$geno)[try.input$try.ord[1,1]],")", sep=""))
    axis(1, at=round(new.dist,1), lwd.ticks = .75, cex.axis=.75, las=2)
    text(new.dist, y=rep(1,length(new.dist)), labels=new.map$seq.num, cex=.7)
    points(new.dist[pos], 1.5, col=2, cex=1.5, pch=25, bg = 2)
    text(x=new.dist[1]-(max(new.dist)/40), y=1 ,"Markers",  adj=c(1,0.5))
    text(x=new.dist[1]-(max(new.dist)/40), y=0 ,"Distance",  adj=c(1,0.2))
    par(op)
  }
  rf.graph.table(new.map, inter=FALSE, axis.cex = .75, main="")
  title(main = "LOD (above diag.) and Recombination Fraction Matrix", cex.main=.9, line=15.4)
}

# end of file
