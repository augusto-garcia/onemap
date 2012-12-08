#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: compare.R                                                     #
# Contains: compare, print.compare                                    #
#                                                                     #
# Written by Gabriel R A Margarido & Marcelo Mollinari                #
# copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## This function evaluates all n!/2 possible orders for n markers
compare <- function(input.seq,n.best=50,tol=10E-4,verbose=FALSE) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(sQuote(deparse(substitute(input.seq)))," is not an object of class 'sequence'")
  if(length(input.seq$seq.num) > 5) cat("WARNING: this operation may take a VERY long time\n")
  flush.console()
  if(length(input.seq$seq.num) > 10) {
    cat("\nIt is not wisely trying to use 'compare' with more than 10 markers \n")
    ANSWER <- readline("Are you sure you want to proceed? [y or n]\n")
    while(substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 1) != "y")
      ANSWER <- readline("\nPlease answer: 'y' or 'n' \n")
    if (substr(ANSWER, 1, 1) == "n") stop("Execution stopped!")
  }
  if(length(input.seq$seq.num) == 2) return(map(input.seq, tol=tol)) ## nothing to be done for 2 markers
  else {
    ## allocating variables
    rf.init <- vector("list",length(input.seq$seq.num)-1)
    phase.init <- vector("list",length(input.seq$seq.num)-1)
    best.ord <- matrix(NA,(n.best+1),length(input.seq$seq.num))
    best.ord.rf <- matrix(NA,(n.best+1),length(input.seq$seq.num)-1)
    best.ord.phase <- matrix(NA,(n.best+1),length(input.seq$seq.num)-1)
    best.ord.like <- best.ord.LOD <- rep(-Inf,(n.best+1))
    
    ## 'phases' gathers information from two-point analyses
    list.init <- phases(input.seq)
    
    ## 'perm.pars' generates all n!/2 orders
    all.ord <- perm.pars(input.seq$seq.num)
    cat("\nComparing",nrow(all.ord),"orders:     \n\n")
    if (verbose){ 
      for(i in 1:nrow(all.ord)){
        ## print output for each order
        cat("Order", i, ":", all.ord[i,], "\n")
        flush.console()  
        ## get initial values for the HMM
        all.match <- match(all.ord[i,],input.seq$seq.num)
        for(j in 1:(length(input.seq$seq.num)-1)){
          if(all.match[j] > all.match[j+1]){
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j]-2)+all.match[j+1]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j]-2)+all.match[j+1]]]
          }
          else {
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j+1]-2)+all.match[j]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j+1]-2)+all.match[j]]]
          }
        }    
        Ph.Init <- comb.ger(phase.init)
        Rf.Init <- comb.ger(rf.init)
        if(nrow(Ph.Init)>1){
          ##Removing ambigous phases
          rm.ab<-rem.amb.ph(M=Ph.Init, w=input.seq, seq.num=all.ord[i,])
          Ph.Init <- Ph.Init[rm.ab,]
          Rf.Init <- Rf.Init[rm.ab,]
          if(class(Ph.Init)=="integer"){
            Ph.Init<-matrix(Ph.Init,nrow=1)
            Rf.Init<-matrix(Rf.Init,nrow=1)
          }
        }
        for(j in 1:nrow(Ph.Init)){
          ## estimate parameters
          final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,all.ord[i,]],
                                 type=get(input.seq$data.name, pos=1)$segr.type.num[all.ord[i,]],
                                 phase=Ph.Init[j,],
                                 rec=Rf.Init[j,],
                                 verbose=FALSE,
                                 tol=tol)
          best.ord[(n.best+1),] <- all.ord[i,]
          best.ord.rf[(n.best+1),] <- final.map$rf
          best.ord.phase[(n.best+1),] <- Ph.Init[j,]
          best.ord.like[(n.best+1)] <- final.map$loglike
          
          ## arrange orders according to the likelihood
          like.order <- order(best.ord.like, decreasing=TRUE)
          best.ord <- best.ord[like.order,]
          best.ord.rf <- best.ord.rf[like.order,]
          best.ord.phase <- best.ord.phase[like.order,]
          best.ord.like <- sort(best.ord.like, decreasing=TRUE)
        }
      }
    }
    else{
      count <- 0
      pb <- txtProgressBar(style=3)
      setTxtProgressBar(pb, 0)
      
      ## nc<-NA
      ## out.pr <- seq(from=1,to=nrow(all.ord), length.out=20)
      cat("    ")
      for(i in 1:nrow(all.ord)){
        
      ##                                    # print output for each order
      ##    if (sum(i == round(out.pr))){
      ##      cat(rep("\b",nchar(nc)+1),sep="")
      ##      nc<-round(i*100/nrow(all.ord))
      ##      cat(nc,"%", sep="")
      ##      flush.console()
      ##    }      
        ## get initial values for the HMM
        all.match <- match(all.ord[i,],input.seq$seq.num)
        for(j in 1:(length(input.seq$seq.num)-1)){
          if(all.match[j] > all.match[j+1]){
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j]-2)+all.match[j+1]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j]-2)+all.match[j+1]]]
          }
          else {
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j+1]-2)+all.match[j]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j+1]-2)+all.match[j]]]
          }
        }    
        Ph.Init <- comb.ger(phase.init)
        Rf.Init <- comb.ger(rf.init)
        if(nrow(Ph.Init)>1){
          ##Removing ambigous phases
          rm.ab<-rem.amb.ph(M=Ph.Init, w=input.seq, seq.num=all.ord[i,])
          Ph.Init <- Ph.Init[rm.ab,]
          Rf.Init <- Rf.Init[rm.ab,]
          if(class(Ph.Init)=="integer"){
            Ph.Init<-matrix(Ph.Init,nrow=1)
            Rf.Init<-matrix(Rf.Init,nrow=1)
          }
        }
        for(j in 1:nrow(Ph.Init)){
          ## estimate parameters
          final.map <- est.map.c(geno=get(input.seq$data.name, pos=1)$geno[,all.ord[i,]],
                                 type=get(input.seq$data.name, pos=1)$segr.type.num[all.ord[i,]],
                                 phase=Ph.Init[j,],
                                 rec=Rf.Init[j,],
                                 verbose=FALSE,
                                 tol=tol)
          best.ord[(n.best+1),] <- all.ord[i,]
          best.ord.rf[(n.best+1),] <- final.map$rf
          best.ord.phase[(n.best+1),] <- Ph.Init[j,]
          best.ord.like[(n.best+1)] <- final.map$loglike
          
          ## arrange orders according to the likelihood
          like.order <- order(best.ord.like, decreasing=TRUE)
          best.ord <- best.ord[like.order,]
          best.ord.rf <- best.ord.rf[like.order,]
          best.ord.phase <- best.ord.phase[like.order,]
          best.ord.like <- sort(best.ord.like, decreasing=TRUE)       
        }
        count<-count+1
        setTxtProgressBar(pb, count/nrow(all.ord))
      }
      close(pb)
    }
    ## cat("\nFinished\n\n")
    cat("\n")
    best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),4)
    
    structure(list(best.ord = best.ord, best.ord.rf = best.ord.rf,
                   best.ord.phase = best.ord.phase, best.ord.like = best.ord.like,
                   best.ord.LOD = best.ord.LOD, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "compare")
  }
}

## print method for object class 'compare'
print.compare <-
  function(x,...) {
    FLAG<-0
    if(class(get(x$data.name, pos=1)) != "outcross") FLAG<-1   
    phases.char <- c("CC","CR","RC","RR")
    n.ord <- max(which(head(x$best.ord.LOD,-1) != -Inf))
    unique.orders <- unique(x$best.ord[1:n.ord,])
    n.ord.nest <- nrow(unique.orders)
    phases.nested <- vector("list",n.ord.nest)
    LOD <- vector("list",n.ord.nest)
    
    for (i in 1:n.ord.nest) {
      same.order <- which(apply(x$best.ord[1:n.ord,],1,function(x) all(x==unique.orders[i,])))
      ifelse(length(same.order)==1,phases.nested[[i]] <- t(as.matrix(x$best.ord.phase[same.order,])),phases.nested[[i]] <- x$best.ord.phase[same.order,])
      LOD[[i]] <- x$best.ord.LOD[same.order]
    }
    skip <- c(nchar(n.ord.nest),max(nchar(unique.orders[1,])+2))
    cat("\nNumber of orders:",n.ord,"\n")
    if(FLAG==0){ ## outcrossing
      leng.print <- nchar(paste("order ",format(n.ord.nest,width=skip[1]),":  ",paste(format(unique.orders[1,],width=skip[2]),collapse=""),"     ",format(11.11,digits=2,format="f",width=6),"     ",format(11.11,digits=2,format="f",width=6),"\n",sep=""))
       cat(paste("Best ",n.ord.nest," unique orders",paste(rep(" ",leng.print-37),collapse=""),"LOD    Nested LOD","\n",sep=""))
      cat(paste(rep("-",leng.print),collapse=""),"\n")
    }
    else if(FLAG==1){ ## other
      leng.print <- nchar(paste("order ",format(n.ord.nest,width=skip[1]),":  ",paste(format(unique.orders[1,],width=skip[2]),collapse=""),"     ",format(11.11,digits=2,format="f",width=6),"\n",sep=""))
      cat(paste("Best ",n.ord.nest," unique orders",paste(rep(" ",leng.print-25),collapse=""),"LOD","\n",sep=""))
      cat(paste(rep("-",leng.print),collapse=""),"\n")
    }
    else stop ("Should not get here!")
    
    if(FLAG==0){ ## outcrossing
      for (i in 1:n.ord.nest) {
        cat(paste("order ",format(i,width=skip[1]),":  ",paste(format(unique.orders[i,],width=skip[2]),collapse=""),"\n",sep=""))
        for (j in 1:dim(phases.nested[[i]])[1]) {
          cat(paste("\t",paste(rep(" ",1+skip[1]+skip[2]),collapse=""),paste(format(phases.char[phases.nested[[i]][j,]],width=skip[2]),collapse=""),"     ",formatC(round(LOD[[i]][j],2),digits=2,format="f",width=6),"     ",formatC(round(LOD[[i]][j]-LOD[[i]][1],2),digits=2,format="f",width=6),"\n",sep=""))
        }
        cat(paste(rep("-",leng.print),collapse=""))
        cat("\n")
      }   
    }
    else if(FLAG==1){ ## other
      for (i in 1:n.ord.nest) {
        cat(paste("order ",format(i,width=skip[1]),":  ",paste(format(unique.orders[i,],width=skip[2]),collapse=""), "     ",formatC(round(LOD[[i]][1],2),digits=2,format="f",width=6),  "\n",sep=""))
      }
      cat(paste(rep("-",leng.print),collapse=""))
      cat("\n")
    }
    else stop ("Should not get here!") 
  }
## end of file
