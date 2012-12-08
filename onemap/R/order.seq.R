#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: order.seq.R                                                   #
# Contains: order.seq, print.order, draw.order                        #
#                                                                     #
# Written by Gabriel R A Margarido & Marcelo Mollinari                #
# copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/04/2011                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## This function automates linkage map construction in two steps:
## first, it applies the 'compare' algorithm to a subset of markers;
## second, it adds markers sequentially with the 'try' function
order.seq <- function(input.seq, n.init=5, subset.search=c("twopt", "sample"),
                      subset.n.try=30, subset.THRES=3, twopt.alg= c("rec", "rcd", "ser", "ug"),
                      THRES=3, touchdown=FALSE, draw.try=FALSE, wait=0, tol=10E-2) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(n.init < 2) stop("'n.init' must be greater than or equal to 2")
  if(!is.logical(touchdown)) stop("'touchdown' must be logical")
  if(!touchdown && THRES <= 10E-10) stop("Threshold must be greater than 0 if 'touchdown' is FALSE")
  if(touchdown && THRES <= (1 + 10E-10)) stop("Threshold must be greater than 1 if 'touchdown' is TRUE")
  if(wait < 0){
    warning("'wait' shouldn't be < 0!")
    wait<-0
  }
  if(draw.try == FALSE) wait<-0
  if(length(input.seq$seq.num) <= n.init) {
    ## in this case, only the 'compare' function is used
    cat("   Length of sequence ",deparse(substitute(input.seq))," is less than n.init \n   Returning the best order using compare function:\n")
    ifelse(length(input.seq$seq.num) == 2, seq.ord <- map(input.seq,tol=10E-5), seq.ord <- make.seq(compare(input.seq=input.seq,tol=10E-5),1))
    seq.ord<-map(seq.ord, tol=10E-5)
    structure(list(ord=seq.ord, mrk.unpos=NULL, LOD.unpos=NULL, THRES=THRES,
                   ord.all=seq.ord, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "order")
  }
  else {
    ## here, the complete algorithm will be applied
    if(class(get(input.seq$data.name, pos=1)) == "f2.onemap") FLAG <- "f2"
    else if (class(get(input.seq$data.name, pos=1)) == "bc.onemap" || class(get(input.seq$data.name, pos=1)) == "riself.onemap" || class(get(input.seq$data.name, pos=1)) == "risib.onemap") FLAG <- "bc"
    else if (class(get(input.seq$data.name, pos=1)) == "outcross") FLAG <- "outcross"
    else stop("Invalid cross type\n")
    cross.type<-substr(class(get(input.seq$data.name, pos=1)), 1, nchar(class(get(input.seq$data.name, pos=1)))-7)
    ## select the order in which markers will be added
    if(FLAG == "bc" || FLAG == "f2"){
      subset.search <- match.arg(subset.search)
      if(subset.search == "twopt"){
        cat("\nCross type: ", cross.type, "\nChoosing initial subset using 'two-point' approach\n")
        twopt.alg <- match.arg(twopt.alg)
        tpt.type <- switch(EXPR=twopt.alg,
                           'rec'={
                             seq.rec<-record(input.seq=input.seq,tol=0.1)$seq.num ##ordering using RECORD algorithm
                             seq.init<-seq.rec[unique(round(seq(from=1, to=length(seq.rec), length.out=n.init)))] ##taking equally spaced markers
                           },
                           'rcd'={
                             seq.rcd<-rcd(input.seq=input.seq,tol=0.1)$seq.num ##ordering using RCD algorithm
                             seq.init<-seq.rcd[unique(round(seq(from=1, to=length(seq.rcd), length.out=n.init)))] ##taking equally spaced markers
                           },
                           'ser'={
                             seq.ser<-seriation(input.seq=input.seq,tol=0.1)$seq.num ##ordering using SERIATION algorithm
                             seq.init<-seq.ser[unique(round(seq(from=1, to=length(seq.ser), length.out=n.init)))] ##taking equally spaced markers
                           },
                           'ug'={
                             seq.ug<-ug(input.seq=input.seq,tol=0.1)$seq.num ##ordering using UG algorithm
                             seq.init<-seq.ug[unique(round(seq(from=1, to=length(seq.ug), length.out=n.init)))] ##taking equally spaced markers
                           })
        ##if(is.null(tpt.type)) stop("Invalid two point method")
        seq.rest <- input.seq$seq.num[-pmatch(seq.init, input.seq$seq.num)] ##the rest of the markers 
        seq.mis <- apply(as.matrix(get(input.seq$data.name, pos=1)$geno[,seq.rest]), 2, function(x) sum(x==0)) ##checking missing markers for the rest
        names(seq.mis)<-colnames(get(input.seq$data.name, pos=1)$geno)[seq.rest]
        
        if(FLAG == "bc") {
          rest.ord <- pmatch(names(seq.mis), colnames(get(input.seq$data.name, pos=1)$geno))
          seq.work <- pmatch(c(seq.init,rest.ord), input.seq$seq.num)
        }
        else if(FLAG == "f2"){
          seq.type <- get(input.seq$data.name, pos=1)$segr.type.num[seq.rest] ##checking maker type (4: co-dominant; 5: dominant)
          names(seq.type) <- names(seq.mis)
          tp.ord <- sort(seq.type) ##sorting by type, automatically co-dominant markers (4) come first
          rest.ord <- c(sort(seq.mis[pmatch(names(tp.ord)[tp.ord==4],names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord==5],names(seq.mis))]))##within each type of marker, sort by missing 
          rest <- pmatch(names(rest.ord), colnames(get(input.seq$data.name, pos=1)$geno)) ##matching names with the positions on the raw data
          seq.work <- pmatch(c(seq.init,rest), input.seq$seq.num) ##sequence to work with
        }
        else stop("Invalid cross type\n")
      }
      else if(subset.search=="sample"){
        cat("\nCross type: ", cross.type, "\nChoosing initial subset using the 'sample' approach\n")
        LOD.test <- i <- 0
        while(abs(LOD.test) < abs(subset.THRES) && i < subset.n.try){
          smp.seq <- make.seq(get(input.seq$twopt), sample(input.seq$seq.num, size=n.init), twopt=input.seq$twopt)
          res.test <- compare(smp.seq)
          LOD.test <- res.test$best.ord.LOD[2]
          i < -i+1
        }     
        if(abs(LOD.test) >= abs(subset.THRES)){
          seq.init <- res.test$best.ord[1,] ##best order based on 'compare'
          seq.rest <- input.seq$seq.num[-pmatch(seq.init, input.seq$seq.num)] ##the rest of the markers 
          seq.mis <- apply(as.matrix(get(input.seq$data.name, pos=1)$geno[,seq.rest]), 2, function(x) sum(x==0)) ##checking missing markers for the rest
           names(seq.mis)<-colnames(get(input.seq$data.name, pos=1)$geno)[seq.rest]
          if(FLAG == "bc"){
            rest.ord <- pmatch(names(seq.mis), colnames(get(input.seq$data.name, pos=1)$geno))
            seq.work <- pmatch(c(seq.init,rest.ord), input.seq$seq.num)
          }
          else if(FLAG == "f2"){
            seq.type <- get(input.seq$data.name, pos=1)$segr.type.num[seq.rest] ##checking maker type (4: co-dominant; 5: dominant)
            names(seq.type) <- names(seq.mis)
            tp.ord <- sort(seq.type) ##sorting by type, automatically co-dominant markers (4) come first
            rest.ord <- c(sort(seq.mis[pmatch(names(tp.ord)[tp.ord==4],names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord==5],names(seq.mis))]))##within each type of marker, sort by missing 
            rest <- pmatch(names(rest.ord),colnames(get(input.seq$data.name, pos=1)$geno)) ##matching names with the positions on the raw data
            seq.work <- pmatch(c(seq.init,rest), input.seq$seq.num) ##sequence to work with
          }
          else stop("Invalid cross type\n")
        }
        else stop("Cannot find any subset using 'subset.n.try'=", subset.n.try, " and 'subset.THRES'= ", subset.THRES,"\n")
      }
      ## else stop("Invalid subset search\n")
    }
    else if(FLAG == "outcross") {
      cat("\nCross type: outcross\nUsing segregation types of the markers to choose initial subset\n")
      segregation.types <- get(input.seq$data.name, pos=1)$segr.type.num[input.seq$seq.num]
      if(sum(segregation.types == 7) > sum(segregation.types == 6)) segregation.types[segregation.types == 6] <- 8 ## if there are more markers of type D2 than D1, try to map those first
      seq.work <- order(segregation.types)
      seq.init <- input.seq$seq.num[seq.work[1:n.init]] 
    }
    else stop("Invalid cross type")   
    ##apply the 'compare' step to the subset of initial markers
    seq.ord <- compare(input.seq=make.seq(get(input.seq$twopt), seq.init, twopt=input.seq$twopt), n.best=50)
    
    ## 'try' to map remaining markers
    input.seq2 <- make.seq(seq.ord,1)
    cat ("\n\nRunning try algorithm\n")
    for (i in (n.init+1):length(input.seq$seq.num)){
      time.elapsed<-system.time(seq.ord <- try.seq(input.seq2,input.seq$seq.num[seq.work[i]],tol=tol, draw.try=draw.try))[3]
      if(time.elapsed < wait)
        Sys.sleep(wait - time.elapsed)
      if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < -THRES))
        input.seq2 <- make.seq(seq.ord,which.max(seq.ord$LOD))
    }
    
    ## markers that do not meet the threshold remain unpositioned
    mrk.unpos <- input.seq$seq.num[which(is.na(match(input.seq$seq.num, input.seq2$seq.num)))]
    LOD.unpos <- NULL
    cat("\nLOD threshold =",THRES,"\n\nPositioned markers:", input.seq2$seq.num, "\n\n")
    cat("Markers not placed on the map:", mrk.unpos, "\n")
    
    if(touchdown && length(mrk.unpos) > 0) {
      ## here, a second round of the 'try' algorithm is performed, if requested
      cat("\n\n\nTrying to map remaining markers with LOD threshold ",THRES-1,"\n")
      for (i in mrk.unpos) {
        time.elapsed<-system.time(seq.ord <- try.seq(input.seq2,i,tol=tol, draw.try=draw.try))[3]
        if(time.elapsed < wait)
          Sys.sleep(wait - time.elapsed)
        if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < (-THRES+1)))
          input.seq2 <- make.seq(seq.ord,which.max(seq.ord$LOD))
      }
      
      ## markers that do not meet this second threshold still remain unpositioned 
      mrk.unpos <- input.seq$seq.num[which(is.na(match(input.seq$seq.num, input.seq2$seq.num)))]
      cat("\nLOD threshold =",THRES-1,"\n\nPositioned markers:", input.seq2$seq.num, "\n\n")
      cat("Markers not placed on the map:", mrk.unpos, "\n")
    }
    
    if(length(mrk.unpos) > 0) {
      ## LOD-Scores are calculated for each position, for each unmapped marker, if any
      LOD.unpos <- matrix(NA,length(mrk.unpos),(length(input.seq2$seq.num)+1))
      j <- 1
      cat("\n\nCalculating LOD-Scores\n")
      for (i in mrk.unpos){
        LOD.unpos[j,] <- try.seq(input.seq=input.seq2,mrk=i,tol=tol)$LOD
        j <- j+1
      }
    }
    else mrk.unpos <- NULL

    ## to end the algorithm, possibly remaining markers are 'forced' into the map
    input.seq3 <- input.seq2
    if(!is.null(mrk.unpos)) {
      cat("\n\nPlacing remaining marker(s) at most likely position\n")
      
      ## these markers are added from the least to the most doubtful
      which.order <- order(apply(LOD.unpos,1,function(x) max(x[-which(x==0)[1]])))
      
      for (i in mrk.unpos[which.order]) {        
        time.elapsed<-system.time(seq.ord <- try.seq(input.seq3,i,tol,draw.try=draw.try))[3]
        if(time.elapsed < wait)
          Sys.sleep(wait - time.elapsed)  
        input.seq3 <- make.seq(seq.ord,which(seq.ord$LOD==0)[sample(sum(seq.ord$LOD==0))[1]])
      }
    }
    cat("\nEstimating final genetic map using tol = 10E-5.\n\n")
    input.seq2<-map(input.seq2, tol=10E-5)
    input.seq3<-map(input.seq3, tol=10E-5)
    if (draw.try == TRUE)
      draw.order(input.seq3)
    structure(list(ord=input.seq2, mrk.unpos=mrk.unpos, LOD.unpos=LOD.unpos, THRES=THRES,
                   ord.all=input.seq3, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "order")
  }
}

print.order <- function(x,...) {
  cat("\nBest sequence found.")
  ## print the 'safe' order
  print(x$ord)
  if(!is.null(x$mrk.unpos)) {
    ## print LOD-Score information for unpositioned markers
    cat("\n\nThe following markers could not be uniquely positioned.\n")
    cat("Printing most likely positions for each unpositioned marker:\n")
    
    size1 <- max(3,max(nchar(x$mrk.unpos)))
    mrk.unpos.pr <- format(x$mrk.unpos,width=size1)
    size2 <- max(nchar(x$ord$seq.num))
    seq.pr <- format(x$ord$seq.num,width=size2)
    
######limit <- (x$THRES-2)/2 ## previously used limit
    
    cat("\n")
    cat(paste(rep("-",size2+4+length(mrk.unpos.pr)*(size1+3)),collapse=""),"\n")
    cat("| ",rep("",size2),"|")
###### MAYBE WE SHOULD PUT A LIMIT TO THE NUMBER OF UNPOSITIONED MARKERS
    for(j in 1:length(mrk.unpos.pr)) {
      cat(rep("",max(0,3-size1)+1),mrk.unpos.pr[j],"|")
    }
    cat("\n")
    cat(paste("|",paste(rep("-",size2+2),collapse=""),"|",sep=""))
    cat(paste(rep(paste(paste(rep("-",size1+2),collapse=""),"|",sep=""),length(mrk.unpos.pr)),collapse=""),"\n")
    cat("| ",rep("",size2),"|")
    for(j in 1:length(x$mrk.unpos)) {
      if(x$LOD.unpos[j,1] > -0.0001) cat(rep("",max(0,3-size1)+1),"*** |")
      else if(x$LOD.unpos[j,1] > -1.0) cat(rep("",max(0,3-size1)+1),"**  |")
      else if(x$LOD.unpos[j,1] > -2.0) cat(rep("",max(0,3-size1)+1),"*   |")
      else cat(rep("",max(0,3-size1)+1),"    |")
    }
    cat("\n")
    for(i in 1:length(seq.pr)) {
      cat("|",seq.pr[i],"|")
      cat(paste(rep(paste(paste(rep(" ",size1+2),collapse=""),"|",sep=""),length(mrk.unpos.pr)),collapse=""),"\n")
      cat(paste("|",paste(rep(" ",size2+2),collapse=""),"|",sep=""))
      for(j in 1:length(x$mrk.unpos)) {
        if(x$LOD.unpos[j,i+1] > -0.0001) cat(rep("",max(0,3-size1)+1),"*** |")
        else if(x$LOD.unpos[j,i+1] > -1.0) cat(rep("",max(0,3-size1)+1),"**  |")
        else if(x$LOD.unpos[j,i+1] > -2.0) cat(rep("",max(0,3-size1)+1),"*   |")
        else cat(rep("",max(0,3-size1)+1),"    |")
      }
      cat("\n")
    }
    cat(paste(rep("-",size2+4+length(mrk.unpos.pr)*(size1+3)),collapse=""),"\n")
    cat("\n")
    cat("'***' indicates the most likely position(s) (LOD = 0.0)\n\n")
    cat("'**' indicates very likely positions (LOD > -1.0)\n\n")
    cat("'*' indicates likely positions (LOD > -2.0)\n\n")
  }
}

draw.order<-function(map.input){
  layout(matrix(c(1,2),2,1), heights = c(.8,2.5))
  op<-par(mar=c(6,5,4,2), cex=.75, xpd=TRUE)
  new.dist<-cumsum(c(0, kosambi(map.input$seq.rf)))
  new.dist.len<-length(new.dist)
  plot(x=new.dist, rep(1,new.dist.len), xlab="", ylab="", axes=FALSE, type="n", main="Final Genetic Map")
  text(new.dist, y=rep(1,length(new.dist)), labels=map.input$seq.num, cex=.7)
  axis(1, at=round(new.dist,1), lwd.ticks = .75, cex.axis=.7, las=2)
  text(x=new.dist[1]-(max(new.dist)/40), y=1 ,"Markers",  adj=c(1,0.5))
  text(x=new.dist[1]-(max(new.dist)/40), y=0 ,"Distance",  adj=c(1,0.2))
  par(op)
  rf.graph.table(map.input, inter=FALSE, axis.cex = .75, main="")
  title(main = "LOD (above diag.) and Recombination Fraction Matrix", cex.main=.9, line=15.4)
}
## end of file
