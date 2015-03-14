#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: rf.graph.table.R                                              #
# Contains: rf.graph.table, plotFunction.out and draw.rf.inter        #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 03/05/2009                                           #
# Last update: 09/25/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

rf.graph.table <- function(input.seq, scale=1, axis.cex=1, main=NULL, inter=TRUE) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(input.seq$seq.phases[1] == -1 || input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like)) stop("You must estimate parameters before running 'rf.graph.table' ") 
  ## making a list with necessary information  
  size <- length(input.seq$seq.num)
  max.rf <- 0.5
  LOD<-list(CC=matrix(NA,size,size),
            CR=matrix(NA,size,size),
            RC=matrix(NA,size,size),
            RR=matrix(NA,size,size))
  mat <- matrix(0,size,size)
  mat.rf <- matrix(0,size,size)
  for (i in 2:size) {
    for (j in 1:(i-1)) {
      big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[j])
      small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[j])
      current <- get(input.seq$twopt, pos=1)$analysis[acum(big-2)+small,,]
      probab <- which(current[,2]>(max(current[,2]-0.005)) & current[,2]<=max(current[,2]))
      goodness <- numeric(4)
      for (a in probab) {
        if (current[a,1] <= max.rf) goodness[a] <- 1
        else goodness[a] <- 0
      }
      goodness[-probab] <- 0
      phase <- which(goodness==1)
      if (length(phase)==0){
        mat[j,i] <- 0.0
        mat.rf[j,i] <- mat.rf[i,j] <- mat[i,j] <- 0.5
      }
      else {
        mat[j,i] <- max(current[phase,2])
        mat.rf[j,i] <- mat.rf[i,j] <- mat[i,j] <- min(current[phase,1])
      }
      LOD$CC[i,j] <- LOD$CC[j,i] <- current[1,2]
      LOD$CR[i,j] <- LOD$CR[j,i] <- current[2,2]
      LOD$RC[i,j] <- LOD$RC[j,i] <- current[3,2]
      LOD$RR[i,j] <- LOD$RR[j,i] <- current[4,2]
    }
  }

  ##Scaling the LODs to print them properly
  mat[row(mat) > col(mat) & mat > 0.5] <- 0.5 ## if there are recombinations greater than 0.5 (for numrical convergence problems), assuming 0.5
  mat[row(mat) < col(mat)][mat[row(mat) < col(mat)] < 10E-2]<-10E-2
  log.LOD<-log(mat[row(mat) < col(mat)]) + abs(min(log(mat[row(mat) < col(mat)]), na.rm=TRUE))
  mat[row(mat) < col(mat)] <- 0.5 - log.LOD/(2*max(log.LOD, na.rm=TRUE))
  diag(mat)<-NA

  ##Write multipoint estimates
  for (i in 1:(size-1)){
    mat[i+1,i] <- input.seq$seq.rf[i]
    mat.rf[i+1,i] <- mat.rf[i,i+1] <- input.seq$seq.rf[i]
  }
  colnames(mat) <- get(input.seq$twopt, pos=1)$marnames[input.seq$seq.num]

  ##Write NAs in two-point recombination fractions between markers of type D1 and D2
  types <- get(input.seq$data.name, pos=1)$segr.type[input.seq$seq.num]
  which.D1D2<-outer((substr(types, 1,2)=="D1"),(substr(types, 1,2)=="D2"))
  which.D1D2<-which.D1D2+t(which.D1D2)
  which.D1D2[which.D1D2==1]<-NA
  which.D1D2[which.D1D2==0]<-1
  diag.si<-rbind(1:(ncol(which.D1D2)-1),2:ncol(which.D1D2))
  for(i in 1:(ncol(which.D1D2)-1)) which.D1D2[diag.si[1,i],diag.si[2,i]] <- which.D1D2[diag.si[2,i],diag.si[1,i]] <- 1
  mat<-mat*which.D1D2
  missing<-100*apply(get(input.seq$data.name, pos=1)$geno[,input.seq$seq.num],2, function(x) sum(x==0))/get(input.seq$data.name, pos=1)$n.ind
  
  ##info.graph contains all information necessary to plot the graphics 
  info.graph <- list(mat=mat,
                     mat.rf=mat.rf,
                     seq.num=input.seq$seq.num,
                     n=ncol(mat),
                     names=colnames(mat),
                     types=types,
                     LOD=LOD,
                     missing=missing,
                     data.name=input.seq$data.name)
  if (inter==FALSE)
    plotFunction.out(info.graph=info.graph, cex=axis.cex, main=main)

  else draw.rf.inter(info.graph=info.graph,scale=scale,cex=axis.cex)
}

##This function plots the recombination fraction NOT using interactive Tcl-Tk interface
plotFunction.out <- function(info.graph, cex, main)
  {
##layout(matrix(1:2,1,2), widths = c(10,2))
params <- par(bg="white", plt=c(0.1,.95, 0.1, 0.9),xpd=TRUE)
y.adj<- .8/(-2+info.graph$n*2)
##print(info.graph$mat)
if(is.null(main)) main="LOD (above diag.) and Recombination Fraction Matrix"
image(info.graph$mat, axes=FALSE,
      col=rainbow(n=500, start=min(info.graph$mat,na.rm=TRUE)*1.3, end=max(info.graph$mat,na.rm=TRUE)*1.3))
title(main = main, cex.main=.75)
x<-seq(from=0, to=1, length=info.graph$n)
text(x=x, y=y.adj+rep(-diff(x)[1],info.graph$n),info.graph$seq.num, srt=90, cex=cex, adj=1)
text(y=x, x=y.adj+rep(-diff(x)[1],info.graph$n),info.graph$seq.num, cex=cex, adj=1)
##par(cex.axis=.8)
##plot(x=rep(1,51), y=seq(from=0, to=.5, by=0.01), xlim=c(1,1.5), cex=2, col=rainbow(51, start=0, end=.65), pch=15, axes=FALSE, xlab="", ylab="")
##axis(4, pos=c(1.1,0))
par(params)
}

##This function plots the recombination fraction using interactive Tcl-Tk interface
draw.rf.inter<-function(info.graph, scale, cex){
  ##Identical above, but used to Tcl-Tk plot
  plotFunction <- function()
    {
     ## layout(matrix(1:2,1,2), widths = c(10,2))
      params <- par(bg="white", plt=c(0.1,.95, 0.1, 0.9),xpd=TRUE)
      y.adj<- .8/(-2+info.graph$n*2)
      image(info.graph$mat, axes=FALSE, main="",
            col=rainbow(n=500, start=min(info.graph$mat,na.rm=TRUE)*1.3, end=max(info.graph$mat,na.rm=TRUE)*1.3))
      title(main = "LOD (above diag.) and Recombination Fraction Matrix", cex.main=(2+scale)/4, line=.6)
      x<-seq(from=0, to=1, length=info.graph$n)
      text(x=x, y=y.adj+rep(-diff(x)[1],info.graph$n),info.graph$seq.num, srt=90, cex=cex, adj=1)
      text(y=x, x=y.adj+rep(-diff(x)[1],info.graph$n),info.graph$seq.num, cex=cex, adj=1)
     ## par(cex.axis=.8)
     ## plot(x=rep(1,51), y=seq(from=0, to=.5, by=0.01), xlim=c(1,1.5), cex=scale*1.2, col=rainbow(51, start=0, end=.65), pch=15, axes=FALSE, xlab="", ylab="")
     ## axis(4, pos=c(1.1,0))
      par(params)
    }
  
  ## Getting the mouse coords with TclTk
  OnLeftClick <- function(x,y)  {
    ##parPlotSize<-c(0.10/1.2, 0.95/1.2, 0.10, 0.90)
    parPlotSize<-c(0.10, 0.95, 0.10, 0.90)
    usrCoords<-rep(c(-1/(-2+info.graph$n*2), 1+1/(-2+info.graph$n*2)),2)
    xClick <- x
    yClick <- y
    xCoords<-(seq(from=0,to=1,by=info.graph$n))
    yCoords<-(seq(from=0,to=1,by=info.graph$n))
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
    xMin <- parPlotSize[1] * width
    xMax <- parPlotSize[2] * width
    yMin <- parPlotSize[3] * height
    yMax <- parPlotSize[4] * height
    rangeX <- usrCoords[2] - usrCoords[1]
    rangeY <- usrCoords[4] - usrCoords[3]
    imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
    imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
    xClick <- as.numeric(xClick)+0.5
    yClick <- as.numeric(yClick)+0.5
    yClick <- height - yClick
    xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
    yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
    y<-seq(from=-1/(2*(info.graph$n-1)), to=1+(1/(2*(info.graph$n-1))), by=1/(info.graph$n-1))
    ##Printing information about selected markers
    x.n<-sum(xPlotCoord > y)
    y.n<-sum(yPlotCoord > y)
    mkx.n<-info.graph$seq.num[x.n]
    mky.n<-info.graph$seq.num[y.n]
    if(mkx.n==mky.n){
      if(class(get(info.graph$data.name, pos=1))=="outcross"){
        ##information for message box (class 'outcross')
        msg <- paste("Marker name: \n    ", info.graph$names[x.n],
                     "\n\nMarker number:\n    ", mkx.n,
                     "\n\nMarker type: \n    ", info.graph$types[x.n],
                     "\n\n", format(info.graph$missing[x.n], digits=2), "% of missing data for this marker",
                     sep="")
      }
      else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap"){
        ##getting type of marker ('another classes')
        if(info.graph$types[x.n]=="B3.7")  mkt<-"AA : AB : BB (1:2:1) "
        else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == 1)  mkt<-" Not BB : BB (3:1) "
        else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == -1)  mkt<-" Not AA : AA (3:1) "
        else if(info.graph$types[x.n]=="D1.10")  mkt<-"AA : AB (1:1)"
        else stop ("invalid type of marker at marker ", info.graph$names[x.n])
        ##information for message box
        msg <- paste("Marker name: \n    ", info.graph$names[x.n],
                     "\n\nMarker number:\n    ", mkx.n,
                     "\n\nMarker type: \n    ", mkt,
                     "\n\n", format(info.graph$missing[x.n], digits=2), "% of missing data for this marker",
                     sep="")
        
      }
      else if (class(get(info.graph$data.name, pos=1))=="riself.onemap" || class(get(info.graph$data.name, pos=1))=="risib.onemap"){
        ##getting type of marker ('another classes')
        if(info.graph$types[x.n]=="D1.10")  mkt<-"AA : BB (1:1)"
        else stop ("invalid type of marker at marker ", info.graph$names[x.n])
        ##information for message box
        msg <- paste("Marker name: \n    ", info.graph$names[x.n],
                     "\n\nMarker number:\n    ", mkx.n,
                     "\n\nMarker type: \n    ", mkt,
                     "\n\n", format(info.graph$missing[x.n], digits=2), "% of missing data for this marker",
                     sep="")
        
      }
      else stop("invalid type of cross")
      mbval<- tkmessageBox(title="Labeling Marker",message=msg,type="ok",icon="question")
    }
    else{
      if(x.n==(y.n+1) || y.n==(x.n+1)){
        if(class(get(info.graph$data.name, pos=1))=="outcross"){
          ##information for message box (class 'outcross')
          msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                       "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                       "\n\nMarker types: \n    ", info.graph$types[y.n], " and ", info.graph$types[x.n],
                       "\n\nMultipoint recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                       "\n\nLOD-Scores of the linkage phases:",
                       "\n    CC: ", format(info.graph$LOD$CC[x.n, y.n], digits=2),
                       "\n    CR: ", format(info.graph$LOD$CR[x.n, y.n], digits=2),
                       "\n    RC: ", format(info.graph$LOD$RC[x.n, y.n], digits=2),
                       "\n    RR: ", format(info.graph$LOD$RR[x.n, y.n], digits=2),
                       sep="")
        }
        else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap"){
           ##getting type of marker
          if(info.graph$types[x.n]=="B3.7")  mktx<-"AA : AB : BB (1:2:1) "
          else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == 1)  mktx<-" Not BB : BB (3:1) "
          else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == -1)  mktx<-" Not AA : AA (3:1) "
          else if(info.graph$types[x.n]=="D1.10")  mktx<-"AA : AB (1:1)"
          else stop ("invalid type of marker at marker ", info.graph$names[x.n])
          if(info.graph$types[y.n]=="B3.7")  mkty<-"AA : AB : BB (1:2:1) "
          else if(info.graph$types[y.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mky.n] == 1)  mkty<-" Not BB : BB (3:1) "
          else if(info.graph$types[y.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mky.n] == -1)  mkty<-" Not AA : AA (3:1) "
          else if(info.graph$types[y.n]=="D1.10")  mkty<-"AA : AB (1:1)"
          else stop ("invalid type of marker at marker ", info.graph$names[y.n])
          
          if(info.graph$LOD$CC[x.n, y.n] > info.graph$LOD$RR[x.n, y.n]) LODScore<-info.graph$LOD$CC[x.n, y.n]
          else LODScore<-info.graph$LOD$RR[x.n, y.n]
          
          ##information for message box (another classes)
          msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                       "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                       "\n\nMarker types: \n    ", mkty, " \n    and \n    ", mktx,
                       "\n\nMultipoint recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                       "\n\nLOD-Score: \n    ", format(LODScore, digits=2),
                       sep="")
        }

        else if (class(get(info.graph$data.name, pos=1))=="riself.onemap" || class(get(info.graph$data.name, pos=1))=="risib.onemap"){
          ##getting type of marker
          if(info.graph$types[x.n]=="D1.10")  mktx<-"AA : BB (1:1)"
          else stop ("invalid type of marker at marker ", info.graph$names[x.n])
          if(info.graph$types[y.n]=="D1.10")  mkty<-"AA : BB (1:1)"
          else stop ("invalid type of marker at marker ", info.graph$names[y.n])
           
          if(info.graph$LOD$CC[x.n, y.n] > info.graph$LOD$RR[x.n, y.n]) LODScore<-info.graph$LOD$CC[x.n, y.n]
          else LODScore<-info.graph$LOD$RR[x.n, y.n]
          
          ##information for message box (another classes)
          msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                       "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                       "\n\nMarker types: \n    ", mkty, " \n    and \n    ", mktx,
                       "\n\nMultipoint recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                       "\n\nLOD-Score: \n    ", format(LODScore, digits=2),
                       sep="")
        }
        else stop("invalid type of cross") 
      }
      else{
        if((substr(info.graph$types[x.n],1,2)=="D1" && substr(info.graph$types[y.n],1,2)=="D2") || (substr(info.graph$types[x.n],1,2)=="D2" && substr(info.graph$types[y.n],1,2)=="D1")){
          msg <- paste("Marker names: \n    ", info.graph$names[x.n], "\n    and \n    ", info.graph$names[y.n],
                       "\n\nMarkers of type \n    ", info.graph$types[x.n], " and ", info.graph$types[y.n],
                       "\n\nImpossible to estimate recombination fraction via two-point.",
                       sep="")
        }
        else{
          if(class(get(info.graph$data.name, pos=1))=="outcross"){
            ##information for message box (class 'outcross')
            msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                         "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                         "\n\nMarker types: \n    ", info.graph$types[y.n], " and ", info.graph$types[x.n], 
                         "\n\nTwo-point recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                         "\n\nLOD-Scores of the linkage phases:",
                         "\n    CC: ", format(info.graph$LOD$CC[x.n, y.n], digits=2),
                         "\n    CR: ", format(info.graph$LOD$CR[x.n, y.n], digits=2),
                         "\n    RC: ", format(info.graph$LOD$RC[x.n, y.n], digits=2),
                         "\n    RR: ", format(info.graph$LOD$RR[x.n, y.n], digits=2),
                         sep="")
          }
          else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap"){
            ##getting type of marker
            if(info.graph$types[x.n]=="B3.7")  mktx<-"AA : AB : BB (1:2:1) "
            else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == 1)  mktx<-" Not BB : BB (3:1) "
            else if(info.graph$types[x.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mkx.n] == -1)  mktx<-" Not AA : AA (3:1) "
            else if(info.graph$types[x.n]=="D1.10")  mktx<-"AA : AB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[x.n])
            if(info.graph$types[y.n]=="B3.7")  mkty<-"AA : AB : BB (1:2:1) "
            else if(info.graph$types[y.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mky.n] == 1)  mkty<-" Not BB : BB (3:1) "
            else if(info.graph$types[y.n]=="C.8"  &  get(info.graph$data.name, pos=1)$phase[mky.n] == -1)  mkty<-" Not AA : AA (3:1) "
            else if(info.graph$types[y.n]=="D1.10")  mkty<-"AA : AB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[y.n])
            
            if(info.graph$LOD$CC[x.n, y.n] > info.graph$LOD$RR[x.n, y.n]) LODScore<-info.graph$LOD$CC[x.n, y.n]
            else LODScore<-info.graph$LOD$RR[x.n, y.n]
            
            ##information for message box (another classes)
            msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                         "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                         "\n\nMarker types: \n    ", mkty, "\n    and \n    ",  mktx, 
                         "\n\nTwo-point recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                         "\n\nLOD-Score: \n    ", format(LODScore, digits=2),
                         sep="")
            
          }
          else if (class(get(info.graph$data.name, pos=1))=="riself.onemap" || class(get(info.graph$data.name, pos=1))=="risib.onemap"){
            ##getting type of marker
            if(info.graph$types[x.n]=="D1.10")  mktx<-"AA : BB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[x.n])
            if(info.graph$types[y.n]=="D1.10")  mkty<-"AA : BB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[y.n])
            
            if(info.graph$LOD$CC[x.n, y.n] > info.graph$LOD$RR[x.n, y.n]) LODScore<-info.graph$LOD$CC[x.n, y.n]
            else LODScore<-info.graph$LOD$RR[x.n, y.n]
            
            ##information for message box (another classes)
            msg <- paste("Marker names: \n    ", info.graph$names[y.n], "\n    and \n    ", info.graph$names[x.n],
                         "\n\nMarker numbers:\n    ",mky.n," and ",mkx.n,
                         "\n\nMarker types: \n    ", mkty, "\n    and \n    ",  mktx, 
                         "\n\nTwo-point recombination fraction:\n    rf = ",format(info.graph$mat.rf[x.n, y.n], digits=4),
                         "\n\nLOD-Score: \n    ", format(LODScore, digits=2),
                         sep="")
            
          }
          else stop("invalid type of cross") 

        }
      }
      mbval<- tkmessageBox(title="Labeling recombination fraction",message=msg,type="ok",icon="question")
    }   
  }
  
  tt<-tktoplevel()
  tkwm.title(tt,"Click on a pixel to label it")
  img <- tkrplot(tt,fun=plotFunction,hscale=scale,vscale=scale)
  tkgrid(img)
  tkbind(img, "<Button-1>", OnLeftClick)
  tkconfigure(img,cursor="hand2")
}


