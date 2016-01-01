#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: rf.graph.table.R                                              ##
## Contains: rf.graph.table, plotFunction.out and draw.rf.inter        ##
##                                                                     ##
## Written by Marcelo Mollinari                                        ##
## copyright (c) 2009, Marcelo Mollinari                               ##
##                                                                     ##
## First version: 03/05/2009                                           ##
## Last update: 09/25/2009                                             ##
## Description was modified by Augusto Garcia on 2015/07/25            ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' Plots pairwise recombination fractions and LOD Scores in a heatmap
##'
##' Plots a matrix of pairwise recombination fractions (under the diagonal) and
##' LOD Scores (upper the diagonal) using a color scale. Any value of the
##' matrix can be easily accessed using an interactive Tcl-Tk interface,
##' helping users to check for possible problems.
##'
##' The color scale varies from red (small distances or big LODs) to dark blue.
##' When clicking on a cell, a dialog box is displayed with some information
##' about corresponding markers for that cell (line \eqn{\times} column). They are:
##' \eqn{i}) the name of the markers; \eqn{ii}) the number of
##' the markers on the data set; \eqn{iii}) the segregation types; \eqn{iv})
##' the recombination fraction between the markers and \eqn{v}) the LOD-Score
##' for each possible linkage phase calculated via two-point analysis. For
##' neighbor markers, the multipoint recombination fraction is printed;
##' otherwise, the two-point recombination fraction is printed. For markers of
##' type \code{D1} and \code{D2}, it's impossible to calculate recombination
##' fraction via two-point analysis and, therefore the corresponding cell will
##' be empty. For cells on the diagonal of the matrix, the name, the number and
##' the type of the marker are printed, as well as the percentage of missing
##' data for that marker.
##'
##' @param input.seq an object of class \code{sequence} with a predefined
##' order.
##' @param scale controls the plot size. If \code{inter == FALSE} this value is
##' not used.
##' @param axis.cex the magnification to be used for axis annotation.
##' @param main the title for no interactive plot, i.e. it is only used if
##' \code{inter == FALSE}.
##' @param inter logical. If \code{TRUE}, an interactive graphic is plotted.
##' Otherwise, a default graphic device is used.
##' @param mrk.names logical. If \code{TRUE}, displays the names of the markers.
##' @param colorkey logical. If \code{TRUE}, a colorkey is plotted
##'     along horizontal axis, indicating recombination fraction, and
##'     along vertical axis, indicating the LOD Score.
##' @author Marcelo Mollinari, \email{mmollina@@gmail.com}
##' @keywords utilities
##' @examples
##'
##' ##outcross example
##'   data(example.out)
##'   twopt <- rf.2pts(example.out)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rcd <- rcd(LG1)
##'   rf.graph.table(LG1.rcd, inter=FALSE)
##' 
##' \dontrun{
##'   ##Now, using interactive Tcl-Tk
##'   rf.graph.table(LG1.rcd, scale=2, inter=TRUE)
##'
##'   ##F2 example
##'   data(fake.f2.onemap)
##'   twopt <- rf.2pts(fake.f2.onemap)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'
##'   ##"pre-allocate" an empty list of length groups$n.groups (3, in this case)
##'   maps.list<-vector("list", groups$n.groups)
##'
##'   for(i in 1:groups$n.groups){
##'     ##create linkage group i
##'     LG.cur <- make.seq(groups,i)
##'     ##ordering
##'     map.cur<-order.seq(LG.cur, subset.search = "sample")
##'     ##assign the map of the i-th group to the maps.list
##'     maps.list[[i]]<-make.seq(map.cur, "force")
##'   }
##'   ##Plot LOD/recombination fraction matrices for each group
##'   op <- par(mfrow = c(1, 3))
##'   for(i in 1:groups$n.groups)
##'     rf.graph.table(maps.list[[i]], axis.cex=.7, main=paste("Group", i),inter=FALSE)
##'   par(op)
##' }
##'

rf.graph.table <- function(input.seq, scale=1, axis.cex=1, main=NULL, inter=TRUE, mrk.names=FALSE, colorkey = TRUE) {
    ## checking for correct objects
    if(!any(class(input.seq)=="sequence"))
        stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")

    if(is(get(input.seq$data.name), "outcross"))
    {
        if(input.seq$seq.phases[1] == -1 || input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like))
            stop("You must estimate parameters before running 'rf.graph.table' ")
        ## making a list with necessary information
        n.mrk <- length(input.seq$seq.num)
        LOD <- lapply(get(input.seq$twopt)$analysis,
                      function(x, w){
                          m<-matrix(0, length(w), length(w))
                          for(i in 1:(length(w)-1)){
                              for(j in (i+1):length(w)){
                                  z<-sort(c(w[i],w[j]))
                                  m[j,i]<-m[i,j]<-x[z[1], z[2]]
                              }
                          }
                          return(m)
                      }, input.seq$seq.num
                      )
        mat<-t(get_mat_rf_out(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
        mat.rf<-t(get_mat_rf_out(input.seq, LOD=FALSE,  max.rf = 0.501, min.LOD = -0.1))
    }
    else
    {
        if(input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like))
            stop("You must estimate parameters before running 'rf.graph.table' ")
        ## making a list with necessary information
        n.mrk <- length(input.seq$seq.num)
        
        LOD<-matrix(0, length(input.seq$seq.num), length(input.seq$seq.num))
        for(i in 1:(length(input.seq$seq.num)-1)){
            for(j in (i+1):length(input.seq$seq.num)){
                z<-sort(c(input.seq$seq.num[i],input.seq$seq.num[j]))
                LOD[j,i]<-LOD[i,j]<-get(input.seq$twopt)$analysis[z[1], z[2]]
            }
        }
        mat<-t(get_mat_rf_in(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
        mat.rf<-t(get_mat_rf_in(input.seq, LOD=FALSE,  max.rf = 0.501, min.LOD = -0.1))
    }
    
    ##Scaling the LODs to print them properly 
    range.LOD<-range(as.dist(t(mat)), na.rm=TRUE)
    range.rf<-range(as.dist(mat), na.rm=TRUE)
    mat[row(mat) > col(mat) & mat > 0.5] <- 0.5 ## if there are recombinations greater than 0.5 (for numrical convergence problems), assuming 0.5
    mat[row(mat) < col(mat)][mat[row(mat) < col(mat)] < 10E-2]<-10E-2
    min.scale<-abs(min(log(mat[row(mat) < col(mat)]), na.rm=TRUE))
    log.LOD<-log(mat[row(mat) < col(mat)]) + min.scale
    max.scale<-2*max(log.LOD, na.rm=TRUE)  
    mat[row(mat) < col(mat)] <- 0.5 - log.LOD/max.scale
    diag(mat)<-NA
    range.scaled.LOD<-range(as.dist(t(mat)), na.rm=TRUE)
    
    ##Write multipoint estimates
    for (i in 1:(n.mrk-1)){
        mat[i+1,i] <- input.seq$seq.rf[i]
        mat.rf[i+1,i] <- mat.rf[i,i+1] <- input.seq$seq.rf[i]
    }
    colnames(mat) <- colnames(get(input.seq$data.name, pos=1)$geno)[input.seq$seq.num]
    
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
                       min.scale=min.scale,
                       max.scale=max.scale,
                       range.rf=range.rf,
                       range.LOD=range.LOD,
                       range.scaled.LOD=range.scaled.LOD,
                       missing=missing,
                       data.name=input.seq$data.name,
                       mrk.names=mrk.names,
                       colorkey=colorkey)
    if (inter==FALSE)
        plotFunction.out(info.graph=info.graph, cex=axis.cex, main=main)
    else draw.rf.inter(info.graph=info.graph,scale=scale,cex=axis.cex)
}

##This function plots the recombination fraction WITHOUT using interactive Tcl-Tk interface
plotFunction.out <- function(info.graph, cex, main){
    if(!info.graph$colorkey){
        params <- par(bg="white",
                      plt=c(0.1,.95, 0.1, 0.9),
                      xpd=TRUE)
    }
  else{
      layout(matrix(c(4,1,3,2),2,2), widths = c(10,2), heights=c(2,10))
      params <- par(bg="white",
                    plt=c(0.1,.95, 0.1, 0.95),
                    xpd=TRUE)
  }
    y.adj<- .8/(-2+info.graph$n*2)
    if(is.null(main)) main="LOD (above diag.) and Recombination Fraction Matrix"
    image(info.graph$mat, axes=FALSE,
          col=rev(tim.colors(200)),
          useRaster=TRUE)
    if(!info.graph$colorkey)
        title(main = main, cex.main=cex)
    x<-seq(from=0, to=1, length=info.graph$n)
    if (info.graph$mrk.names == TRUE) { 
        text(x = x, y = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$names,
             srt = 90, cex = cex, adj = 1)
        text(y = x, x = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$names, 
             cex = cex, adj = 1)
    } else {
        text(x = x, y = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$seq.num, 
             srt = 90, cex = cex, adj = 1)
        text(y = x, x = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$seq.num, 
             cex = cex, adj = 1)
    }
    if(info.graph$colorkey){
        par(cex.axis=cex, las=1)
        plot(x=rep(1,200), y=seq(from=info.graph$range.rf[1], to=info.graph$range.rf[2], length.out=200), xlim=c(1,1.5), cex=2,
             col=rev(tim.colors(200)),
             pch=15, axes=FALSE, xlab="", ylab="")
        mtext("rec. frac.",1, adj = 0)
        Axis(x=info.graph$range.rf,
             at=round(seq(from=info.graph$range.rf[1], to=info.graph$range.rf[2], length.out=10),2),
             side=4,
             pos=c(1.1,0))
        plot(0,type="n",axes=FALSE, xlab="", ylab="")#to fill the upper right corner
        plot(x=seq(from=info.graph$range.LOD[1], to=info.graph$range.LOD[2], length.out=200), y=rep(1,200),
             ylim=c(1,1.5),
             cex=2,
             col=tim.colors(200),
             pch=15,
             axes=FALSE,
             xlab="",
             ylab="")
        mtext("LOD",2, padj = 4.2)
        tk.pos<-seq(from=info.graph$range.scaled.LOD[1], to=info.graph$range.scaled.LOD[2], length.out=10)
        tk.lab<-rev(round(exp((.5-tk.pos)*info.graph$max.scale-info.graph$min.scale),1))
        Axis(x=info.graph$range.LOD,
             at=round(seq(from=info.graph$range.LOD[1], to=info.graph$range.LOD[2], length.out=10)),
             labels=tk.lab,
             side=3,
             pos=c(1.1,0))
        text(mean(info.graph$range.LOD), 1.5, main, cex=cex*1.5, pos=1)
    }
    par(params)
}

##This function plots the recombination fraction using interactive Tcl-Tk interface
draw.rf.inter<-function(info.graph, scale, cex){
    plotFunction <- function()
    {
        if(info.graph$colorkey){
            layout(matrix(c(4,1,3,2),2,2), widths = c(10,2), heights=c(2,10))
            params <- par(bg="white", plt=c(0.1,.95, 0.1, 0.95),xpd=TRUE)
        }
      else{
          params <- par(bg="white", plt=c(0.1,.95, 0.1, 0.9),xpd=TRUE)
      } 
        y.adj<- .8/(-2+info.graph$n*2)
        image(info.graph$mat, axes=FALSE,
              col=rev(tim.colors(200)), #rainbow(n=500, start=min(info.graph$mat, na.rm=TRUE)*1.3, end=max(info.graph$mat, na.rm=TRUE)*1.3)
              useRaster=TRUE)
        if(!info.graph$colorkey) title(main = "LOD (above diag.) and Recombination Fraction Matrix", cex.main=(4+scale)/4, line=.6)
        x<-seq(from=0, to=1, length=info.graph$n)
        if (info.graph$mrk.names == TRUE) { 
            text(x = x, y = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$names,
                 srt = 90, cex = cex, adj = 1)
            text(y = x, x = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$names, 
                 cex = cex, adj = 1)
        } else {
            text(x = x, y = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$seq.num, 
                 srt = 90, cex = cex, adj = 1)
            text(y = x, x = y.adj + rep(-diff(x)[1], info.graph$n), info.graph$seq.num, 
                 cex = cex, adj = 1)
        }
        if(info.graph$colorkey){
            par(cex.axis=cex, las=1)
            plot(x=rep(1,200), y=seq(from=info.graph$range.rf[1], to=info.graph$range.rf[2], length.out=200), xlim=c(1,1.5), cex=2,
                 col=rev(tim.colors(200)),#rainbow(n=200, start=info.graph$range.rf[1]*1.3, end=info.graph$range.rf[2]*1.3),
                 pch=15, axes=FALSE, xlab="", ylab="")
            mtext("rec. frac.",1, adj = 0)
            Axis(x=info.graph$range.rf,
                 at=round(seq(from=info.graph$range.rf[1], to=info.graph$range.rf[2], length.out=10),2),
                 side=4,
                 pos=c(1.1,0))
            plot(0,type="n",axes=FALSE, xlab="", ylab="")#to fill the upper right corner
            plot(x=seq(from=info.graph$range.LOD[1], to=info.graph$range.LOD[2], length.out=200), y=rep(1,200),
                 ylim=c(1,1.5),
                 ##xlim=c(0,info.graph$range.LOD[2]+5),
                 cex=2,
                 col=tim.colors(200),#rev(rainbow(n=200, start=info.graph$range.scaled.LOD[1]*1.3, end=info.graph$range.scaled.LOD[2]*1.3)),
                 pch=15,
                 axes=FALSE,
                 xlab="",
                 ylab="")
            mtext("LOD", 2, padj = 4.2)
            tk.pos<-seq(from=info.graph$range.scaled.LOD[1], to=info.graph$range.scaled.LOD[2], length.out=10)
            tk.lab<-rev(round(exp((.5-tk.pos)*info.graph$max.scale-info.graph$min.scale),1))
            Axis(x=info.graph$range.LOD,
                 at=round(seq(from=info.graph$range.LOD[1], to=info.graph$range.LOD[2], length.out=10)),
                 labels=tk.lab,
                 side=3,
                 pos=c(1.1,0))
            text(mean(info.graph$range.LOD), 1.5, "LOD (above diag.) and Recombination Fraction Matrix", cex=round((4+scale)/4), vfont=c("sans serif", "bold italic"), pos=1)
        }
        par(params)
    }
    ## Getting the mouse coords with TclTk
    OnLeftClick <- function(x,y)  {
        if(info.graph$colorkey) parPlotSize<-c(0.10/1.2, 0.95/1.2, 0.10/1.2, 0.95/1.2)
        else parPlotSize<-c(0.10, 0.95, 0.10, 0.90)
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
            else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap")
            {    
                ##getting type of marker ('another classes')
                if(info.graph$types[x.n]=="A.H.B")  mkt<-"AA : AB : BB (1:2:1) "
                else if(info.graph$types[x.n]=="D.B")  mkt<-" Not BB : BB (3:1) "
                else if(info.graph$types[x.n]=="C.A")  mkt<-" Not AA : AA (3:1) "
                else if(info.graph$types[x.n]=="A.H")  mkt<-"AA : AB (1:1)"
                else if(info.graph$types[x.n]=="M.X")  mkt<- "Mixed: Dominant & Co-dominant"
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
                if(info.graph$types[x.n]=="A.B")  mkt<-"AA : BB (1:1)"
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
                             "\n    CC: ", round(info.graph$LOD$CC[x.n, y.n], digits=1),
                             "\n    CR: ", round(info.graph$LOD$CR[x.n, y.n], digits=1),
                             "\n    RC: ", round(info.graph$LOD$RC[x.n, y.n], digits=1),
                             "\n    RR: ", round(info.graph$LOD$RR[x.n, y.n], digits=1),
                             sep="")
            }
            else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap"){
            ##getting type of marker
            if(info.graph$types[x.n]=="A.H.B")  mktx<-"AA : AB : BB (1:2:1) "
            else if(info.graph$types[x.n]=="D.B")  mktx<-" Not BB : BB (3:1) "
            else if(info.graph$types[x.n]=="C.A")  mktx<-" Not AA : AA (3:1) "
            else if(info.graph$types[x.n]=="A.H")  mktx<-"AA : AB (1:1)"
            else if(info.graph$types[x.n]=="M.X")  mktx<- "Mixed: Dominant & Co-dominant"
            else stop ("invalid type of marker at marker ", info.graph$names[x.n])
            if(info.graph$types[y.n]=="A.H.B")  mkty<-"AA : AB : BB (1:2:1) "
            else if(info.graph$types[y.n]=="D.B")  mkty<-" Not BB : BB (3:1) "
            else if(info.graph$types[y.n]=="C.A")  mkty<-" Not AA : AA (3:1) "
            else if(info.graph$types[y.n]=="A.H")  mkty<-"AA : AB (1:1)"
            else if(info.graph$types[y.n]=="M.X")  mkty<- "Mixed: Dominant & Co-dominant"
            else stop ("invalid type of marker at marker ", info.graph$names[y.n])
            LODScore<-info.graph$LOD[x.n, y.n]
                
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
            if(info.graph$types[x.n]=="A.B")  mktx<-"AA : BB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[x.n])
            if(info.graph$types[y.n]=="A.B")  mkty<-"AA : BB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[y.n])  
            LODScore<-info.graph$LOD[x.n, y.n]
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
          if((substr(info.graph$types[x.n],1,2)=="D1" && substr(info.graph$types[y.n],1,2)=="D2") ||
             (substr(info.graph$types[x.n],1,2)=="D2" && substr(info.graph$types[y.n],1,2)=="D1")){
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
                             "\n    CC: ", round(info.graph$LOD$CC[x.n, y.n], digits=1),
                             "\n    CR: ", round(info.graph$LOD$CR[x.n, y.n], digits=1),
                             "\n    RC: ", round(info.graph$LOD$RC[x.n, y.n], digits=1),
                             "\n    RR: ", round(info.graph$LOD$RR[x.n, y.n], digits=1),
                             sep="")
            }
          else if (class(get(info.graph$data.name, pos=1))=="f2.onemap" || class(get(info.graph$data.name, pos=1))=="bc.onemap"){
              ##getting type of marker
              if(info.graph$types[x.n]=="A.H.B")  mktx<-"AA : AB : BB (1:2:1) "
        else if(info.graph$types[x.n]=="D.B")  mktx<-" Not BB : BB (3:1) "
        else if(info.graph$types[x.n]=="C.A")  mktx<-" Not AA : AA (3:1) "
        else if(info.graph$types[x.n]=="A.H")  mktx<-"AA : AB (1:1)"
          else if(info.graph$types[x.n]=="M.X")  mktx<- "Mixed: Dominant & Co-dominant"
        else stop ("invalid type of marker at marker ", info.graph$names[x.n])
              
              if(info.graph$types[y.n]=="A.H.B")  mkty<-"AA : AB : BB (1:2:1) "
        else if(info.graph$types[y.n]=="D.B")  mkty<-" Not BB : BB (3:1) "
        else if(info.graph$types[y.n]=="C.A")  mkty<-" Not AA : AA (3:1) "
        else if(info.graph$types[y.n]=="A.H")  mkty<-"AA : AB (1:1)"
          else if(info.graph$types[y.n]=="M.X")  mkty<- "Mixed: Dominant & Co-dominant"
        else stop ("invalid type of marker at marker ", info.graph$names[y.n])

              LODScore<-info.graph$LOD[x.n, y.n]
              
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
              if(info.graph$types[x.n]=="A.B")  mktx<-"AA : BB (1:1)"
            else stop ("invalid type of marker at marker ", info.graph$names[x.n])
              if(info.graph$types[y.n]=="A.B")  mkty<-"AA : BB (1:1)"
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
##end of rf.graph.table
