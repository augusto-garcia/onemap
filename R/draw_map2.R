#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: draw_map2.R                                                   #
# Contains: draw_map2                                                 #
#                                                                     #
# Written by Getulio Caixeta Ferreira                                 #
# copyright (c) 2018, Getulio C Ferreira                              #
#                                                                     #
# First version: 12/12/2018                                           #
# Last update: 12/12/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################



##' Draw a linkage map
##'
##' Provides a simple draw of a linkage map.
##'
##'@importFrom graphics segments
##'@importFrom grDevices bmp dev.off jpeg pdf png postscript tiff
##'
##' @param ... map(s). Object(s) of class \code{sequence} and/or \code{data.frame}. If \code{data.frame}, it must have two columns: column 1: marker id; column 2: position (cM) (\code{numeric}).
##' @param tag name(s) of the marker(s) to highlight. If "all", all markers will be highlighted. Default is \code{NULL}.
##' @param id logical. If \code{TRUE} (default), shows name(s) of tagged marker(s).
##' @param pos logical. If \code{TRUE} (default), shows position(s) of tagged marker(s).
##' @param cex.label the magnification used for label(s) of tagged marker(s). If \code{NULL} (default), the cex will be automatically calculated to avoid overlapping.
##' @param main an overall title for the plot. Default is \code{NULL}.
##' @param group.names name(s) to identfy the group(s). If \code{NULL} (default), the name(s) of the sequence(s) will be used.
##' @param centered logical. If \code{TRUE}, the group(s) will be aligned in the center. If \code{FALSE} (default), the group(s) will be aligned at the top.
##' @param y.axis logical. If \code{TRUE} (default), shows y axis. If centered = \code{TRUE}, the y axis will always be hidden.
##' @param space numerical. Spacing between groups. If \code{NULL} (default), the spacing will be automatically calculated to avoid overlapping.
##' @param col.group the color used for group(s).
##' @param col.mark the color used for marker(s).
##' @param col.tag the color used for highlighted marker(s) and its/theirs label(s).
##' @param output the name of the output file. The file format can be specified by adding its extension. Available formats: 'bmp', 'jpeg', 'png', 'tiff', 'pdf' and 'eps' (default).
##' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
##' @keywords rqtl
##' @examples
##'
##' \dontrun{
##'  data("onemap_example_out")
##'  twopt <- rf_2pts(onemap_example_out)
##'  lg<-group(make_seq(twopt, "all"))
##'  seq1<-make_seq(order_seq(input.seq= make_seq(lg,1),twopt.alg = "rcd"), "force")
##'  seq2<-make_seq(order_seq(input.seq= make_seq(lg,2),twopt.alg = "rcd"), "force")
##'  seq3<-make_seq(order_seq(input.seq= make_seq(lg,3),twopt.alg = "rcd"), "force")
##'  draw_map2(seq1,seq2,seq3,tag = c("M1","M2","M3","M4","M5"))
##'  
##'  data("onemap_example_f2")
##'  twopt <- rf_2pts(onemap_example_f2)
##'  lg<-group(make_seq(twopt, "all"))
##'  seq<-list(
##'  make_seq(order_seq(input.seq= make_seq(lg,1),twopt.alg = "rcd"), "force"),
##'  make_seq(order_seq(input.seq= make_seq(lg,2),twopt.alg = "rcd"), "force"),
##'  make_seq(order_seq(input.seq= make_seq(lg,3),twopt.alg = "rcd"), "force")
##'  )
##'  draw_map2(seq,tag = "all",group.names = c("Chr 1","Chr 2","Chr 3"),main="Linkage Map")
##'
##' }
##'@export
draw_map2<-function(...,tag=NULL,id=TRUE,pos =TRUE,cex.label=NULL,main=NULL,group.names=NULL,centered=F,y.axis=TRUE,space=NULL,col.group=NULL,col.mark=NULL,col.tag=NULL,output=NULL){
  #check input
  input<-list(...)
  if(length(input)==0) stop("argument '...' missing, with no default")
  map.data<-list()
  for(i in seq(input)) map.data<-c(map.data,if(is(input[[i]], "list")) input[i] else do.call(c,input[i]))
  if(!all(sapply(map.data, function(x) is(x, c("sequence","data.frame"))))) stop(paste("\n'",names(map.data)[!sapply(map.data, function(x) is(x, c("sequence","data.frame")))],"' is not an object of class 'sequence' or 'data.frame",sep=""))

  #sequence to data.frame
  for(i in seq_along(map.data)){
    if(is(map.data[[i]], "sequence")){
      if(is.character(map.data[[i]]$data.name)){
        if (!map.data[[i]]$data.name %in% ls(.GlobalEnv)) stop(paste("Object data missing:", map.data[[i]]$data.name))
        map.data[[i]]$data.name <- get(map.data[[i]]$data.name, envir = .GlobalEnv)
      }
      map.data[[i]] <- data.frame(marker=colnames(map.data[[i]]$data.name$geno)[map.data[[i]]$seq.num], pos=c(0,cumsum(kosambi(map.data[[i]]$seq.rf))), chr=i)
    } else {map.data[[i]] <- cbind(map.data[[i]], chr=i)}
  }

  #check data.frame
  for(i in seq_along(map.data)){
    if(ncol(map.data[[i]]) != 3) stop(paste("\n", names(map.data)[i],"has incorrect number of columns."))
    if(any(duplicated(map.data[[i]][,1]))) stop(paste("\n", names(map.data)[i],"has duplicated marker names."))
    if(!is.numeric(map.data[[i]][,2])) stop(paste("\n", names(map.data)[i],", position column (2) must be numeric"))
    colnames(map.data[[i]]) <- c("marker", "pos", "chr")
    map.data[[i]]$marker <- as.character(map.data[[i]]$marker)
    map.data[[i]] <- map.data[[i]][order(map.data[[i]]$pos),]
  }

  nchr<-length(map.data)
  max.pos<-ceiling(max(do.call("rbind",map.data)$pos)/10)*10
  for(i in seq_along(map.data))map.data[[i]]$coord<--map.data[[i]]$pos*100/max.pos
  
  if(is.null(main)) main<-""
  if(is.null(group.names)) group.names <- seq(map.data)
  if(is.null(tag)) pos <- id <- F
  if("all"%in%tag) tag<-unique(unlist(lapply(map.data,function(x) x$marker),use.names = F))
  if(is.null(space)) space<-1+pos+id
  if(centered==T){
    for(i in seq_along(map.data)){map.data[[i]]$coord<-map.data[[i]]$coord-(100+min(map.data[[i]]$coord))/2}
    y.axis<-F
  }
  if(y.axis==T){
    yl<-"Distance (cM)"
    mleft<-5
  } else {
    yl<-""
    mleft<-1
  }
  if(is.null(col.group)) col.group<-"grey85"
  if(is.null(col.mark)) col.mark<-"#cc662f"
  if(is.null(col.tag)) col.tag<-"#003350"
  if(is.null(output)) output<-"map"
  
  
  # Split output
  if(strsplit(output,"")[[1]][length(strsplit(output,"")[[1]])]=="/") output<-paste(output,"map",sep = "")
  output<-strsplit(output,"/",T)[[1]]
  if(length(output)>1){
    output.dir<-paste(output[-length(output)],collapse = "/")
    output<-output[length(output)]
  } else output.dir<-"."
  if(!dir.exists(output.dir)) stop("\nInvalid directory")
  output<-strsplit(output,".",T)[[1]]
  if(output[length(output)]%in%c("bmp","jpeg","png","tiff","pdf","eps")){
    output.ext<-output[length(output)]
    output<-paste(output[1:(length(output)-1)],collapse = ".")
  } else{
    output.ext<-"eps"
    output<-paste(output[1:length(output)],collapse = ".")
  }
  n<-0
  if(paste(output,output.ext,sep = ".")%in%list.files(output.dir)){
    repeat{
      n<-n+1
      output2<-paste(output,"(",n,")",sep = "")
      if(paste(output2,output.ext,sep = ".")%in%list.files(output.dir)==F){
        output<-output2
        break
      }
    }
  }
  
  # Preparing plot area
  if(output.ext=="bmp"){
    bmp(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15,width = nchr*(1+space)+(mleft+1)/2,res = 300,units = "cm")
  } else  if(output.ext=="jpeg"){
    jpeg(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15,width = nchr*(1+space)+(mleft+1)/2,res = 300,units = "cm")
  } else if(output.ext=="png"){
    png(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15,width = nchr*(1+space)+(mleft+1)/2,res = 300,units = "cm")
  } else if(output.ext=="tiff"){
    tiff(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15,width = nchr*(1+space)+(mleft+1)/2,res = 300,units = "cm")
  } else if(output.ext=="pdf"){
    pdf(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15/2.54,width = (nchr*(1+space)+(mleft+1)/2)/2.54)
  } else if(output.ext=="eps"){
    postscript(paste(output.dir,"/",output,".",output.ext,sep=""),height = 15/2.54,width = (nchr*(1+space)+(mleft+1)/2)/2.54,paper = "special",horizontal = F,onefile = F)
  }
  
  par(mar=c(1,mleft,2,1))
  plot(NA,NA,xlim=c(0,nchr*(1+space)),ylim = c(-100,5),type = "n",axes = F,ylab = yl,main=main)
  text(seq((1+space)/2,(nchr-.5)*(1+space),1+space),5,group.names)
  if(y.axis==T){
    axis(2,at=seq(0,-100, -1000/max.pos),labels = seq(0,max.pos, 10), las = 2)
  }
  #Drawing Group(s)
  segments(seq((1+space)/2,(nchr-.5)*(1+space),1+space),tapply(do.call("rbind",map.data)$coord,do.call("rbind",map.data)$chr,max),
           seq((1+space)/2,(nchr-.5)*(1+space),1+space),tapply(do.call("rbind",map.data)$coord,do.call("rbind",map.data)$chr,min),
           lwd = 20,col = col.group)
  segments((do.call("rbind",map.data)$chr-.5)*(1+space)-.25,do.call("rbind",map.data)$coord,
           (do.call("rbind",map.data)$chr-.5)*(1+space)+.25,do.call("rbind",map.data)$coord,
           lwd = 2,col = col.mark)
  
  #Labeling
  tag.data<-do.call("rbind",map.data)[do.call("rbind",map.data)$marker%in%tag,]
  tag.data$coord2<-tag.data$coord
  tag.data<-split(tag.data,tag.data$chr)
  if(length(tag.data)>0){
    if(is.null(cex.label)){
      if(max(table(do.call("rbind",tag.data)$chr))<=30){cex.label<-1
      } else{cex.label=30/max(table(do.call("rbind",tag.data)$chr))}
    }
    if(length(tag.data)>0){
      dif.min<--floor(1000/3*cex.label)/100
      for(i in seq_along(tag.data)){
        if(nrow(tag.data[[i]])>30/cex.label){tag.data[[i]]$coord2<-seq(0,-100,-100/(nrow(tag.data[[i]])-1))
        } else{
          repeat{
            dif<-round(diff(tag.data[[i]]$coord2),2)
            over.dif.min<-c(dif>dif.min,F)
            if(T%in%over.dif.min){
              tag.data[[i]]$coord2[c(over.dif.min)]<-tag.data[[i]]$coord2[c(over.dif.min)]+(-dif.min+dif[over.dif.min])/2
              tag.data[[i]]$coord2[c(F,over.dif.min[-length(over.dif.min)])]<-tag.data[[i]]$coord2[c(F,over.dif.min[-length(over.dif.min)])]-(-dif.min+dif[over.dif.min])/2
            } else break
          }
          repeat{
            if(T%in%(tag.data[[i]]$coord2>0)){
              dif<-round(diff(tag.data[[i]]$coord2),2)
              over.dif.min<-(dif<dif.min)
              next.spot<-(1:length(over.dif.min))[over.dif.min][1]
              if(length(next.spot)==0){tag.data[[i]]$coord2<-tag.data[[i]]$coord2-tag.data[[i]]$coord2[1]
              } else if(is.na(next.spot)){tag.data[[i]]$coord2<-tag.data[[i]]$coord2-tag.data[[i]]$coord2[1]
              } else tag.data[[i]]$coord2[1:next.spot]<-tag.data[[i]]$coord2[1:next.spot]-min(c(tag.data[[i]]$coord2[1],dif.min-dif[next.spot]))
            }
            if(T%in%(tag.data[[i]]$coord2>0)==F) break
          }
          repeat{
            if(T%in%(tag.data[[i]]$coord2<(-100))){
              dif<-round(diff(tag.data[[i]]$coord2),2)
              over.dif.min<-(dif<dif.min)
              next.spot<-(1:length(over.dif.min))[over.dif.min][length((1:length(over.dif.min))[over.dif.min])]+1
              if(length(next.spot)==0){tag.data[[i]]$coord2<-tag.data[[i]]$coord2-100-tag.data[[i]]$coord2[length(tag.data[[i]]$coord2)]
              } else if(is.na(next.spot)){tag.data[[i]]$coord2<-tag.data[[i]]$coord2-100-tag.data[[i]]$coord2[length(tag.data[[i]]$coord2)]
              } else tag.data[[i]]$coord2[length(tag.data[[i]]$coord2):next.spot]<-tag.data[[i]]$coord2[length(tag.data[[i]]$coord2):next.spot]+min(c(-100-tag.data[[i]]$coord2[length(tag.data[[i]]$coord2)],dif.min-dif[next.spot-1]))
            }
            if(T%in%(tag.data[[i]]$coord2<(-100))==F) break
          }
        }
      }
    }
    tag.data<-do.call("rbind",tag.data)
    segments((tag.data$chr-.5)*(1+space)-.25,tag.data$coord,
             (tag.data$chr-.5)*(1+space)+.25,tag.data$coord,lwd = 2,col = col.tag)
    if(id==T){
      segments((tag.data$chr-.5)*(1+space)+.275,tag.data$coord,
               (tag.data$chr-.5)*(1+space)+.5,tag.data$coord2,lwd = .5,col = col.tag)
      text((tag.data$chr-.5)*(1+space)+.4,tag.data$coord2,
           pos = 4,tag.data$marker,col = col.tag,cex = cex.label)
    }
    if(pos==T){
      segments((tag.data$chr-.5)*(1+space)-.275,tag.data$coord,
               (tag.data$chr-.5)*(1+space)-.5,tag.data$coord2,lwd = .5,col = col.tag)
      text((tag.data$chr-.5)*(1+space)-.4,tag.data$coord2,
           pos = 2,round(tag.data$pos,1),col = col.tag,cex = cex.label)
    }
  }
  dev.off()
  cat("Completed\nOutput file: ")
  if(output.dir==".") cat(getwd(),"/",output,".",output.ext,"\n",sep = "") else cat(output.dir,"/",output,".",output.ext,"\n",sep = "")
}
##end of file