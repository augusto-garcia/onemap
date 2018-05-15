#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: draw_map.R                                                    #
# Contains: draw_map                                                  #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2010, Marcelo Mollinari                               #
#                                                                     #
# First version: 11/30/2010                                           #
# Last update: 02/19/2011                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################



##' Draw a genetic map
##'
##' Provides a simple draw of a genetic map.
##'
##'@importFrom graphics abline axis lines plot points text
##'
##' @param map.list a map, i.e. an object of class \code{sequence} with a
##' predefined order, linkage phases, recombination fraction and likelihood;
##' also it could be a list of maps.
##' @param horizontal if \code{TRUE}, indicates that the map should be plotted
##' horizontally. Default is \code{FALSE}
##' @param names if \code{TRUE}, displays the names of the markers. Default is
##' \code{FALSE}
##' @param grid if \code{TRUE}, displays a grid in the background. Default is
##' \code{FALSE}
##' @param cex.mrk the magnification to be used for markers.
##' @param cex.grp the magnification to be used for group axis annotation.
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @keywords rqtl
##' @examples
##'
##' \dontrun{
##'  #outcross example
##'   data(example_out)
##'   twopt <- rf_2pts(example_out)
##'   lg<-group(make_seq(twopt, "all"))
##'   maps<-vector("list", lg$n.groups)
##'   for(i in 1:lg$n.groups)
##'      maps[[i]]<- make_seq(order_seq(input.seq= make_seq(lg,i),twopt.alg =
##'    "rcd"), "force")
##'   draw_map(maps, grid=TRUE)
##'   draw_map(maps, grid=TRUE, horizontal=TRUE)
##'
##'   #F2 example
##'   data(onemap_example_f2)
##'   twopt<-rf_2pts(onemap_example_f2)
##'   lg<-group(make_seq(twopt, "all"))
##'   maps<-vector("list", lg$n.groups)
##'   for(i in 1:lg$n.groups)
##'      maps[[i]]<- make_seq(order_seq(input.seq= make_seq(lg,i),twopt.alg =
##'    "rcd"), "force")
##'   draw_map(maps, grid=TRUE)
##'   draw_map(maps, grid=TRUE, horizontal=TRUE)
##'
##' }
##'@export
draw_map<-function(map.list, horizontal=FALSE, names=FALSE, grid=FALSE, cex.mrk=1, cex.grp=.75){
  ## checking for correct object
  if(!any(class(map.list)=="list" | class(map.list)=="sequence")) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequnece'")

  ## if map.list is just a single chormosome, convert it  into a list
  if(class(map.list)=="sequence") map.list<-list(map.list)
  j<-1

  ##converting to data.frame
  out<-data.frame()
  pos<-NULL #to satisfy codetools
  marker<-NULL #to satisfy codetools
  for(i in length(map.list):1){
    if(!any(class(map.list[[i]])=="sequence")) stop("Object ", i , " in map.list is not an object of class 'sequnece'")
    if(is.null(map.list[[i]]$seq.like))  stop("Parameters are not estimated for object ", i, " in map.list")
    map<-cumsum(c(0,get(get(".map.fun", envir=.onemapEnv))(map.list[[i]]$seq.rf)))
    marnames<-colnames(get(map.list[[i]]$data.name, pos=1)$geno)[map.list[[i]]$seq.num]
    out<-rbind(out, data.frame(dist=map, pos=j,marker=marnames))
    j<-j+1
  }
  x<-tapply(out$dist, out$pos, max)
  y<-unlist(unique(out[2]))

  ##Plotting region
  out.fake <-  data.frame(dist=rep(c(0, max(out$dist)),max(y)+2) , pos=c(0:(max(y)+1)))

  if(horizontal==TRUE){
    plot(out.fake, axes=FALSE, col=0, xlab="Distance (cM)", ylab="", main="Genetic Map")
    points(out[1:2], pch="|", cex=cex.mrk, xlab="Distance (cM)", ylab="", main="Genetic Map")
    axis(1, at = seq(from=0, to=10*round(max(x)/10), by=10) , labels=seq(from=0, to=10*round(max(x)/10), by=10) , cex.axis=.75)
    axis(2, y, paste("Group", rev(y)), lwd=0, las=2, cex.axis=cex.grp)
    if(grid==TRUE)
      abline(v=seq(from=0, to=10*round(max(x)/10), by=10), lty=2, lwd=.5, col=2)
    for(i in y){
      if(names==TRUE) text(x=unlist(subset(out, pos==i, dist)), y=i+max(y)/80, labels=unlist(subset(out, pos==i, marker)), srt=90, cex=cex.mrk*.75,  adj = c(0,0.5))
      lines(x=c(0,x[i]), y=c(y[i],y[i]))
    }
  }
  else{
    plot(-out.fake[2:1], axes=FALSE, col=0, ylab="Distance (cM)", xlab="", main="Genetic Map")
    points(-out[2:1], pch= 95, cex=cex.mrk)
    axis(2, cex.axis=.75, at=-seq(from=0, to=10*round(max(x)/10), by=10), labels = seq(from=0, to=10*round(max(x)/10), by=10), las=2)
    axis(1, -y, paste("Group", rev(y)), lwd=0, las=2, cex.axis=cex.grp)
    if(grid==TRUE)
      abline(h=-seq(from=0, to=10*round(max(x)/10), by=10), lty=2, lwd=.5, col=2)
    for(i in y){
      if(names==TRUE) text(x=-i+max(y)/100, y=-unlist(subset(out, pos==i, dist)), labels=unlist(subset(out, pos==i, marker)), cex=cex.mrk*.75,  adj = c(0,0.5))
      lines(y=c(-0.2,-x[i]+0.2), x=c(-y[i],-y[i]))
    }
  }
}
 #end of file
