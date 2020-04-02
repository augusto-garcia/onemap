#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: order_seq.R                                                   #
# Contains: order_seq, print.order, draw_order                        #
#                                                                     #
# Written by Gabriel R A Margarido & Marcelo Mollinari with minor     #
# changes by Cristiane Taniguti
# copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 08/09/2017                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

## This function automates linkage map construction in two steps:
## first, it applies the 'compare' algorithm to a subset of markers;
## second, it adds markers sequentially with the 'try' function


##' Search for the best order of markers combining compare and try_seq
##' functions
##'
##' For a given sequence of markers, this function first uses the
##' \code{compare} function to create a framework for a subset of informative
##' markers. Then, it tries to map remaining ones using the \code{try_seq}
##' function.
##'
##' For outcrossing populations, the initial subset and the order in which
##' remaining markers will be used in the \code{try_seq} step is given by the
##' degree of informativeness of markers (i.e markers of type A, B, C and D, in
##' this order).
##'
##' For backcrosses, F2s or RILs, two methods can be used for
##' choosing the initial subset: i) \code{"sample"} randomly chooses a number
##' of markers, indicated by \code{n.init}, and calculates the multipoint
##' log-likelihood of the \eqn{\frac{n.init!}{2}}{n.init!/2} possible orders.
##' If the LOD Score of the second best order is greater than
##' \code{subset.THRES}, than it takes the best order to proceed with the
##' \code{try_seq} step. If not, the procedure is repeated. The maximum number
##' of times to repeat this procedure is given by the \code{subset.n.try}
##' argument. ii) \code{"twopt"} uses a two-point based algorithm, given by the
##' option \code{"twopt.alg"}, to construct a two-point based map. The options
##' are \code{"rec"} for RECORD algorithm, \code{"rcd"} for Rapid Chain
##' Delineation, \code{"ser"} for Seriation and \code{"ug"} for Unidirectional
##' Growth. Then, equally spaced markers are taken from this map. The
##' \code{"compare"} step will then be applied on this subset of markers.
##'
##' In both cases, the order in which the other markers will be used in the
##' \code{try_seq} step is given by marker types (i.e. co-dominant before
##' dominant) and by the missing information on each marker.
##'
##' After running the \code{compare} and \code{try_seq} steps, which result in
##' a "safe" order, markers that could not be mapped are "forced" into the map,
##' resulting in a map with all markers positioned.
##'
##'@importFrom graphics abline axis layout lines par plot points text title
##'
##' @param input.seq an object of class \code{sequence}.
##' @param n.init the number of markers to be used in the \code{compare} step
##' (defaults to 5).
##' @param subset.search a character string indicating which method should be
##' used to search for a subset of informative markers for the
##' \code{\link{compare}} step. It is used for backcross, \eqn{F_2}{F_2} or RIL
##' populations, but not for outcrosses. See the \code{Details} section.
##' @param subset.n.try integer. The number of times to repeat the subset
##' search procedure. It is only used if \code{subset.search=="sample"}. See
##' the \code{Details} section.
##' @param subset.THRES numerical. The threshold for the subset search
##' procedure. It is only used if \code{subset.search=="sample"}. See the
##' \code{Details} section.
##' @param twopt.alg a character string indicating which two-point algorithm
##' should be used if \code{subset.search=="twopt"}. See the \code{Details}
##' section.
##' @param THRES threshold to be used when positioning markers in the
##' \code{try_seq} step.
##' @param touchdown logical. If \code{FALSE} (default), the \code{try_seq}
##' step is run only once, with the value of \code{THRES}. If \code{TRUE},
##' \code{try_seq} runs with \code{THRES} and then once more, with
##' \code{THRES-1}. The latter calculations take longer, but usually are able
##' to map more markers.
##' @param tol tolerance number for the C routine, i.e., the value used to
##' evaluate convergence of the EM algorithm.
##' @return An object of class \code{order}, which is a list containing the
##' following components: \item{ord}{an object of class \code{sequence}
##' containing the "safe" order.} \item{mrk.unpos}{a \code{vector} with
##' unpositioned markers (if they exist).} \item{LOD.unpos}{a \code{matrix}
##' with LOD-Scores for unmapped markers, if any, for each position in the
##' "safe" order.} \item{THRES}{the same as the input value, just for
##' printing.} \item{ord.all}{an object of class \code{sequence} containing the
##' "forced" order, i.e., the best order with all markers.}
##' \item{data.name}{name of the object of class \code{onemap} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
##' 2-point analyses.}
##' @author Gabriel R A Margarido, \email{gramarga@@usp.br} and Marcelo
##' Mollinari, \email{mmollina@@gmail.com}
##' @seealso \code{\link[onemap]{make_seq}}, \code{\link[onemap]{compare}} and
##' \code{\link[onemap]{try_seq}}.
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci with
##' dominant and missing markers in various crosses from two inbred lines.
##' \emph{Genetica} 101: 47-58.
##'
##' Lander, E. S. and Green, P. (1987). Construction of multilocus genetic
##' linkage maps in humans. \emph{Proc. Natl. Acad. Sci. USA} 84: 2363-2367.
##'
##' Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J., Lincoln,
##' S. E. and Newburg, L. (1987) MAPMAKER: An interactive computer package for
##' constructing primary genetic linkage maps of experimental and natural
##' populations. \emph{Genomics} 1: 174-181.
##'
##' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A. A. F.
##' (2009) Evaluation of algorithms used to order markers on genetics maps.
##' \emph{Heredity} 103: 494-502.
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##'
##' Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage mapping of
##' sex-specific differences. \emph{Genetical Research} 79: 85-96
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'   #outcross example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG2 <- make_seq(groups,2)
##'   LG2.ord <- order_seq(LG2,touchdown=TRUE)
##'   LG2.ord
##'   make_seq(LG2.ord) # get safe sequence
##'   make_seq(LG2.ord,"force") # get forced sequence
##'
##'   #F2 example
##'   data(onemap_example_f2)
##'   twopt <- rf_2pts(onemap_example_f2)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG3 <- make_seq(groups,3)
##'   LG3.ord <- order_seq(LG3, subset.search = "twopt", twopt.alg = "rcd", touchdown=TRUE)
##'   LG3.ord
##'   make_seq(LG3.ord) # get safe sequence
##'   ord.1<-make_seq(LG3.ord,"force") # get forced sequence
##'
##'   LG3.ord.s <- order_seq(LG3, subset.search = "sample", touchdown=TRUE)
##'   LG3.ord.s
##'   make_seq(LG3.ord) # get safe sequence
##'   ord.2<-make_seq(LG3.ord,"force") # get forced sequence
##'
##'   rbind(ord.1$seq.num, ord.2$seq.num) # probably, the same order for
##'   this dataset
##' }
##'@export
order_seq <- function(input.seq, n.init=5, subset.search=c("twopt", "sample"),
                      subset.n.try=30, subset.THRES=3, twopt.alg= c("rec", "rcd", "ser", "ug"),
                      THRES=3, touchdown=FALSE, tol=10E-2) {
  ## checking for correct objects
  if(!is(input.seq,"sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(n.init < 2) stop("'n.init' must be greater than or equal to 2")
  if(!is.logical(touchdown)) stop("'touchdown' must be logical")
  if(!touchdown && THRES <= 10E-10) stop("Threshold must be greater than 0 if 'touchdown' is FALSE")
  if(touchdown && THRES <= (1 + 10E-10)) stop("Threshold must be greater than 1 if 'touchdown' is TRUE")
  
  if(length(input.seq$seq.num) <= n.init) {
    ## in this case, only the 'compare' function is used
    cat("   Length of sequence ",deparse(substitute(input.seq))," is less than n.init \n   Returning the best order using compare function:\n")
    ifelse(length(input.seq$seq.num) == 2, seq.ord <- map(input.seq,tol=10E-5), seq.ord <- make_seq(compare(input.seq=input.seq,tol=10E-5),1))
    seq.ord<-map(seq.ord, tol=10E-5)
    structure(list(ord=seq.ord, mrk.unpos=NULL, LOD.unpos=NULL, THRES=THRES,
                   ord.all=seq.ord, data.name=input.seq$data.name, probs = seq.ord$probs, twopt=input.seq$twopt), class = "order")
  }
  else
  {
    ## here, the complete algorithm will be applied
    cross.type <- class(get(input.seq$data.name, pos=1))[2]
    if(cross.type == "f2") FLAG <- "f2"
    else if (cross.type == "backcross" ||
             cross.type == "riself" ||
             cross.type == "risib") FLAG <- "bc"
    else if (cross.type == "outcross") FLAG <- "outcross"
    else stop("Invalid cross type\n")
    ## select the order in which markers will be added
    if(FLAG == "bc"){
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
        else stop("Invalid cross type\n")
      }
      else if(subset.search=="sample"){
        cat("\nCross type: ", cross.type, "\nChoosing initial subset using the 'sample' approach\n")
        LOD.test <- i <- 0
        while(abs(LOD.test) < abs(subset.THRES) && i < subset.n.try){
          smp.seq <- make_seq(get(input.seq$twopt), sample(input.seq$seq.num, size=n.init), twopt=input.seq$twopt)
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
          else stop("Invalid cross type\n")
        }
        else stop("Cannot find any subset using 'subset.n.try'=", subset.n.try, " and 'subset.THRES'= ", subset.THRES,"\n")
      }
      ## else stop("Invalid subset search\n")
    }
    else if(FLAG == "outcross" || FLAG == "f2") {
      cat(paste("\nCross type:", FLAG, "\nUsing segregation types of the markers to choose initial subset\n"))
      segregation.types <- get(input.seq$data.name, pos=1)$segr.type.num[input.seq$seq.num]
      if(sum(segregation.types == 7) > sum(segregation.types == 6)) segregation.types[segregation.types == 6] <- 8 ## if there are more markers of type D2 than D1, try to map those first
      seq.work <- order(segregation.types)
      seq.init <- input.seq$seq.num[seq.work[1:n.init]]
    }
    else stop("Invalid cross type")
    ##apply the 'compare' step to the subset of initial markers
    seq.ord <- compare(input.seq=make_seq(get(input.seq$twopt), seq.init, twopt=input.seq$twopt), n.best=50)
    
    ## 'try' to map remaining markers
    input.seq2 <- make_seq(seq.ord,1)
    cat ("\n\nRunning try algorithm\n")
    for (i in (n.init+1):length(input.seq$seq.num)){
      seq.ord <- try_seq(input.seq2,input.seq$seq.num[seq.work[i]],tol=tol)
      if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < -THRES))
        input.seq2 <- make_seq(seq.ord,which.max(seq.ord$LOD))
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
        seq.ord <- try_seq(input.seq2,i,tol=tol)
        if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < (-THRES+1)))
          input.seq2 <- make_seq(seq.ord,which.max(seq.ord$LOD))
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
        LOD.unpos[j,] <- try_seq(input.seq=input.seq2,mrk=i,tol=tol)$LOD
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
        seq.ord <- try_seq(input.seq3,i,tol)
        input.seq3 <- make_seq(seq.ord,which(seq.ord$LOD==0)[sample(sum(seq.ord$LOD==0))[1]])
      }
    }
    cat("\nEstimating final genetic map using tol = 10E-5.\n\n")
    input.seq2<-map(input.seq2, tol=10E-5)
    input.seq3<-map(input.seq3, tol=10E-5)
    structure(list(ord=input.seq2, mrk.unpos=mrk.unpos, LOD.unpos=LOD.unpos, THRES=THRES,
                   ord.all=input.seq3, data.name=input.seq$data.name, probs2 = input.seq2$probs, 
                   probs3 = input.seq3$probs, twopt=input.seq$twopt), class = "order")
  }
}

##'@export
##'@method print order
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

draw_order<-function(map.input){
  .Defunct(msg = "Defunct since version 2.0.9")
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
  rf_graph_table(map.input, inter=FALSE, main="")
  title(main = "LOD (above diag.) and Recombination Fraction Matrix", cex.main=.9, line=15.4)
}
## end of file
