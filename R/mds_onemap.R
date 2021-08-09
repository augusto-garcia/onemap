#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: mds_onemap.R                                                  ##
## Contains: mds_onemap                                                ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
## copyright (c) 2007-9, Cristiane Taniguti                            ##
##                                                                     ##
## First version: 22/11/2019                                           ## 
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################


#' OneMap interface with MDSMap package with possibilitie of multipoint distances estimation
#'
#' For a given sequence of markers, apply mds method described in Preedy and Hackett (2016)
#' using MDSMap package to ordering markers and estimates the genetic distances with OneMap
#' multipoint approach. Also gives MDSMap input file format for directly analysis in this package.
#'
#' @param input.seq an object of class \code{sequence}
#' @param out.file path to the generated MDSMap input file.
#' @param p Integer - the penalty for deviations from the sphere - higher p
#' forces points more closely onto a sphere.
#' @param displaytext Shows markers names in analysis graphic view
#' @param weightfn Character string specifying the values to use for the weight
#' matrix in the MDS 'lod2' or 'lod'.
#' @param mapfn Character string specifying the map function to use on the
#' recombination fractions 'haldane' is default, 'kosambi' or 'none'.
#' @param ispc Logical determining the method to be used to estimate the map. By default 
#' this is TRUE and the method of principal curves will be used. If FALSE then the 
#' constrained MDS method will be used.
#' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
#' if \code{TRUE} one of the markers is removed and mds is performed again.
#' @param size The center size around which an optimum is to be searched
#' @param overlap The desired overlap between batches
#' @param phase_cores The number of parallel processes to use when estimating
#' the phase of a marker. (Should be no more than 4)
#' @param tol tolerance for the C routine, i.e., the value used to evaluate
#' convergence.
#' 
#' @return An object of class \code{sequence}, which is a list containing the
#' following components: \item{seq.num}{a \code{vector} containing the
#' (ordered) indices of markers in the sequence, according to the input file.}
#' \item{seq.phases}{a \code{vector} with the linkage phases between markers
#' in the sequence, in corresponding positions. \code{-1} means that there are
#' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
#' recombination frequencies between markers in the sequence. \code{-1} means
#' that there are no estimated recombination frequencies.}
#' \item{seq.like}{log-likelihood of the corresponding linkage map.}
#' \item{data.name}{name of the object of class \code{onemap} with the raw
#' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
#' 2-point analyses.}
#' 
#' @details For better description about MDS method, see MDSMap package vignette.
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' @seealso  \url{https://CRAN.R-project.org/package=MDSMap}.
#'
#'
#' @references 
#' Preedy, K. F.  and  Hackett, C. A.  (2016). A rapid marker ordering approach for high-density
#' genetic linkage maps in experimental autotetraploid populations using multidimensional
#' scaling. \emph{Theoretical and Applied Genetics} 129: 2117-2132
#'
#' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A. A. F.
#' (2009) Evaluation of algorithms used to order markers on genetics maps.
#' \emph{Heredity} 103: 494-502.
#'
#' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a) Simultaneous maximum
#' likelihood estimation of linkage and linkage phases in outcrossing species.
#' \emph{Theoretical Population Biology} 61: 349-363.
#'
#' Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage mapping of
#' sex-specific differences. \emph{Genetical Research} 79: 85-96
#'
#'@import MDSMap
#'@importFrom utils write.table
#'@importFrom reshape2 melt
#'
#'@export
mds_onemap <- function(input.seq, 
                       out.file= "out.file", 
                       p = NULL, 
                       ispc=TRUE,
                       displaytext=FALSE, 
                       weightfn='lod2', 
                       mapfn='haldane',
                       rm_unlinked=TRUE, 
                       size = NULL, 
                       overlap = NULL,
                       phase_cores = 1, 
                       tol = 1e-05,
                       hmm=TRUE){
  
  ## checking for correct object
  if(!is(input.seq, "sequence"))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  
  n_ind <- input.seq$data.name$n.ind
  if(is(input.seq$data.name,c("outcross", "f2"))){
    mat<- get_mat_rf_out(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1)
    # Include NA in D1D2 markers
    seg_type <- input.seq$data.name$segr.type.num[input.seq$seq.num]
    for(i in 1:length(seg_type))
      for(j in 1:(length(seg_type)-1))
        if((seg_type[i] == 7 & seg_type[j] == 6) | (seg_type[i] == 6 & seg_type[j] == 7)){
          mat[i,j] <- mat[j,i] <- NA
        }
    
  } else {
    mat<-get_mat_rf_in(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1)
  }
  n_mk <- nrow(mat)
  
  mat.rf <- mat.lod <- matrix(rep(NA, n_mk*n_mk), nrow = n_mk)
  colnames(mat.rf) <- colnames(mat.lod) <- rownames(mat.rf) <- rownames(mat.lod) <- colnames(mat)
  mat.lod[lower.tri(mat.lod)] <- mat[lower.tri(mat)]
  mat.rf[upper.tri(mat.rf)] <- mat[upper.tri(mat)]
  
  df <- melt(mat.lod, na.rm = TRUE)
  df.rf <- melt(mat.rf, na.rm = TRUE)
  
  df <- cbind(df.rf, df$value)
  df <- df[with(df, order(Var1, Var2)),]
  
  n_col <- dim(df)[1]
  vector.file <- apply(df, 1, function(x) paste(x, collapse = " "))
  file.out <- c(paste(n_mk, n_col, collapse = " "), vector.file)
  
  write.table(file.out, file = out.file, col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  
  mds_map <- estimate.map(out.file, p = p, ispc = ispc,
                                   weightfn = weightfn, mapfn = mapfn)
  
  plot(mds_map, displaytext = displaytext)
  
  ord_mds <- match(as.character(mds_map$locimap[,2]), colnames(input.seq$data.name$geno)) 
  seq_mds <- make_seq(input.seq$twopt, ord_mds)
  if(hmm){
    if(phase_cores == 1 | is(input.seq$data.name, c("backcross", "riself", "risib"))){
      mds_map <- map(seq_mds, rm_unlinked = rm_unlinked)
    } else{
      if(is.null(size) | is.null(overlap)){
        stop("If you want to parallelize the HMM in multiple cores (phase_cores != 1) 
             you should also define `size` and `overlap` arguments.")
      } else {
        mds_map <- map_overlapping_batches(input.seq = seq_mds,
                                           size = size, overlap = overlap, 
                                           phase_cores = phase_cores, 
                                           tol=tol, rm_unlinked = rm_unlinked)
      }
    }
    
    if(!is.list(mds_map)) {
      new.seq <- make_seq(input.seq$twopt, mds_map)
      mds_map <- mds_onemap(new.seq, out.file= out.file, 
                            p = NULL, ispc=TRUE,
                            displaytext=displaytext, 
                            weightfn=weightfn, 
                            mapfn=mapfn, 
                            rm_unlinked=rm_unlinked,
                            size = size, 
                            overlap = overlap,
                            phase_cores = phase_cores, 
                            tol = tol)
    }
    return(mds_map)
  } else {
    return(seq_mds)
  }
}


#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@import smacof 
#'@keywords internal
calc.maps.sphere<-function(fname,p=100,weightfn='lod2',mapfn='haldane'){
  lodrf<-calc.pair.rf.lod(fname,weightfn)
  confplotno<-1:lodrf$nloci
  r<-lodrf$rf
  lod<-lodrf$lod
  M<-dmap(r,mapfn)
  nloci=length(confplotno)
  smacofsym<-smacofSym(M,ndim=2,weightmat=lod,itmax=100000)
  smacofsphere<-smacofSphere(M,ndim=2,algorithm="dual",weightmat=lod,penalty=p,itmax=1000000,mod=10,verbose=FALSE)
  mapsphere<-map.to.interval(smacofsphere,nloci)
  length<-mapsphere$chromlength[nloci]
  distmap<-outer(mapsphere$maporder,mapsphere$maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(mapsphere$maporder,mapsphere$maporder,Vectorize(function(i,j)lod[i,j]))
  #stressratio=smacofsphere$stress/smacofsym$stress
  locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  sr=smacofsphere$stress/smacofsym$stress
  ssphere=smacofsphere$stress
  ssym=smacofsym$stress
  nnfit<-calc.nnfit(distmap,lodmap,mapsphere$chromlength)
  locimap<-data.frame(confplotno=confplotno[mapsphere$maporder],
                      locus=locikey$locus[mapsphere$maporder],position=mapsphere$chromlength,
                      nnfit=nnfit$pointfits,row.names=1:nloci)
  removedloci<-n
  
  retlist<-list(smacofsym=smacofsym,smacofsphere=smacofsphere,mapsphere=mapsphere,distmap=distmap,
                lodmap=lodmap,locimap=locimap,length=length,removed=n,locikey=locikey,stressratio=sr,
                ssphere=ssphere,ssym=ssym,meannnfit=nnfit$meanfit)
  class(retlist) <- "onemap.spheremap"
  retlist
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@import smacof
#'@import princurve
#'@keywords internal
calc.maps.pc<-function(fname,spar=NULL,ndim=2,weightfn='lod2',mapfn='haldane'){
  lodrf<-calc.pair.rf.lod(fname,weightfn)
  confplotno<-1:lodrf$nloci
  
  r<-lodrf$rf
  lod<-lodrf$lod
  
  M<-dmap(r,mapfn)
  nloci=length(confplotno)
  
  smacofsym<-smacofSym(M,ndim=ndim,weightmat=lod,itmax=100000)
  pc1<-principal_curve(smacofsym$conf,maxit=150,spar=spar,smoother="smooth_spline")
  scale<-sum(smacofsym$delta)/sum(smacofsym$dhat) 
  # Configuration dissim are based on the normalized observed diss - dhat. 
  # True observed dissimilarities are delta
  maporder<-pc1$ord
  estpos<-pc1$lambda[maporder]*scale*100
  # gives the estimated length from the beginning of the line
  rownames<-lodrf$locinames[maporder]
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder, Vectorize(function(i,j)lod[i,j]))
  rownames(distmap)<-rownames;colnames(distmap)<-rownames
  rownames(lodmap)<-rownames;colnames(lodmap)<-rownames
  locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  nnfit<-calc.nnfit(distmap,lodmap,estpos)
  locimap<-data.frame(confplotno=confplotno[maporder],locus=locikey$locus[maporder],position=estpos,nnfit=nnfit$pointfits,row.names=1:nloci)
  removedloci<-n
  
  retlist<-list(smacofsym=smacofsym,pc=pc1,distmap=distmap,lodmap=lodmap,locimap=locimap,length=max(estpos),removed=n,locikey=locikey,meannnfit=nnfit$meanfit)
  if(ndim == 2) {
    class(retlist) <- "onemap.pcmap"
  } else {
    class(retlist) <- "onemap.pcmap3d"
  }
  return(retlist)
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@importFrom utils read.table
#'@importFrom stats relevel
#'@importFrom reshape cast
#'@keywords internal
calc.pair.rf.lod<-function(fname,weightfn='lod',...){
  if(!file.exists(fname)) {
    fname2 <- paste(fname,'.txt',sep="")
    if(file.exists(fname2)) {
      fname <- fname2
    }
  }  
  nloci<-scan(fname,what=integer(),nmax=1)
  d<-read.table(fname,skip=1,header=FALSE)
  names(d)<-c("name1","name2","rfreq","lodscore")
  if(weightfn=='lod2') d$lodscore<-d$lodscore^2
  dd<-d
  
  missing1<-with(d,unique(as.character(name1[!name1%in%name2])))
  missing2<-with(d,as.character(unique(name2[!name2%in%name1])))
  
  if(length(missing1)>1){
    dd$name1<-as.character(dd$name1);dd$name2<-as.character(dd$name2)
    for(i in 2:length(missing1))dd<-rbind(dd,list(missing1[1],missing1[i],0,0))
  }
  if(length(missing2)>1){
    dd$name1<-as.character(dd$name1);dd$name2<-as.character(dd$name2)
    for(i in 2:length(missing2))dd<-rbind(dd,list(missing2[i],missing2[1],0,0))
  }
  dd$name1<-factor(dd$name1,unique(as.character(dd$name1)))
  dd$name1<-relevel(dd$name1,missing1[1])
  dd$name2<-factor(dd$name2,levels=c(as.character(levels(dd$name1)[2:length(levels(dd$name1))]),as.character(missing2[1])))
  
  d<-dd
  b<-matrix(0,ncol=nloci,nrow=nloci)
  temp<-cast(name1~name2,data=d,value="rfreq",add.missing=TRUE,fill=0)
  tt<-as.matrix(temp[,2:(nloci)])
  colnames(tt)<-names(temp)[2:nloci]
  rownames(tt)<-temp$name1
  b[upper.tri(b)]<-tt[upper.tri(tt,diag=TRUE)]
  rfmat<-b+t(b)
  colnames(rfmat)<-c(rownames(tt)[1],colnames(tt))
  rownames(rfmat)<-c(rownames(tt),colnames(tt)[nloci-1])
  rm(temp,tt,b)
  b<-matrix(0,ncol=nloci,nrow=nloci)  
  temp<-cast(name1~name2,data=d,value="lodscore",add.missing=TRUE,fill=0)
  tt<-as.matrix(temp[,2:(nloci)])  
  colnames(tt)<-names(temp)[2:nloci]
  rownames(tt)<-temp$name1
  b[upper.tri(b)]<-tt[upper.tri(tt,diag=TRUE)]
  lodmat<-b+t(b)  
  colnames(lodmat)<-c(rownames(tt)[1],colnames(tt))
  rownames(lodmat)<-c(rownames(tt),colnames(tt)[nloci-1])
  rfmat[rfmat>0.499999]<-0.499999
  rm(temp,tt,b)
  diag(rfmat)<-NA
  diag(lodmat)<-NA
  lodmat[lodmat<0]<-0  
  list(rf=rfmat,lod=lodmat,nloci=nloci,locinames=rownames(rfmat))
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
estimate.map<-function(fname,p=NULL,ispc=TRUE,ndim=2,weightfn='lod2',mapfn='haldane',D1lim=NULL,D2lim=NULL,D3lim=NULL){
  if(ispc==FALSE){
    map<-calc.maps.sphere(fname, p, weightfn=weightfn, mapfn=mapfn)
  } else {
    map<-calc.maps.pc(fname, spar=p, ndim=ndim, weightfn=weightfn, mapfn=mapfn)
  }
  
  write(paste('Stress:',map$smacofsym$stress),"")
  write(paste('Mean Nearest Neighbour Fit:',map$meannnfit),"")
  return(map)
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@importFrom utils read.csv write.table
#'@keywords internal
calc.nnfit.from.file<-function(estmap,fname,mapfn='haldane',header=FALSE){
  estmap<- read.csv(estmap,header=header)
  lodrf<-calc.pair.rf.lod(fname)
  r<-lodrf$rf
  lod<-lodrf$lod
  M<-dmap(r,mapfn)
  lnames<-colnames(M)
  names<-estmap[,1]
  maporder<- sapply(1:length(names),function(i)which(lnames==names[i]))
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder,Vectorize(function(i,j)lod[i,j]))
  nnfit<-calc.nnfit(distmap,lodmap,estmap[,2])
  newmap<-data.frame(name=estmap[,1],position=estmap[,2],nnfit=nnfit$pointfits)
  if(!is.null(fname)) write.table(newmap,file=fname,sep=',')
  nnfit
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
dmap<-function(rf,mapfn="haldane"){
  if (mapfn=="haldane") return(-0.5*log(1-2*rf)) 
  if (mapfn=="kosambi") return(0.25*log((1+2*rf)/(1-2*rf)))
  if (mapfn=="none") return (rf)
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
convert.polar<-function(mdsobject,nloci){
  conf=mdsobject$conf
  l<-dim(conf)[1]
  start<-l+1-nloci
  if(start>1){
    x<-conf[start:l,1]-conf[1,1]
    y<-conf[start:l,2]-conf[1,2]
  } else {
    x<-conf[start:l,1]
    y<-conf[start:l,2]
  }
  yadd<-ifelse(y<=0,2*pi,0)
  xadd<-ifelse(x<=0,pi,yadd)
  theta<-atan(y/x)+xadd
  newtheta<-sort(theta)
  diff=newtheta[2:(length(newtheta))]-newtheta[1:(length(newtheta)-1)]
  maxd<-max(diff)
  rotation<-ifelse(maxd>(pi/3),-min(newtheta[(which(diff==maxd)+1):length(theta)]),0)
  
  rtheta<-(theta+rotation)%%(2*pi)
  radius<-sqrt(x^2+y^2)
  list(theta=(rtheta-min(rtheta))%%(2*pi),radius=radius)  
}


#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@importFrom stats median
#'@keywords internal
map.to.interval<-function(mdsobject,nloci){
  pol<-convert.polar(mdsobject,nloci) #detrend
  lin<-pol$theta
  radmed<-median(pol$radius)
  scale<-sum(mdsobject$delta[lower.tri(mdsobject$delta)])/sum(mdsobject$dhat) 
  # configuration dissim are based on the normalized observed diss - dhat. 
  # True observed dissimilarities are delta
  rlin<-rank(lin,ties.method="random")
  path<-sapply(1:nloci,function(i)return(lin[which(rlin==i)]))
  maporder<-sapply(1:nloci,function(i)return(which(rlin==i)))
  thetalength<-path-path[1]
  chromlength<-scale*radmed*thetalength*100
  locilength<-chromlength[rlin]
  list(chromlength=chromlength,order=rlin,locilength=locilength,maporder=maporder)
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
get.nearest.informative<-function(loci,lodmap){
  #split matrix by loci
  neighbours<-NULL
  
  if(loci>1) {
    locileft<-lodmap[loci,(loci-1):1]
    if(length(which(locileft!=0))>0)    neighbours<-loci-min(which(locileft!=0))
  }
  if(loci<dim(lodmap)[2]){
    lociright<-lodmap[loci,(loci+1):dim(lodmap)[2]]
    if(length(which(lociright!=0))>0)  neighbours<-c(neighbours,loci+min(which(lociright!=0)))
  }
  neighbours
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
calc.nnfit.loci<-function(loci,distmap,lodmap,estmap){
  nns<-get.nearest.informative(loci,lodmap)
  obs<-distmap[loci,nns]
  est<-estmap[loci]-estmap[nns]
  nn.fit<-sum(abs(obs-est))
  nn.fit
}


#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'@keywords internal
calc.nnfit<-function(distmap,lodmap,estmap){
  pointfits<-unlist(lapply(1:dim(distmap)[2],calc.nnfit.loci,distmap=distmap,lodmap=lodmap,estmap=estmap))
  fit<-sum(pointfits)
  list(fit=fit,pointfits=pointfits,meanfit=mean(pointfits))
}


#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@import graphics
#'@export
plot.onemap.pcmap <- function (x,D1lim=NULL,D2lim=NULL,displaytext=TRUE,...){
  par(mfrow=c(1,2))
  with(x,{
    if (displaytext==TRUE) {
      labels=locikey$locus
    } else {
      labels=locikey$confplotno
    }
    plot(smacofsym$conf,type="n",main='MDS with principal curve',xlim=D1lim,ylim=D2lim,xlab='Dim 1',ylab='Dim 2')
    text(smacofsym$conf,labels=labels,cex=0.8)
    lines(pc)
    if (displaytext==TRUE)  {
      labels1=locimap$locus
    } else  {
      labels1=locimap$confplotno
    }
    plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
    text(locimap$position,locimap$nnfit,labels1)
  })
}


#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@import graphics
#'@export
plot.onemap.pcmap3d <- function (x,D1lim=NULL,D2lim=NULL,D3lim=NULL,displaytext=TRUE,...) {
  par(mfrow=c(2,2))
  with(x,{
    if (displaytext==TRUE) {
      labels=locikey$locus 
    } else {
      labels=locikey$confplotno
    }
    par(mfrow=c(2,2))
    plot(smacofsym$conf[,'D1'],smacofsym$conf[,'D2'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 2',xlim=D1lim,ylim=D2lim)
    text(smacofsym$conf[,'D1'],smacofsym$conf[,'D2'],labels=labels,cex=0.8)
    lines(pc$s[,'D1'][pc$ord],pc$s[,'D2'][pc$ord])
    plot(smacofsym$conf[,'D1'],smacofsym$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 3',xlim=D1lim,ylim=D3lim)
    text(smacofsym$conf[,'D1'],smacofsym$conf[,'D3'],labels=labels,cex=0.8)
    lines(pc$s[,'D1'][pc$ord],pc$s[,'D3'][pc$ord])
    plot(smacofsym$conf[,'D2'],smacofsym$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 2',ylab='Dimension 3',xlim=D2lim,ylim=D3lim)
    text(smacofsym$conf[,'D2'],smacofsym$conf[,'D3'],labels=labels,cex=0.8)
    lines(pc$s[,'D2'][pc$ord],pc$s[,'D3'][pc$ord])
    if (displaytext==TRUE) {
      labels1=locimap$locus
    } else {
      labels1=locimap$confplotno
    }
    plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
    text(locimap$position,locimap$nnfit,labels1)
  })
}

#'@author Katharine F. Preedy, \email{katharine.preedy@bioss.ac.uk}
#'
#'@import graphics
#'@export
plot.onemap.spheremap <- function (x,displaytext=TRUE,...) {
  
  par(mfrow=c(2,2))
  with(x,{
    if (displaytext==TRUE) {
      labels=locikey$locus 
    } else {
      labels=locikey$confplotno
    }
    plot(c(0,1),c(0,1),type='n',axes=F,xlab="",ylab="")
    text(0.5,0.7,paste('Sym Stress =',round(ssym,digits=4)))
    text(0.5,0.55,paste('Sphere Stress/Sym Stress =',round(stressratio,digits=4)))
    text(0.5,0.4,paste('Sphere Stress =',round(ssphere,digits=4)))
    
    plot(smacofsym,plot.type="confplot",type="n",main='Unconstrained',label.conf=list(label=FALSE,pos=1,col=1))
    text(smacofsym$conf,labels=labels,cex=0.8)
    xlower=min(smacofsym$conf[,1],smacofsphere$conf[,1])
    xupper=max(smacofsym$conf[,1],smacofsphere$conf[,1])
    ylower=min(smacofsym$conf[,2],smacofsphere$conf[,2])
    yupper=max(smacofsym$conf[,2],smacofsphere$conf[,2])
    plot(smacofsym,plot.type="confplot",type="n",main='Unconstrained + Spherical',label.conf=list(label=FALSE,pos=1,col=1),xlim=c(xlower,xupper),ylim=c(ylower,yupper))
    text(smacofsym$conf,labels=labels,cex=0.8)
    l=dim(smacofsphere$conf)[1]-1
    text(smacofsphere$conf[1:l+1,],labels=labels,cex=0.8,col="red") 
    if (displaytext==TRUE)  {
      labels1=locimap$locus 
    } else {
      labels1=locimap$confplotno
    }
    plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
    text(locimap$position,locimap$nnfit,labels1)
  })
}
