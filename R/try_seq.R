#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: try_seq.R                                                     ##
## Contains: try_seq, print.try, draw.try, try_seq_by_seq              ##
##                                                                     ##
## Written by Marcelo Mollinari and Cristiane Taniguti                 ##
## copyright (c) 2009, Marcelo Mollinari                               ##
##                                                                     ##
## First version: 02/27/2009                                           ##
## Description of the function was modified by augusto.garcia@usp.br   ##
## on 2015/07/25
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     #
#######################################################################


##' Try to map a marker into every possible position between markers
##' in a given map
##'
##' For a given linkage map, tries do add an additional unpositioned
##' marker.  This function estimates parameters for all possible maps
##' including the new marker in all possible positions, while keeping
##' the original linkage map unaltered.
##'
##' @aliases try_seq
##'
##' @param input.seq an object of class \code{sequence} with a
##'     predefined order.
##'
##' @param mrk the index of the marker to be tried, according to the
##'     input file.
##'
##' @param tol tolerance for the C routine, i.e., the value used to
##'     evaluate convergence.
##'
##' @param pos defines in which position the new marker \code{mrk}
##'     should be placed for the diagnostic graphic. If \code{NULL}
##'     (default), the marker is placed on the best position i.e. the
##'     one which results LOD = 0.00
##'
##' @param verbose if \code{FALSE} (default), simplified output is
##'     displayed.  if \code{TRUE}, detailed output is displayed.
##'
##' @return An object of class \code{try}, which is a list containing
##'     the following components: \item{ord}{a \code{list} containing
##'     results for every linkage map estimated.  These results
##'     include linkage phases, recombination frequencies and
##'     log-likelihoods.} \item{LOD}{a \code{vector} with LOD-Scores
##'     for each position where the additional marker is placed. This
##'     Score is based on the best combination of linkage phases for
##'     each map.}  \item{try.ord}{a \code{matrix} with the orders of
##'     all linkage maps.}  \item{data.name}{name of the object of
##'     class \code{onemap} with the raw data.} \item{twopt}{name of
##'     the object of class \code{rf_2pts} with the 2-point analyses.}
##'
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##'
##' @seealso \code{\link[onemap]{make_seq}} and
##'     \code{\link[onemap]{compare}}.
##'
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S.,
##'     Yandell, B.  (2008) \emph{qtl: Tools for analyzing QTL
##'     experiments} R package version 1.09-43
##'
##' Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci
##'     with dominant and missing markers in various crosses from two
##'     inbred lines.  \emph{Genetica} 101: 47-58.
##'
##' Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J.,
##'     Lincoln, S. E. and Newburg, L. (1987) MAPMAKER: An interactive
##'     computer package for constructing primary genetic linkage maps
##'     of experimental and natural populations. \emph{Genomics} 1:
##'     174-181.
##'
##' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
##'     A. A. F.  (2009) Evaluation of algorithms used to order
##'     markers on genetic maps.  \emph{Heredity} 103: 494-502
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a)
##'     Simultaneous maximum likelihood estimation of linkage and
##'     linkage phases in outcrossing species.  \emph{Theoretical
##'     Population Biology} 61: 349-363.
##'
##' Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage
##'     mapping of sex-specific differences. \emph{Genetical Research}
##'     79: 85-96
##'
##' @keywords utilities
##' @examples
##'
##' \donttest{
##'   #outcrossing example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   markers <- make_seq(twopt,c(2,3,12,14))
##'   markers.comp <- compare(markers)
##'   base.map <- make_seq(markers.comp,1)
##'
##'   extend.map <- try_seq(base.map,30)
##'   extend.map
##'   print(extend.map,5) # best position
##'   print(extend.map,4) # second best position
##'
##' }
##'@export
try_seq<-function(input.seq,mrk,
                  tol=10E-2,
                  pos= NULL,
                  verbose=FALSE)
{
  if(inherits(input.seq$data.name, c("outcross","f2")))
    return(try_seq_outcross(input.seq=input.seq,
                            mrk=mrk, tol=tol,
                            pos=pos, verbose=verbose))
  else
    return(try_seq_inbred_bc(input.seq=input.seq,
                             mrk=mrk, tol=tol,
                             pos=pos, verbose=verbose))
}

## Try to map a marker into every possible position between markers
## in a given map (for crosses derived from f2 inbred lines)
try_seq_inbred_f2 <- function(input.seq,mrk,tol=10E-2,pos= NULL,verbose=TRUE)
{
  # checking for correct objects
  if(!inherits(input.seq,"sequence"))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(input.seq$seq.rf[1] == -1 ||
     is.null(input.seq$seq.like))
    stop("You must run 'compare' or 'map' before the 'try_seq' function")
  if(mrk > input.seq$data.name$n.mar)
    stop(deparse(substitute(mrk))," exceeds the number of markers in object ", input.seq$data.name)
  ## allocate variables
  ord.rf <- matrix(NA, length(input.seq$seq.num)+1, length(input.seq$seq.num))
  ord.like <- rep(NA, length(input.seq$seq.num)+1)
  mark.max<-max(nchar(colnames(input.seq$data.name$geno)))
  num.max<-nchar(ncol(input.seq$data.name$geno))
  ## create first order
  try.ord <- c(mrk,input.seq$seq.num)
  if(verbose) cat("TRY", 1,": ", c(mrk,input.seq$seq.num),"\n")
  flush.console()
  seq.temp<-make_seq(input.seq$twopt, arg=try.ord)
  seq.temp$twopt<-input.seq$twopt
  rf.temp<-get_vec_rf_in(seq.temp, acum=FALSE)
  ## estimate parameters for all possible linkage phases for this order
  final.map<-est_map_hmm_f2(geno=t(input.seq$data.name$geno[,try.ord]),
                            error=input.seq$data.name$error[try.ord + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord)),],
                            rf.vec=rf.temp,
                            verbose=verbose,
                            tol=tol)
  ord.rf[1,] <- final.map$rf
  ord.like[1] <- final.map$loglike
  
  ## positioning between markers of the given sequence
  for(i in 1:(length(input.seq$seq.num)-1))
  {
    ## create intermediate orders
    try.ord <- rbind(try.ord,
                     c(input.seq$seq.num[1:i],
                       mrk,
                       input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
    if(verbose)
      cat("TRY", i+1, ": ", try.ord[i+1,], "\n")
    flush.console()
    seq.temp<-make_seq(input.seq$twopt, arg=try.ord[i+1,])
    seq.temp$twopt<-input.seq$twopt
    rf.temp<-get_vec_rf_in(seq.temp, acum=FALSE)
    ## estimate parameters for all possible linkage phases for this order
    final.map<-est_map_hmm_f2(geno=t(input.seq$data.name$geno[,try.ord[i+1,]]),
                              error=input.seq$data.name$error[try.ord[i+1,] + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord[i+1,])),],
                              rf.vec=rf.temp,
                              verbose=verbose,
                              tol=tol)
    ord.rf[i+1,] <- final.map$rf
    ord.like[i+1] <- final.map$loglike
  }
  ## positioning after the given sequence
  ## create last order
  try.ord <- rbind(try.ord,c(input.seq$seq.num,mrk))
  if(verbose) cat("TRY",length(input.seq$seq.num)+1,": ", c(input.seq$seq.num,mrk) ,"\n")
  flush.console()
  ## estimate parameters for all possible linkage phases for this order
  final.map<-est_map_hmm_f2(geno=t(input.seq$data.name$geno[,try.ord[length(input.seq$seq.num)+1,]]),
                            error=input.seq$data.name$error[try.ord[length(input.seq$seq.num)+1,] + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord[length(input.seq$seq.num)+1,])),],
                            rf.vec=rf.temp,
                            verbose=verbose,
                            tol=tol)
  ord.rf[length(input.seq$seq.num)+1,] <- final.map$rf
  ord.like[length(input.seq$seq.num)+1] <- final.map$loglike
  ## calculate LOD-Scores (best linkage phase combination for each position)
  LOD <- (ord.like-max(ord.like))/log(10)
  ord<-vector("list", nrow(try.ord))
  for(i in 1:nrow(try.ord))
  {
    ord[[i]]<-list(rf=matrix(ord.rf[i,], nrow=1),
                   phase=matrix(rep(NA,ncol(ord.rf)), nrow=1),
                   like=ord.like[i])
    
  }
  w<-structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "try")
  w
}

## Try to map a marker into every possible position between markers
## in a given map (for crosses derived from bc and rils inbred lines)
try_seq_inbred_bc <- function(input.seq,mrk,tol=10E-2,pos= NULL,verbose=TRUE)
{
  # checking for correct objects
  if(!inherits(input.seq,"sequence"))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(input.seq$seq.rf[1] == -1 ||
     is.null(input.seq$seq.like))
    stop("You must run 'compare' or 'map' before the 'try_seq' function")
  if(mrk > input.seq$data.name$n.mar)
    stop(deparse(substitute(mrk))," exceeds the number of markers in object ", input.seq$data.name)
  ## allocate variables
  ord.rf <- matrix(NA, length(input.seq$seq.num)+1, length(input.seq$seq.num))
  ord.like <- rep(NA, length(input.seq$seq.num)+1)
  mark.max<-max(nchar(colnames(input.seq$data.name$geno)))
  num.max<-nchar(ncol(input.seq$data.name$geno))
  ## create first order
  try.ord <- c(mrk,input.seq$seq.num)
  if(verbose) cat("TRY", 1,": ", c(mrk,input.seq$seq.num),"\n")
  flush.console()
  seq.temp<-make_seq(input.seq$twopt, arg=try.ord)
  seq.temp$twopt<-input.seq$twopt
  rf.temp<-get_vec_rf_in(seq.temp, acum=FALSE)
  ## estimate parameters for all possible linkage phases for this order
  final.map<-est_map_hmm_bc(geno=t(input.seq$data.name$geno[,try.ord]),
                            error=input.seq$data.name$error[try.ord + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord)),],
                            rf.vec=rf.temp,
                            verbose=verbose,
                            tol=tol)
  if(inherits(input.seq$data.name, c("riself","risib")))
    final.map$rf<-adjust_rf_ril(final.map$rf,
                                type=class(input.seq$data.name)[2],
                                expand = FALSE)
  ord.rf[1,] <- final.map$rf
  ord.like[1] <- final.map$loglike
  
  ## positioning between markers of the given sequence
  for(i in 1:(length(input.seq$seq.num)-1))
  {
    ## create intermediate orders
    try.ord <- rbind(try.ord,
                     c(input.seq$seq.num[1:i],
                       mrk,
                       input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
    if(verbose) cat("TRY", i+1, ": ", try.ord[i+1,], "\n")
    flush.console()
    seq.temp<-make_seq(input.seq$twopt, arg=try.ord[i+1,])
    seq.temp$twopt<-input.seq$twopt
    rf.temp<-get_vec_rf_in(seq.temp, acum=FALSE)
    ## estimate parameters for all possible linkage phases for this order
    final.map<-est_map_hmm_bc(geno=t(input.seq$data.name$geno[,try.ord[i+1,]]),
                              error=input.seq$data.name$error[try.ord[i+1,] + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord[i+1,])),],
                              rf.vec=rf.temp,
                              verbose=verbose,
                              tol=tol)
    if(inherits(input.seq$data.name, c("riself", "risib")))
      final.map$rf<-adjust_rf_ril(final.map$rf,
                                  type=class(input.seq$data.name)[2],
                                  expand = FALSE)
    ord.rf[i+1,] <- final.map$rf
    ord.like[i+1] <- final.map$loglike
  }
  ## positioning after the given sequence
  ## create last order
  try.ord <- rbind(try.ord,c(input.seq$seq.num,mrk))
  if(verbose) cat("TRY",length(input.seq$seq.num)+1,": ", c(input.seq$seq.num,mrk) ,"\n")
  flush.console()
  ## estimate parameters for all possible linkage phases for this order
  final.map<-est_map_hmm_bc(geno=t(input.seq$data.name$geno[,try.ord[length(input.seq$seq.num)+1,]]),
                            error=input.seq$data.name$error[try.ord[length(input.seq$seq.num)+1,] + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(try.ord[length(input.seq$seq.num)+1,])),],
                            rf.vec=rf.temp,
                            verbose=verbose,
                            tol=tol)
  if(inherits(input.seq$data.name, c("riself", "risib")))
    final.map$rf<-adjust_rf_ril(final.map$rf,
                                type=class(input.seq$data.name)[2],
                                expand = FALSE)
  ord.rf[length(input.seq$seq.num)+1,] <- final.map$rf
  ord.like[length(input.seq$seq.num)+1] <- final.map$loglike
  ## calculate LOD-Scores (best linkage phase combination for each position)
  LOD <- (ord.like-max(ord.like))/log(10)
  ord<-vector("list", nrow(try.ord))
  for(i in 1:nrow(try.ord))
  {
    ord[[i]]<-list(rf=matrix(ord.rf[i,], nrow=1),
                   phase=matrix(rep(NA,ncol(ord.rf)), nrow=1),
                   like=ord.like[i])
  }
  w<-structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=input.seq$data.name, twopt=input.seq$twopt), class = "try")
  w
}

## Try to map a marker into every possible position between markers
## in a given map (for outcrosses)
try_seq_outcross<- function(input.seq,mrk,
                            tol=10E-2,
                            pos= NULL,
                            verbose=TRUE)
{
  ## checking for correct objects
  if(!inherits(input.seq,"sequence"))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  if(input.seq$seq.phases[1] == -1 ||
     input.seq$seq.rf[1] == -1 ||
     is.null(input.seq$seq.like))
    stop("You must run 'compare' or 'map' before the 'try_seq' function")
  if(mrk > input.seq$data.name$n.mar)
    stop(deparse(substitute(mrk))," exceeds the number of markers in object ", input.seq$data.name)
  
  # allocate variables
  rf.init <- vector("list",length(input.seq$seq.num))
  phase.init <- vector("list",length(input.seq$seq.num))
  ord <- list(list(matrix(NA,16,length(input.seq$seq.num)),
                   matrix(NA,16,length(input.seq$seq.num)),
                   rep(-Inf,16)))
  names(ord[[1]]) <- c("rf","phase","like")
  ord <- rep(ord,length(input.seq$seq.num)+1)
  best.seq <- rep(-Inf,length(input.seq$seq.num)+1)
  store_probs <- rep(as.list(NA), length(input.seq$seq.num)+1)
  store_probs <- lapply(store_probs, function(x) as.list(rep(x, 2)))
  
  ### positioning before the given sequence
  # get two-point information
  list.init <- phases(make_seq(input.seq$twopt,c(mrk,input.seq$seq.num[1]),twopt=input.seq$twopt))
  rf.init[[1]] <- list.init$rf.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) rf.init[[j+1]] <- input.seq$seq.rf[j]
  phase.init[[1]] <- list.init$phase.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) phase.init[[j+1]] <- input.seq$seq.phases[j]
  Ph.Init <- comb_ger(phase.init)
  Rf.Init <- comb_ger(rf.init)
  mark.max<-max(nchar(colnames(input.seq$data.name$geno)))
  num.max<-nchar(ncol(input.seq$data.name$geno))
  
  # create first order
  try.ord <- c(mrk,input.seq$seq.num)
  if(verbose) cat("TRY", 1,": ", c(mrk,input.seq$seq.num),"\n")
  flush.console()
  
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem_amb_ph(M=Ph.Init, w=input.seq, seq.num=c(mrk,input.seq$seq.num))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(is(Ph.Init, "numeric") || is(Ph.Init,"integer")){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est_map_hmm_out(geno=t(input.seq$data.name$geno[,c(mrk,input.seq$seq.num)]),
                                 error=input.seq$data.name$error[c(mrk,input.seq$seq.num) + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(c(mrk,input.seq$seq.num))),],
                                 type=input.seq$data.name$segr.type.num[c(mrk,input.seq$seq.num)],
                                 phase=Ph.Init[j,],
                                 rf.vec=Rf.Init[j,],
                                 verbose=verbose,
                                 tol=tol)
    ord[[1]]$rf[j,] <- final.map$rf
    ord[[1]]$phase[j,] <- Ph.Init[j,]
    ord[[1]]$like[j] <- final.map$loglike
    store_probs[[1]][[j]] <- final.map$probs 
    best.seq[1] <- max(best.seq[1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[1]]$like, decreasing=TRUE)
  ord[[1]]$rf <- ord[[1]]$rf[ord.ind,]
  ord[[1]]$phase <- ord[[1]]$phase[ord.ind,]
  ord[[1]]$like <- ord[[1]]$like[ord.ind]
  store_probs[[1]] <- store_probs[[1]][ord.ind]
  
  ### positioning between markers of the given sequence
  for(i in 1:(length(input.seq$seq.num)-1)) {
    # get two-point information
    list.init <- phases(make_seq(input.seq$twopt,c(input.seq$seq.num[i],mrk,input.seq$seq.num[i+1]),twopt=input.seq$twopt))
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
    Ph.Init <- comb_ger(phase.init)
    Rf.Init <- comb_ger(rf.init)
    
    # create intermediate orders
    try.ord <- rbind(try.ord,c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
    if(verbose) cat("TRY", i+1,": ",c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]) ,"\n")
    flush.console()
    
    if(nrow(Ph.Init)>1){
      ##Removing ambigous phases
      rm.ab<- rem_amb_ph(M=Ph.Init, w=input.seq, seq.num=c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]))
      Ph.Init <- Ph.Init[rm.ab,]
      Rf.Init <- Rf.Init[rm.ab,]
      if(is(Ph.Init, "numeric") || is(Ph.Init,"integer")){
        Ph.Init<-matrix(Ph.Init,nrow=1)
        Rf.Init<-matrix(Rf.Init,nrow=1)
      }
    }
    ## estimate parameters for all possible linkage phases for the current order
    for(j in 1:nrow(Ph.Init)) {
      final.map <- est_map_hmm_out(geno=t(input.seq$data.name$geno[,c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)])]),
                                   error=input.seq$data.name$error[c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]) + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)]))),],
                                   type=input.seq$data.name$segr.type.num[c(input.seq$seq.num[1:i], mrk, input.seq$seq.num[(i+1):length(input.seq$seq.num)])],
                                   phase=Ph.Init[j,],
                                   rf.vec=Rf.Init[j,],
                                   verbose=verbose,
                                   tol=tol)
      ord[[i+1]]$rf[j,] <- final.map$rf
      ord[[i+1]]$phase[j,] <- Ph.Init[j,]
      ord[[i+1]]$like[j] <- final.map$loglike
      best.seq[i+1] <- max(best.seq[i+1],final.map$loglike)
      store_probs[[i+1]][[j]] <- final.map$probs 
    }
    # sort linkage phases by log-likelihood
    ord.ind <- order(ord[[i+1]]$like, decreasing=TRUE)
    ord[[i+1]]$rf <- ord[[i+1]]$rf[ord.ind,]
    ord[[i+1]]$phase <- ord[[i+1]]$phase[ord.ind,]
    ord[[i+1]]$like <- ord[[i+1]]$like[ord.ind]
    store_probs[[i+1]] <- store_probs[[i+1]][ord.ind]
  }
  
  ### positioning after the given sequence
  # get two-point information
  list.init <- phases(make_seq(input.seq$twopt,c(input.seq$seq.num[length(input.seq$seq.num)],mrk),twopt=input.seq$twopt))
  rf.init[[(length(input.seq$seq.num))]] <- list.init$rf.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) rf.init[[j]] <- input.seq$seq.rf[j]
  phase.init[[(length(input.seq$seq.num))]] <- list.init$phase.init[[1]]
  for(j in 1:(length(input.seq$seq.num)-1)) phase.init[[j]] <- input.seq$seq.phases[j]
  Ph.Init <- comb_ger(phase.init)
  Rf.Init <- comb_ger(rf.init)
  
  # create last order
  try.ord <- rbind(try.ord,c(input.seq$seq.num,mrk))
  if(verbose) cat("TRY",length(input.seq$seq.num)+1,": ", c(input.seq$seq.num,mrk) ,"\n")
  flush.console()
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem_amb_ph(M=Ph.Init, w=input.seq, seq.num=c(input.seq$seq.num,mrk))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(is(Ph.Init, "numeric") || is(Ph.Init,"integer")){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est_map_hmm_out(geno=t(input.seq$data.name$geno[,c(input.seq$seq.num,mrk)]),
                                 error=input.seq$data.name$error[c(input.seq$seq.num,mrk) + rep(c(0:(input.seq$data.name$n.ind-1))*input.seq$data.name$n.mar, each=length(c(input.seq$seq.num,mrk))),],
                                 type=input.seq$data.name$segr.type.num[c(input.seq$seq.num,mrk)],
                                 phase=Ph.Init[j,],
                                 rf.vec=Rf.Init[j,],
                                 verbose=verbose,
                                 tol=tol)
    ord[[length(input.seq$seq.num)+1]]$rf[j,] <- final.map$rf
    ord[[length(input.seq$seq.num)+1]]$phase[j,] <- Ph.Init[j,]
    ord[[length(input.seq$seq.num)+1]]$like[j] <- final.map$loglike
    store_probs[[length(input.seq$seq.num)+1]][[j]] <- final.map$probs
    best.seq[length(input.seq$seq.num)+1] <- max(best.seq[length(input.seq$seq.num)+1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[length(input.seq$seq.num)+1]]$like, decreasing=TRUE)
  ord[[length(input.seq$seq.num)+1]]$rf <- ord[[length(input.seq$seq.num)+1]]$rf[ord.ind,]
  ord[[length(input.seq$seq.num)+1]]$phase <- ord[[length(input.seq$seq.num)+1]]$phase[ord.ind,]
  ord[[length(input.seq$seq.num)+1]]$like <- ord[[length(input.seq$seq.num)+1]]$like[ord.ind]
  store_probs[[length(input.seq$seq.num)+1]] <- store_probs[[length(input.seq$seq.num)+1]][ord.ind]
  
  # calculate LOD-Scores (best linkage phase combination for each position)
  LOD <- (best.seq-max(best.seq))/log(10)
  structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=input.seq$data.name, twopt=input.seq$twopt, probs = store_probs), class = "try")
}

##' Print method for object class 'try'
##'
##' @aliases print.try
##'
##' @param x an object of class \code{try}.
##'
##' @param j if \code{NULL} (default), output is a summary of the
##'     results for all possible positions of the additional
##'     marker. Otherwise, an integer makes detailed output to be
##'     printed for the corresponding position. This integer must be
##'     less than or equal to the length of the original sequence plus
##'     1.  @param \dots further arguments, passed to other
##'     methods. Currently ignored.
##' @param ... currently ignored
##' @return \code{NULL}
##' @keywords internal
##' @method print try
##' @export
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

# Defunct draw.try function:
# The diagnostic graphic is made of three figures: i) the top figure
# represents the new genetic map obtained with the insertion of the
# new marker \code{mrk} on position \code{pos}. If \code{pos = NULL}
# (default), the marker is placed on the best position i.e. the one
# which results LOD = 0.00, which is indicated by a red triangle;
# ii) the left bottom figure represents the base map (contained in
# \code{input.seq}) on x-axis and the LOD-Scores of the linkage maps
# obtained with the new marker \code{mrk} tested at the beginning,
# between and at the end of the base map. Actually, it is a graphic
# representation of the \code{LOD} vector (see \code{Value}
# section). The red triangle indicates the best position where the
# new marker \code{mrk} should be placed; iii) the right bottom
# figure is the non-interactive \code{\link[onemap]{rf_graph_table}}
# function for the new genetic map. It plots a matrix of pairwise
# recombination fractions (under the diagonal) and LOD Scores (upper
# the diagonal) using a color scale.

draw.try<-function(base.input, try.input, pos=NULL){
  .Defunct(msg = "Defunct since version 2.0.9")
}

#' Run try_seq considering previous sequence
#' 
#' It uses try_seq function repeatedly trying to positioned each marker 
#' in a vector of markers into a already ordered sequence.
#' Each marker in the vector \code{"markers"} is kept in the sequence 
#' if the difference of LOD and total group size of the models 
#' with and without the marker are below the thresholds \code{"lod.thr"} and \code{"cM.thr"}. 
#' 
#' @param sequence object of class sequence with ordered markers
#' @param markers vector of integers defining the marker numbers to be inserted in the \code{sequence}
#' @param cM.thr number defining the threshold for total map size increase when inserting a single marker
#' @param lod.thr the difference of LODs between model before and after inserting the marker need to have 
#' value higher than the value defined in this argument 
##' @param verbose A logical, if TRUE it output progress status
##' information.
#' 
#' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{onemap} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
##' 2-point analyses.}
#' 
#' @export
try_seq_by_seq <- function(sequence, markers, cM.thr=10, lod.thr=-10, verbose=TRUE){
  
  if(!inherits(sequence, c("sequence"))) stop("Input object must be of class sequence")
  
  seq_now <- sequence
  for(i in 1:length(markers)){
    try_edit <- try_seq(seq_now, markers[i])  
    pos <- which(try_edit$LOD == 0)[1]
    new_map <- make_seq(try_edit, pos)
    size_new <- cumsum(kosambi(new_map$seq.rf))[length(new_map$seq.rf)]
    size_old <- cumsum(kosambi(seq_now$seq.rf))[length(seq_now$seq.rf)]
    lod_new <- new_map$seq.like
    lod_old <- seq_now$seq.like
    diff_size <- size_new - size_old
    diff_lod <- lod_new - lod_old
    if(diff_size < cM.thr & diff_lod > lod.thr){
      seq_now <- new_map
      if(verbose) cat("Marker", markers[i], "was included \n")
    } 
  }
  return(seq_now)
}
