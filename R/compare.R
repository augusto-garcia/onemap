#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: compare.R                                                     ##
## Contains: compare, print.compare                                    ##
##                                                                     ##
## Written by Gabriel R A Margarido & Marcelo Mollinari                ##
## copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       ##
##                                                                     ##
## First version: 02/27/2009                                           ##
## Last update: 07/25/2015 (only documentation, by Augusto Garcia)     ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     #
#######################################################################

##' Compare all possible orders (exhaustive search) for a given sequence of
##' markers
##'
##' For a given sequence with \eqn{n}{n} markers, computes the multipoint
##' likelihood of all \eqn{\frac{n!}{2}}{n!/2} possible orders.
##'
##' Since the number \eqn{\frac{n!}{2}}{n!/2} is large even for moderate values
##' of \eqn{n}{n}, this function is to be used only for sequences with
##' relatively few markers. If markers were genotyped in an outcross population,
##' linkage phases need to be estimated and therefore more states need to be
##' visited in the Markov chain; when segregation types are D1, D2 and C,
##' computation can required a very long time (specially when markers linked in
##' repulsion are involved), so we recomend to use this function up to 6 or 7 markers.
##' For inbred-based populations, up to 10 or 11 markers can be ordered with this function,
##' since linkage phase are known.
##' The multipoint likelihood is calculated according to Wu et al.
##' (2002b) (Eqs. 7a to 11), assuming that the recombination fraction is the
##' same in both parents. Hidden Markov chain codes adapted from Broman et al.
##' (2008) were used.
##' 
##' @importFrom methods is
##' @importFrom utils setTxtProgressBar txtProgressBar tail head
##'
##' @param input.seq an object of class \code{sequence}.
##' @param n.best the number of best orders to store in object (defaults to
##' 50).
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @param verbose if \code{FALSE} (default), simplified output is displayed.
##' if \code{TRUE}, detailed output is displayed.
##' @return An object of class \code{compare}, which is a list containing the
##' following components: \item{best.ord}{a \code{matrix} containing the best
##' orders.} \item{best.ord.rf}{a \code{matrix} with recombination frequencies
##' for the corresponding best orders.} \item{best.ord.phase}{a \code{matrix}
##' with linkage phases for the best orders.} \item{best.ord.like}{a
##' \code{vector} with log-likelihood values for the best orders.}
##' \item{best.ord.LOD}{a \code{vector} with LOD Score values for the best
##' orders.} \item{data.name}{name of the object of class \code{onemap} with
##' the raw data.} \item{twopt}{name of the object of class \code{rf_2pts} with
##' the 2-point analyses.}
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{marker_type}} for details about segregation
##' types and \code{\link[onemap]{make_seq}}.
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci with
##' dominant and missing markers in various crosses from two inbred lines.
##' \emph{Genetica} 101: 47-58.
##'
##' Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J., Lincoln,
##' S. E. and Newburg, L. (1987) MAPMAKER: An interactive computer package for
##' constructing primary genetic linkage maps of experimental and natural
##' populations. \emph{Genomics} 1: 174-181.
##'
##' Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A.  A. F.
##' (2009) Evaluation of algorithms used to order markers on genetics maps.
##' _Heredity_ 103: 494-502.
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
##'   #outcrossing example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   markers <- make_seq(twopt,c(12,14,15,26,28))
##'   (markers.comp <- compare(markers))
##'   (markers.comp <- compare(markers,verbose=TRUE))
##'
##'   #F2 example
##'   data(onemap_example_f2)
##'   twopt <- rf_2pts(onemap_example_f2)
##'   markers <- make_seq(twopt,c(17,26,29,30,44,46,55))
##'   (markers.comp <- compare(markers))
##'   (markers.comp <- compare(markers,verbose=TRUE))
##'
##' }
##'
##'@export

compare<- function(input.seq,n.best=50,tol=10E-4,verbose=FALSE) {
    if(is(get(input.seq$data.name), "outcross"))
        return(compare_outcross(input.seq=input.seq,n.best=n.best,tol=tol,verbose=verbose))
    else
        return(compare_inbred(input.seq=input.seq,n.best=n.best,tol=tol,verbose=verbose))
}

## Compare all possible orders (exhaustive search) for a given sequence of
## markers (for outcrosses)
compare_outcross<- function(input.seq,n.best=50,tol=10E-4,verbose=FALSE)
{
    ## checking for correct objects
    if(!any(class(input.seq)=="sequence"))
        stop(sQuote(deparse(substitute(input.seq)))," is not an object of class 'sequence'")
    if(length(input.seq$seq.num) > 5)
        cat("WARNING: this operation may take a VERY long time\n")
    utils::flush.console()
    if(length(input.seq$seq.num) > 10) {
        cat("\nIt is not wise trying to use 'compare' with more than 10 markers \n")
        ANSWER <- readline("Are you sure you want to proceed? [y or n]\n")
        while(substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 1) != "y")
            ANSWER <- readline("\nPlease answer: 'y' or 'n' \n")
        if (substr(ANSWER, 1, 1) == "n") stop("Execution stopped!")
    }
    if(length(input.seq$seq.num) == 2)
        return(map(input.seq, tol=tol)) ## nothing to be done for 2 markers
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

        ## 'perm_pars' generates all n!/2 orders
        all.ord <- perm_pars(input.seq$seq.num)
        cat("\nComparing",nrow(all.ord),"orders:     \n\n")
        if (verbose){
            for(i in 1:nrow(all.ord)){
                ## print output for each order
                cat("Order", i, ":", all.ord[i,], "\n")
                utils::flush.console()
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
                Ph.Init <- comb_ger(phase.init)
                Rf.Init <- comb_ger(rf.init)
                if(nrow(Ph.Init)>1){
                    ##Removing ambigous phases
                    rm.ab<-rem_amb_ph(M=Ph.Init, w=input.seq, seq.num=all.ord[i,])
                    Ph.Init <- Ph.Init[rm.ab,]
                    Rf.Init <- Rf.Init[rm.ab,]
                    if(class(Ph.Init)=="integer"){
                        Ph.Init<-matrix(Ph.Init,nrow=1)
                        Rf.Init<-matrix(Rf.Init,nrow=1)
                    }
                }
                for(j in 1:nrow(Ph.Init)){
                    ## estimate parameters
                    final.map <- est_map_hmm_out(geno=t(get(input.seq$data.name, pos=1)$geno[,all.ord[i,]]),
		    error=t(get(input.seq$data.name, pos=1)$error[,all.ord[i,]]),
                                                 type=get(input.seq$data.name, pos=1)$segr.type.num[all.ord[i,]],
                                                 phase=Ph.Init[j,],
                                                 rf.vec=Rf.Init[j,],
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
            ## print output for each order
            ##    if (sum(i == round(out.pr))){
            ##      cat(rep("\b",nchar(nc)+1),sep="")
            ##      nc<-round(i*100/nrow(all.ord))
            ##      cat(nc,"%", sep="")
            ##      utils::flush.console()
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
            Ph.Init <- comb_ger(phase.init)
            Rf.Init <- comb_ger(rf.init)
            if(nrow(Ph.Init)>1){
                ##Removing ambigous phases
                rm.ab<-rem_amb_ph(M=Ph.Init, w=input.seq, seq.num=all.ord[i,])
                Ph.Init <- Ph.Init[rm.ab,]
                Rf.Init <- Rf.Init[rm.ab,]
                if(class(Ph.Init)=="integer"){
                    Ph.Init<-matrix(Ph.Init,nrow=1)
                    Rf.Init<-matrix(Rf.Init,nrow=1)
                }
            }
            for(j in 1:nrow(Ph.Init)){
                ## estimate parameters
                final.map <- est_map_hmm_out(geno=t(get(input.seq$data.name, pos=1)$geno[,all.ord[i,]]),
		error=t(get(input.seq$data.name, pos=1)$error[,all.ord[i,]]),
                                             type=get(input.seq$data.name, pos=1)$segr.type.num[all.ord[i,]],
                                             phase=Ph.Init[j,],
                                             rf.vec=Rf.Init[j,],
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
        cat("\n")
        best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),4)
        structure(list(best.ord = best.ord,
                       best.ord.rf = best.ord.rf,
                       best.ord.phase = best.ord.phase,
                       best.ord.like = best.ord.like,
                       best.ord.LOD = best.ord.LOD,
                       data.name=input.seq$data.name,
                       twopt=input.seq$twopt), class = "compare")

    }
}

## Compare all possible orders (exhaustive search) for a given sequence of
## markers (crosses derived from inbred lines)
compare_inbred<- function(input.seq,n.best=50,tol=10E-4,verbose=FALSE) {
    ## checking for correct objects
    if(!any(class(input.seq)=="sequence"))
        stop(sQuote(deparse(substitute(input.seq)))," is not an object of class 'sequence'")
    if(length(input.seq$seq.num) > 5)
        cat("WARNING: this operation may take a VERY long time\n")
    utils::flush.console()
    if(length(input.seq$seq.num) > 10) {
        cat("\nIt is not wisely trying to use 'compare' with more than 10 markers \n")
        ANSWER <- readline("Are you sure you want to proceed? [y or n]\n")
        while(substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 1) != "y")
            ANSWER <- readline("\nPlease answer: 'y' or 'n' \n")
        if (substr(ANSWER, 1, 1) == "n") stop("Execution stopped!")
    }
    if(length(input.seq$seq.num) == 2)
        return(map(input.seq, tol=tol)) ## nothing to be done for 2 markers
    else
    {
        ## allocating variables
        all.ord <- perm_pars(input.seq$seq.num)
        rf.init <- matrix(NA,nrow(all.ord),length(input.seq$seq.num)-1)
        best.ord <- matrix(NA,(n.best+1),length(input.seq$seq.num))
        best.ord.rf <- matrix(NA,(n.best+1),length(input.seq$seq.num)-1)
        best.ord.like <- best.ord.LOD <- rep(-Inf,(n.best+1))
        cat("\nComparing",nrow(all.ord),"orders:     \n\n")
        if(!verbose)
        {
            count <- 0
            pb <- txtProgressBar(style=3)
            setTxtProgressBar(pb, 0)
            cat("    ")
        }
        for(i in 1:nrow(all.ord))
        {
            ## print output for each order
            if (verbose) cat("Order", i, ":", all.ord[i,], "\n")
            utils::flush.console()
            seq.temp<-make_seq(get(input.seq$twopt), arg=all.ord[i,])
            seq.temp$twopt<-input.seq$twopt
            rf.temp<-get_vec_rf_in(seq.temp, acum=FALSE)
            final.map<-est_map_hmm_f2(geno=t(get(input.seq$data.name, pos=1)$geno[,all.ord[i,]]),
	    error=t(get(input.seq$data.name, pos=1)$error[,all.ord[i,]]),
                                      rf.vec=rf.temp,
                                      verbose=FALSE,
                                      tol=tol)
            best.ord[(n.best+1),] <- all.ord[i,]
            best.ord.rf[(n.best+1),] <- final.map$rf
            best.ord.like[(n.best+1)] <- final.map$loglike
            ## arrange orders according to the likelihood
            like.order <- order(best.ord.like, decreasing=TRUE)
            best.ord <- best.ord[like.order,]
            best.ord.rf <- best.ord.rf[like.order,]
            best.ord.like <- sort(best.ord.like, decreasing=TRUE)
            if(!verbose)
            {
                count<-count+1
                setTxtProgressBar(pb, count/nrow(all.ord))
            }
        }
        close(pb)
        ## cat("\nFinished\n\n")
        cat("\n")
        best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),4)
        structure(list(best.ord = best.ord,
                       best.ord.rf = best.ord.rf,
                       best.ord.phase = matrix(1, nrow(best.ord.rf), ncol(best.ord.rf)),
                       best.ord.like = best.ord.like,
                       best.ord.LOD = best.ord.LOD,
                       data.name=input.seq$data.name,
                       twopt=input.seq$twopt), class = "compare")
    }
}

## print method for object class 'compare'
##'@export
##'@method print compare
print.compare <- function(x,...) {
        FLAG<-0
        if(!is(get(x$data.name, pos=1), "outcross")) FLAG<-1

        phases.char <- c("CC","CR","RC","RR")
        n.ord <- max(which(head(x$best.ord.LOD,-1) != -Inf))
        unique.orders <- unique(x$best.ord[1:n.ord,])
        n.ord.nest <- nrow(unique.orders)
        phases.nested <- vector("list",n.ord.nest)
        LOD <- vector("list",n.ord.nest)
        if(!FLAG){
            for (i in 1:n.ord.nest) {
                same.order <- which(apply(x$best.ord[1:n.ord,],1,function(x) all(x==unique.orders[i,])))
                ifelse(length(same.order)==1,phases.nested[[i]] <- t(as.matrix(x$best.ord.phase[same.order,])),phases.nested[[i]] <- x$best.ord.phase[same.order,])
                LOD[[i]] <- x$best.ord.LOD[same.order]
            }
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
            cat(paste("order ",format(i,width=skip[1]),":  ",paste(format(unique.orders[i,],width=skip[2]),collapse=""), "     ",formatC(round(x$best.ord.LOD[i],2),digits=2,format="f",width=6),  "\n",sep=""))
        }
        cat(paste(rep("-",leng.print),collapse=""))
        cat("\n")
    }
    else stop ("Should not get here!")
}
## end of file
