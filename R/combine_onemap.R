#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: combine_onemap.R                                              #
# Contains: combine_onemap split_onemap                               #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2016, Gabriel R A Margarido                           #
#                                                                     #
# First version: 01/11/2016                                           #
# Last update: 01/11/2016                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Combine OneMap datasets
##'
##' Merge two or more OneMap datasets from the same cross type. Creates an
##' object of class \code{onemap}.
##'
##' Given a set of OneMap datasets, all from the same cross type (full-sib,
##' backcross, F2 intercross or recombinant inbred lines obtained by self-
##' or sib-mating), merges marker and phenotype information to create a
##' single \code{onemap} object.
##'
##' If sample IDs are present in all datasets (the standard new format), not
##' all individuals need to be genotyped in all datasets - the merged dataset
##' will contain all available information, with missing data elsewhere. If
##' sample IDs are missing in at least one dataset, it is required that all
##' datasets have the same number of individuals, and it is assumed that they
##' are arranged in the same order in every dataset.
##'
##' @param ... Two or more \code{onemap} dataset objects of the same cross
##' type.
##' @return An object of class \code{onemap}, i.e., a list with the following
##' components: \item{geno}{a matrix with integers indicating the genotypes
##' read for each marker. Each column contains data for a marker and each row
##' represents an individual.} \item{n.ind}{number of individuals.}
##' \item{n.mar}{number of markers.} \item{segr.type}{a vector with the
##' segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
##' vector with the segregation type of each marker, represented in a
##' simplified manner as integers, i.e. 1 corresponds to markers of type
##' \code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
##' to markers of type \code{"B2.6"}; 4 corresponds to markers of type
##' \code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
##' to markers of type \code{"D1"} and 7 corresponds to markers of type
##' \code{"D2"}. Markers for F2 intercrosses are coded as 1; all other crosses
##' are left as \code{NA}.} \item{input}{a string indicating that this is a
##' combined dataset.} \item{n.phe}{number of phenotypes.} \item{pheno}{a
##' matrix with phenotypic values. Each column contains data for a trait and
##' each row represents an individual.}
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{\link[onemap]{read_onemap}} and
##' \code{\link[onemap]{read_mapmaker}}.
##' @references Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993)
##' Constructing genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial
##' and reference manual. \emph{A Whitehead Institute for Biomedical Research
##' Technical Report}.
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords IO
##' @examples
##'
##'     combined_data <- combine_onemap(onemap_data1, onemap_data2)
##'   
##'@export
combine_onemap <- function(...) {
    onemap.objs <- list(...)
    n.objs <- length(onemap.objs)
    if (!n.objs) {
        stop("You must provide a list of OneMap objects as input.")
    }
    for (i in 1:n.objs) {
        if(!is(onemap.objs[[i]], "onemap"))
            stop("All objects must be of class 'onemap'.")
    }
    
    n.mar.all <- sapply(onemap.objs, "[[",3)
    idx <- which(n.mar.all == 0)
    if(length(idx) > 0){
        warning(paste("Object", idx, "has no markers.\n"))
        onemap.objs <- onemap.objs[-idx]
        n.objs <- length(onemap.objs)
    }
    
    if (n.objs == 1) {
        warning("Nothing to merge.")
        return(onemap.objs[[1]])
    }
    
    ## Check if all objects are of the same cross type
    crosstype <- class(onemap.objs[[1]])[2]
    for (i in 2:n.objs) {
        if(!is(onemap.objs[[i]], crosstype))
            stop("All objects must be of the same cross type.")
    }
    
    ## Gather required information from each dataset
    n.mar <- 0
    n.ind <- 0
    sampleIDs <- NULL
    n.phe <- 0
    sampleID.flag <- FALSE
    for (i in 1:n.objs) {
        n.mar <- n.mar + onemap.objs[[i]]$n.mar
        n.phe <- n.phe + onemap.objs[[i]]$n.phe
        
        ## Get unique progeny individuals
        cur.sampleIDs <- rownames(onemap.objs[[i]]$geno)
        sampleIDs <- unique(c(sampleIDs, cur.sampleIDs))
        
        ## Check if sample IDs are missing in this dataset
        if(is.null(cur.sampleIDs)) {
            sampleID.flag <- TRUE
        }
    }
    if (sampleID.flag) {
        ## At least one dataset is missing sample IDs: ignore 'sampleIDs' and assume all objects have the same genotype structure
        n.ind <- onemap.objs[[1]]$n.ind
        for (i in 2:n.objs) {
            if(onemap.objs[[i]]$n.ind != n.ind)
                stop("Sample IDs are missing in at least one dataset. All objects must contain the same number of individuals and in the same order.")
        }
    }   else {
        n.ind <- length(sampleIDs)
    }
    
    ## Allocate
    geno <- matrix(0, nrow = n.ind, ncol = n.mar)
    error <- matrix(1, nrow= n.ind*n.mar, ncol= 4)
    rownames.error <- rep(0,n.ind*n.mar)
    colnames(geno) <- rep(NA, n.mar)
    if (!sampleID.flag) {
        rownames(geno) <- sampleIDs
    }
    segr.type <- rep(NA, n.mar)
    segr.type.num <- rep(NA, n.mar)
    CHROM <- rep(NA, n.mar)
    POS <- rep(NA, n.mar)
    if (n.phe) {
        pheno <- matrix(NA, nrow = n.ind, ncol = n.phe)
    }    else {
        pheno <- NULL
    }
    
    ## Merge data
    mrk.start <- 1
    phe.start <- 1
    
    for (i in 1:n.objs) {
        cur.n.mar <- onemap.objs[[i]]$n.mar
        mrk.end <- mrk.start + cur.n.mar - 1
        if (sampleID.flag) {
            ## We assume all progeny individuals are in the same order
            ind.matches <- 1:n.ind
        }
        else {
            ## Find progeny indices
            ind.matches <- match(rownames(onemap.objs[[i]]$geno), rownames(geno))
        }
        geno[ind.matches, mrk.start:mrk.end] <- onemap.objs[[i]]$geno
        colnames(geno)[mrk.start:mrk.end] <- colnames(onemap.objs[[i]]$geno)
        
        error_idx_cur <- rep(1:n.ind, each=cur.n.mar)
        error_idx_joint <- rep(1:n.ind, each=n.mar)
        for(w in 1:length(ind.matches)){
            error[which(error_idx_joint == ind.matches[w]),][mrk.start:mrk.end,] <- onemap.objs[[i]]$error[which(error_idx_cur==w),]
            rownames.error[which(error_idx_joint == ind.matches[w])[mrk.start:mrk.end]] <- rownames(onemap.objs[[i]]$error)[which(error_idx_cur==w)]
        }
        
        segr.type[mrk.start:mrk.end] <- onemap.objs[[i]]$segr.type
        segr.type.num[mrk.start:mrk.end] <- onemap.objs[[i]]$segr.type.num
        if (!is.null(onemap.objs[[i]]$CHROM)) {
            CHROM[mrk.start:mrk.end] <- onemap.objs[[i]]$CHROM
        }
        if (!is.null(onemap.objs[[i]]$POS)) {
            POS[mrk.start:mrk.end] <- onemap.objs[[i]]$POS
        }
        
        cur.n.phe <- onemap.objs[[i]]$n.phe
        phe.end <- phe.start + cur.n.phe - 1
        if (cur.n.phe) {
            pheno[ind.matches, phe.start:phe.end] <- onemap.objs[[i]]$pheno
            colnames(pheno)[phe.start:phe.end] <- colnames(onemap.objs[[i]]$pheno)
        }
        
        mrk.start <- mrk.start + cur.n.mar
        phe.start <- phe.start + cur.n.phe
    }
    
    rownames(error) <- rownames.error
    if (all(is.na(CHROM))) {
        CHROM <- NULL
    }
    if (all(is.na(POS))) {
        POS <- NULL
    }
    
    ## Return "onemap" object
    input <- "combined"
    onemap.obj <- structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                   segr.type = segr.type, segr.type.num = segr.type.num,
                   n.phe = n.phe, pheno = pheno, CHROM = CHROM, POS = POS,
                   input = input, error=error),
              class = c("onemap", crosstype))
    
    onemap.obj <- rm_dupli_mks(onemap.obj)
    return(onemap.obj)
}

#' Split onemap data sets
#' 
#' Receives one onemap object and a vector with markers names to be 
#' removed from the input onemap object and inserted in a new one. The output
#' is a list containing the two onemap objects.
#' 
#' @param onemap.obj object of class onemap
#' @param mks markers names (vector of characters) or number (vector of integers) to be removed and added to a new onemap object
#' 
#' @return  a list containing in first level the original onemap object without 
#' the indicated markers and the second level the new onemap object with only 
#' the indicated markers 
#' 
#' @export
split_onemap <- function(onemap.obj=NULL, mks=NULL){
    
    if(!is(onemap.obj, c("onemap"))) stop("Input object must be of class onemap")

    if(is(mks, "character")){
        idx.mks <- which(colnames(onemap.obj$geno) %in% mks)
        rev.mks <- which(!colnames(onemap.obj$geno) %in% mks)
        if(any(!mks %in% colnames(onemap.obj$geno))) 
            stop("One or more of the selected markers do not exist in the onemap object")
            
    } else if(is(mks, c("numeric"))){
        idx.mks <- mks
        if(any(mks > onemap.obj$n.mar)) 
            stop("One or more of the selected markers do not exist in the onemap object")
        
        idx.temp <- 1:onemap.obj$n.mar
        rev.mks <- idx.temp[-mks]
    }
    
    new.obj <- onemap.obj
    new.obj$geno <- onemap.obj$geno[,idx.mks]
    new.obj$n.mar <- length(idx.mks)
    new.obj$segr.type <- onemap.obj$segr.type[idx.mks]
    new.obj$segr.type.num <- onemap.obj$segr.type.num[idx.mks]
    new.obj$CHROM <- onemap.obj$CHROM[idx.mks]
    new.obj$POS <- onemap.obj$POS[idx.mks]
    new.obj$error <- onemap.obj$error[idx.mks + rep(c(0:(onemap.obj$n.ind-1))*onemap.obj$n.mar, each=length(idx.mks)),]
    return(new.obj)
}

