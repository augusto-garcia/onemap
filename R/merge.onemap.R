merge.onemap <- function(...) {
    onemap.objs <- list(...)
    n.objs <- length(onemap.objs)
    if (!n.objs) {
        stop("You must provide a list of OneMap objects as input.")
    }
    for (i in 1:n.objs) {
        if(class(onemap.objs[[i]])[1] != "onemap")
            stop("All objects must be of class 'onemap'.")
    }
    if (n.objs == 1) {
        stop("Nothing to merge.")
    }
    
    ## Check if all objects are of the same cross type
    crosstype <- class(onemap.objs[[1]])[2]
    for (i in 2:n.objs) {
        if(class(onemap.objs[[i]])[2] != crosstype)
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
    }
    else {
        n.ind <- length(sampleIDs)
    }

    ## Allocate
    geno <- matrix(0, nrow = n.ind, ncol = n.mar)
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
    }
    else {
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
            ## Find individual indices
            ind.matches <- match(rownames(onemap.objs[[i]]$geno), rownames(geno))
        }
        geno[ind.matches, mrk.start:mrk.end] <- onemap.objs[[i]]$geno
        colnames(geno)[mrk.start:mrk.end] <- colnames(onemap.objs[[i]]$geno)
############################ MAYBE CHECK FOR DUPLICATE MARKER NAMES?
        
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
            pheno[,phe.start:phe.end] <- onemap.objs[[i]]$pheno
            colnames(pheno)[phe.start:phe.end] <- colnames(onemap.objs[[i]]$pheno)
        }

        mrk.start <- mrk.start + cur.n.mar
        phe.start <- phe.start + cur.n.phe
    }
    
    if (all(is.na(CHROM))) {
        CHROM <- NULL
    }
    if (all(is.na(POS))) {
        POS <- NULL
    }
    
    ## Return "onemap" object
    input <- "merged"
    structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                   segr.type = segr.type, segr.type.num = segr.type.num,
                   n.phe = n.phe, pheno = pheno, CHROM = CHROM, POS = POS,
                   input = input),
              class = c("onemap", crosstype))
}
