#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: onemap_read_vcfR.R                                            #
# Contains: onemap_read_vcfR write_onemap_raw                         #
#                                                                     #
# Written by Cristiane Hayumi Taniguti                                #
#                                                                     #
# First version: 09/01/2018                                           #
# Last update: 09/01/2018                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Convert vcfR object to onemap object
##'
##' Converts data from a vcfR package to onemap initial object, while trying to identify 
##' the appropriate marker segregation patterns.
##'
##' Only biallelic SNPs and indels for diploid variant sites are considered.
##'
##' Genotype information on the parents is required for all cross types. For
##' full-sib progenies, both outbred parents must be genotyped. For backcrosses,
##' F2 intercrosses and recombinant inbred lines, the \emph{original inbred
##' lines} must be genotyped. Particularly for backcross progenies, the
##' \emph{recurrent line must be provided as the first parent} in the function
##' arguments.
##'
##' Marker type is determined based on parental genotypes. Variants for which parent
##' genotypes cannot be determined are discarded.
##'
##' Reference sequence ID and position for each variant site are also stored.
##'
##' @param vcfR.object object of class 'vcfR' from vcfR package.
##' @param cross type of cross. Must be one of: \code{"outcross"} for full-sibs;
##' \code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
##' \code{"ri self"} for recombinant inbred lines by self-mating; or
##' \code{"ri sib"} for recombinant inbred lines by sib-mating.
##' @param parent1 \code{string} specifying sample ID of the first parent.
##' @param parent2 \code{string} specifying sample ID of the second parent.
##' @param f1 \code{string} if you are working with f2 intercross or backcross populations you may have f1 parents in you vcf, specify its ID here
##' @param only_biallelic if TRUE (default) only biallelic markers are considered, if FALSE multiallelic markers are included.
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' @examples
##' 
##' \dontrun{
##' vcfR.object <- read.vcfR(system.file("extdata/vcf_example_out.vcf", package = "onemap"))
##' data <- onemap_read_vcfR(vcfR.object=vcfR.object,
##'                  cross="outcross",
##'                  parent1=c("P1"),
##'                  parent2=c("P2"))
##' }
##'                 
##'@export                  

onemap_read_vcfR <- function(vcfR.object=NULL,
                             cross = c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"),
                             parent1 =NULL,
                             parent2 =NULL,
                             f1=NULL,
                             only_biallelic = TRUE){
  
  if (is.null(vcfR.object)) {
    stop("You must specify one vcfR object.")
  }
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify samples as parents 1 and 2.")
  }
  if(!is(vcfR.object,"vcfR")){
    stop("You must specify one vcfR object.")
  }
  
  vcf <- vcfR.object
  n.mk <- dim(vcf@gt)[1]
  n.ind <- dim(vcf@gt)[2]-1
  INDS <- dimnames(vcf@gt)[[2]][-1]
  CHROM <- vcf@fix[,1]
  POS <- as.numeric(vcf@fix[,2])
  
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcf@gt)[[2]]==parent1) 
  P2 <- which(dimnames(vcf@gt)[[2]]==parent2) 
  
  if(is.vector(vcf@gt)){
    jump <- 1
  } else if(dim(vcf@gt)[1] == 0){
    jump <- 1
  } else jump <- 0
  
  if(jump == 1){
    warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
    
    onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcf@gt)[[2]]==parent1) -1 
  P2 <- which(dimnames(vcf@gt)[[2]]==parent2) -1
  
  MKS <- vcf@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcf@fix[,1],"_", vcf@fix[,2])
  
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcf@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcf@gt[,i], split=":"), "[[", GT))
  
  
  if(length(P1)==0 | length(P2)==0) stop("One or both parents names could not be found in your data")
  
  # Bugfix
  if(!only_biallelic & length(grep("[|]", GT_matrix[,c(P1,P2)])) > 0){
    all_data <- GT_matrix
    all_pos <- POS
    all_chrom <- CHROM
    all_mks <- MKS
    temp_pos <- temp_chrom <- temp_mks <- vector()
    temp_matrix <- data.frame()
    contigs <- unique(CHROM)
    # garantee that is the same contig
    for(w in 1:length(contigs)){
      CHROM <- all_chrom
      idx <- which(CHROM == contigs[w]) 
      GT_matrix <- all_data[idx,]
      POS <- all_pos[idx]
      MKS <- all_mks[idx]
      if(is(GT_matrix, "matrix")) phased <- grep("[|]", GT_matrix[,P1]) else phased <- grep("[|]", GT_matrix[P1]) 
      idx <- which(phased[-1] - phased[-length(phased)] ==1)
      idx.tot <- unique(sort(c(idx, idx +1)))
      idx.p1 <- phased[idx.tot]
      if(is(GT_matrix, "matrix")) phased <- grep("[|]", GT_matrix[,P2]) else phased <- grep("[|]", GT_matrix[P2])
      idx <- which(phased[-1] - phased[-length(phased)] ==1)
      idx.tot <- unique(sort(c(idx, idx +1)))
      idx.p2 <- phased[idx.tot]
      idx.tot <- unique(sort(c(idx.p1, idx.p2)))
      
      if(length(idx.tot)>0){
        # Filt NAs unphased heterozygotes
        idx.filt <- which(grepl(GT_matrix[idx.tot,P1], pattern = "[.]") |  grepl(GT_matrix[idx.tot,P2],pattern = "[.]"))
        if(length(idx.filt) > 0) idx.tot <- idx.tot[-idx.filt]
        idx.filt <- which(GT_matrix[idx.tot,P1] == "0/1" |  GT_matrix[idx.tot,P2] == "0/1")
        if(length(idx.filt) > 0) idx.tot <- idx.tot[-idx.filt]
        idx.filt <- which(GT_matrix[idx.tot,P1] == "0/1" |  GT_matrix[idx.tot,P2] == "0/1")
        if(length(idx.filt) > 0) idx.tot <- idx.tot[-idx.filt]
        idx.filt <- which(GT_matrix[idx.tot,P1] == "0|0" &  GT_matrix[idx.tot,P2] == "0|0")
        if(length(idx.filt) > 0) idx.tot <- idx.tot[-idx.filt]
        idx.filt <- which(GT_matrix[idx.tot,P1] == "1|1" |  GT_matrix[idx.tot,P2] == "1|1")
        if(length(idx.filt) > 0) idx.tot <- idx.tot[-idx.filt]
        
        idx <- which(idx.tot[-1] - idx.tot[-length(idx.tot)] ==1)
        idx.tot2 <- unique(sort(c(idx, idx +1)))
        idx.tot <- idx.tot[idx.tot2]
        
        if(length(idx.tot)>0){
          #list with haplo
          mnps.num <- split(idx.tot, cumsum(c(1, diff(idx.tot) != 1)))
          mnps <- lapply(mnps.num, function(x) GT_matrix[x,])
          GT_matrix <- GT_matrix[-idx.tot,]
          pos.mnps <- lapply(mnps.num, function(x) POS[x])
          mk.mnps <- lapply(mnps.num, function(x) MKS[x])
          POS <- POS[-idx.tot]
          CHROM <- CHROM[-idx.tot]
          MKS <- MKS[-idx.tot]
          mnp_matrix <- data.frame()
          mnp_pos <- mnp_chrom <- mnp_mk <- vector()
          for(i in 1:length(mnps)){
            if(sum(mnps[[i]] == ".",na.rm = T) > 0)
              mnps[[i]][which(mnps[[i]] == ".")] <- "./." #Techical issue
            temp <- lapply(apply(mnps[[i]],2, function(x) strsplit(x, "[| /]")), function(x) do.call(rbind, x))
            alleles <- sapply(temp, function(x) apply(x,2, function(y) paste0(y,collapse = "_")))
            p.alleles <- sort(unique(as.vector(alleles[,c(P1,P2)])))
            for(j in 1:length(p.alleles)){
              alleles[which(alleles == p.alleles[j])] <- j -1 # We deal with the progeny missing genotypes after
            }
            # Haplotypes found in progeny that are not present in parents are considered missing data
            mnp_matrix <- rbind.data.frame(mnp_matrix, t(apply(alleles, 2, function(x) paste0(x, collapse = "/"))), stringsAsFactors = F)
            mnp_pos <- c(mnp_pos, min(pos.mnps[[i]]))
            mnp_chrom <- c(mnp_chrom, contigs[[w]])
            mnp_mk <- c(mnp_mk, mk.mnps[[i]][which.min(pos.mnps[[i]])]) 
          }
          mnp_matrix <- as.matrix(mnp_matrix)
          mnp_matrix[grep(mnp_matrix, pattern =  "_")] <- "./."
          POS <- c(POS, mnp_pos)
          idx <- order(POS)
          POS <- POS[idx]
          GT_matrix <- rbind(GT_matrix, mnp_matrix)
          GT_matrix <- GT_matrix[idx,]
          CHROM <- c(CHROM, mnp_chrom)
          MKS <- c(MKS, mnp_mk)
        }
        # Markers not joint will be kept for next steps
        temp_matrix <- rbind.data.frame(temp_matrix, GT_matrix, stringsAsFactors = F)
        temp_pos <- c(temp_pos, POS)
        temp_chrom <- c(temp_chrom, CHROM)
        temp_mks <- c(temp_mks, MKS)
      }
    }
    rm(all_data)
    GT_matrix <- temp_matrix
    POS <- temp_pos
    CHROM <- temp_chrom
    MKS <- temp_mks
  }
  
  if(any(grepl("[|]", GT_matrix))){
    GT_matrix <- gsub("[|]", "/", as.matrix(GT_matrix))
    GT_matrix[which(GT_matrix == "1/0")] <- "0/1"
    GT_matrix[which(GT_matrix == "2/0")] <- "0/2"
    GT_matrix[which(GT_matrix == "3/0")] <- "0/3"
    GT_matrix[which(GT_matrix == "2/1")] <- "1/2"
    GT_matrix[which(GT_matrix == "3/1")] <- "1/3"
    GT_matrix[which(GT_matrix == "3/2")] <- "2/3"
  }
  
  # keep only biallelic
  if(only_biallelic | cross != "outcross"){
    rm_multi <- which(apply(GT_matrix, 1, function(x) any(grepl("2", x))))
    if(length(rm_multi) > 0){
      GT_matrix <- GT_matrix[-rm_multi,]
      CHROM <- CHROM[-rm_multi]
      POS <- POS[-rm_multi]
      MKS <- MKS[-rm_multi]
    }
  }
  n.mk <- nrow(GT_matrix)
  
  mk.type <- mk.type.num <- rep(NA, n.mk)
  if (cross == "outcross"){
    P1_1 <- sapply(strsplit(GT_matrix[,P1], "/"), "[", 1)
    P1_2 <- sapply(strsplit(GT_matrix[,P1], "/"), "[", 2)
    P2_1 <- sapply(strsplit(GT_matrix[,P2], "/"), "[", 1)
    P2_2 <- sapply(strsplit(GT_matrix[,P2], "/"), "[", 2)
    
    # Marker types
    
    GT_parents <- cbind(P1_1, P1_2,P2_1, P2_2)
    idx <- which(P1_1 == "." | P2_1 == "." |  P1_2 == "." | P2_2 == ".")
    GT_parents[idx,] <- NA
    
    idx <- apply(GT_parents, 1, function(x) length(x) == length(unique(x)))
    mk.type[idx] <- "A.1"
    mk.type.num[idx] <- 1
    idx <- apply(GT_parents, 1, function(x) (length(x) -1) == length(unique(x)))
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "A.2"
    mk.type.num[idx][idx.sub] <- 1
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] == P2_2[idx])
    mk.type[idx][idx.sub] <- "D1.9"
    mk.type.num[idx][idx.sub] <- 6
    idx.sub <- which(P1_1[idx] == P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "D2.14"
    mk.type.num[idx][idx.sub] <- 7
    idx <- apply(GT_parents, 1, function(x) (length(x) -2) == length(unique(x)))
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "B3.7"
    mk.type.num[idx][idx.sub] <- 4
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] == P2_2[idx])
    mk.type[idx][idx.sub] <- "D1.10"
    mk.type.num[idx][idx.sub] <- 6
    idx.sub <- which(P1_1[idx] == P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "D2.15"
    mk.type.num[idx][idx.sub] <- 7
    
    # It informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./." |  GT_matrix[,P2] == "." | GT_matrix[,P1] == ".")
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
    
    idx <- which(P1_1 == P1_2 & P2_1 == P2_2)
    
    if (length(idx) > 0)
      cat( length(MKS[idx]), "Markers were removed from the dataset because both of parents are homozygotes, these markers are considered non-informative in outcrossing populations.\n")
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      P1_1 <- P1_1[-rm_mk]
      P1_2 <- P1_2[-rm_mk]
      P2_1 <- P2_1[-rm_mk]
      P2_2 <- P2_2[-rm_mk]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <-CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
    } 
    
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
      return(onemap.obj)
    }     
    # Codification for OneMap
    idx <- which(mk.type=="A.1" | mk.type=="A.2")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    cat <- paste0(P1_2[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 3
    cat <- paste0(P1_2[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 4
    
    idx <- which(mk.type=="B3.7")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    cat <- paste0(P1_2[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 3
    
    idx <- which(mk.type=="D1.10")
    idx.sub <- which(P1_1[idx] == P2_1[idx])
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx.sub <- which(P1_2[idx] == P2_1[idx])
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D1.9")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_2[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D2.15" )
    idx.sub <- which(P1_1[idx] == P2_1[idx])
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_2[idx][idx.sub])
    cat.rev <- paste0(P2_2[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx.sub <- which(P1_1[idx] == P2_2[idx])
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_2[idx][idx.sub])
    cat.rev <- paste0(P2_2[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D2.14")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    
    GT_matrix[grepl("/", GT_matrix)] <- 0
    GT_matrix[grepl("[.]", GT_matrix)] <- 0
  } else if(cross== "f2 intercross"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.B.1"
    mk.type[which(GT_matrix[,P1][i] == "1/1" & GT_matrix[,P2][i] == "0/0")] <- "A.H.B.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.\n") 
    
    idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
    } 
    
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
      return(onemap.obj)
    }
    
    # Codification for OneMap
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2
    
    idx <- which(mk.type=="A.H.B.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 3
    
    idx <- which(mk.type=="A.H.B.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 3
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    mk.type <- mk.type.num <- rep("A.H.B", n.mk)
    mk.type.num[mk.type=="A.H.B"] <- 4
    
  } else if(cross=="f2 backcross"){
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.H.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.\n") 
    
    idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
    } 
    
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
      return(onemap.obj)
    }
    
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2
    
    idx <- which(mk.type=="A.H.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 0
    
    idx <- which(mk.type=="A.H.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 0
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    
    mk.type <- mk.type.num <- rep("A.H", n.mk)
    mk.type.num[mk.type=="A.H"] <- 8
    
  } else if(cross=="ri self" || cross=="ri sib"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.B.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.B.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in RILs populations.\n") 
    
    idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
    if (length(idx) > 0)
      cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
    }
    
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
      return(onemap.obj)
    }
    
    # Onemap codification
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 0
    
    idx <- which(mk.type=="A.B.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 3
    
    idx <- which(mk.type=="A.B.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 3
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    mk.type <- mk.type.num <- rep("A.B", n.mk)
    mk.type.num[mk.type=="A.B"] <- 9
    
  }
  GT_matrix[is.na(GT_matrix)] <- 0
  
  if(is.vector(GT_matrix)){
    jump <- 1
  } else if(dim(GT_matrix)[1]==0){
    jump <- 1
  } else jump <- 0
  
  if(jump == 1){
    warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
    
    onemap.obj <- empty_onemap_obj(vcf, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Removing parents
  if(is.null(f1)){
    GT_matrix <- apply(GT_matrix[,-c(P1,P2)],2,as.numeric)
    rownames(GT_matrix)  <- MKS
    colnames(GT_matrix)  <-  INDS[-c(P1,P2)] 
  } else{
    F1 <- which(dimnames(vcf@gt)[[2]]==f1) - 1
    GT_matrix <- apply(GT_matrix[,-c(P1,P2,F1)],2,as.numeric)
    rownames(GT_matrix)  <- MKS
    colnames(GT_matrix)  <-  INDS[-c(P1,P2,F1)] 
  }
  
  # Removing duplicated markers
  dupli <- MKS[duplicated(MKS)]
  if(length(dupli)>0){
    n.rm.mks <- length(dupli)
    dupli <- unique(dupli)
    cat(paste("There are more than one marker with the same IDs:", paste(MKS[duplicated(MKS)], collapse = " "), "\nOnly the one with less missing data was kept."))
    for(w in 1:length(dupli)){
      temp_GT <- GT_matrix[MKS==dupli[w],]
      mis_count <- apply(temp_GT, 1, function(x) sum(x==0))
      discard <- temp_GT[-which.min(mis_count),]
      if(is(discard, "matrix")){
        for(j in 1:dim(discard)[1]){
          idx <- which(apply(GT_matrix, 1, function(x) all(x == discard[j,])))
          idx <- idx[MKS[idx] == dupli[w]][1]
          GT_matrix <- GT_matrix[-idx,]
          mk.type <- mk.type[-idx]
          mk.type.num <- mk.type.num[-idx] 
          CHROM <- CHROM[-idx]
          POS <- POS[-idx]
          MKS <- MKS[-idx]
        }
      } else {
        idx <- which(apply(GT_matrix, 1, function(x) all(x == discard)))
        idx <- idx[MKS[idx] == dupli[w]][1]
        GT_matrix <- GT_matrix[-idx,]
        mk.type <- mk.type[-idx]
        mk.type.num <- mk.type.num[-idx] 
        CHROM <- CHROM[-idx]
        POS <- POS[-idx]
        MKS <- MKS[-idx]
      }
    }
    n.mk <- n.mk - n.rm.mks
  }
  legacy_crosses <- setNames(c("outcross", "f2", "backcross", "riself", "risib"), 
                             c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"))
  
  onemap.obj <- structure(list(geno= t(GT_matrix),
                               n.ind = dim(GT_matrix)[2],
                               n.mar = n.mk,
                               segr.type = mk.type,
                               segr.type.num = as.numeric(mk.type.num),
                               n.phe = 0,
                               pheno = NULL,
                               CHROM = CHROM,
                               POS = POS,
                               input = "vcfR.object"),
                          class=c("onemap",legacy_crosses[cross]))
  
  new.onemap.obj <- create_probs(onemap.obj, global_error = 10^-5)
  return(new.onemap.obj)
}


##' Convert onemap object to onemap raw file
##' 
##' Converts onemap R object to onemap input file. The input file brings information about the mapping population:
##' First line: cross type, it can be "outcrossing", "f2 intercross", "f2 backcross", "ri self" or "ri sib".
##' Second line:  number of individuals, number of markers, presence (1) or absence (0) of chromossome and position of the markers, and number of phenotypes mesured.
##' Third line: Individuals/sample names; 
##' Followed lines: marker name, marker type and genotypes. One line for each marker.
##' Final lines: chromossome, position and phenotypes informations. 
##' See more about input file format at vignettes.
##' 
##' @param onemap.obj object of class `onemap`
##' 
##' @param file.name a character for the onemap raw file name. Default is "out.raw"
##' 
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' @examples
##' \dontrun{
##' data(onemap_example_out)
##' 
##' write_onemap_raw(onemap_example_out, file.name = "onemap_example_out.raw", cross="outcross")
##' }
##'@export                  
write_onemap_raw <- function(onemap.obj=NULL, 
                             file.name = "out.raw"){
  
  if(is(onemap.obj, "outcross")){
    cross <- "outcross"
  } else if(is(onemap.obj, "f2")){
    cross <- "f2 intercross"
  } else if(is(onemap.obj, "backcross")){
    cross <- "f2 backcross"
  } else if(is(onemap.obj, "riself")){
    cross <- "ri self"
  } else if(is(onemap.obj, "risib")){
    cross <- "ri sib"
  }
  
  fileConn<-file(file.name, "w")
  head1 <- paste("data type", cross)
  head2 <- paste(onemap.obj$n.ind,
                 onemap.obj$n.mar,
                 as.numeric(!is.null(onemap.obj$CHROM)),
                 as.numeric(!is.null(onemap.obj$POS)), 
                 onemap.obj$n.phe)
  ind.names <- rownames(onemap.obj$geno)
  if(is.null(ind.names))
    ind.names <- paste0("ID", 1:onemap.obj$n.ind)
  
  geno.mat <- onemap.obj$geno
  
  if(is(onemap.obj, "outcross")){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    
    idx <- which(onemap.obj$segr.type == "A.1")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ad"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "bc"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "bd"
    
    idx <- which(onemap.obj$segr.type == "A.2")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "ba"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "bc"
    
    idx <- which(onemap.obj$segr.type == "A.3")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "bc"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "A.4")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "o"
    
    idx <- which(onemap.obj$segr.type == "B1.5" | onemap.obj$segr.type == "B2.6" | onemap.obj$segr.type == "B3.7")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D1.9" | onemap.obj$segr.type == "D2.14")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "bc"
    
    idx <- which(onemap.obj$segr.type == "D1.10" | onemap.obj$segr.type == "D2.15")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    
    idx <- which(onemap.obj$segr.type == "D1.11"  | onemap.obj$segr.type == "D2.16")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D1.12" | onemap.obj$segr.type == "D2.17")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    
    idx <- which(onemap.obj$segr.type == "D1.13" | onemap.obj$segr.type == "D2.18" | onemap.obj$segr.type == "C.8")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "o"
  }
  if(is(onemap.obj, c("f2","backcross"))){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    
    idx <- which(onemap.obj$segr.type == "A.H.B" | onemap.obj$segr.type == "A.H")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D.B")
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "d"
    
    idx <- which(onemap.obj$segr.type == "C.A")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 5)] <- "c"
  }
  if(is(onemap.obj, c("riself", "risib"))){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    geno.mat[which(geno.mat == 1)] <- "a"
    geno.mat[which(geno.mat == 3)] <- "b"
  }
  
  mat <- data.frame(paste0(rep("*", onemap.obj$n.mar),colnames(onemap.obj$geno)), 
                    onemap.obj$segr.type, t(geno.mat))
  colnames(mat) <- rownames(mat) <- NULL
  mat <- apply(mat, 1, function(x) paste(x, collapse = " "))
  
  writeLines(c(head1, head2, paste(ind.names, collapse = " "), mat),con =  fileConn)
  
  if(onemap.obj$n.phe > 0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    fen <- paste0("*", paste(colnames(onemap.obj$pheno), apply(onemap.obj$pheno,2, function(x) paste(x, collapse = " "))))
    writeLines( fen, con = fileConn)
  }
  
  if(length(onemap.obj$CHROM)>0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    chrom <- paste(paste0("*", "CHROM"), paste(onemap.obj$CHROM, collapse = " "))
    writeLines(chrom, con = fileConn)
  }
  
  if(length(onemap.obj$POS)>0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    pos <- paste(paste0("*", "POS"), paste(onemap.obj$POS, collapse = " "))
    writeLines(pos, con = fileConn)
  }
  close(fileConn)
}


