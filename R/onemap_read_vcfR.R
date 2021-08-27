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


##' Convert vcf file to onemap object
##'
##' Converts data from a vcf file to onemap initial object, while identify 
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
##' @param vcf string defining the path to VCF file;
##' @param cross type of cross. Must be one of: \code{"outcross"} for full-sibs;
##' \code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
##' \code{"ri self"} for recombinant inbred lines by self-mating; or
##' \code{"ri sib"} for recombinant inbred lines by sib-mating.
##' @param parent1 \code{string} specifying sample ID of the first parent. If f2 backcross population, define here the ID of the backcrossed parent.
##' @param parent2 \code{string} specifying sample ID of the second parent.
##' @param f1 \code{string} if you are working with f2 intercross or backcross populations you may have f1 parents in you vcf, specify its ID here
##' @param only_biallelic if TRUE (default) only biallelic markers are considered, if FALSE multiallelic markers are included.
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' 
##' 
##' @importFrom rebus number_range
##' @importFrom vcfR read.vcfR extract.gt masplit
##' 
##' @examples
##' \dontrun{
##' vcfR.object <- read.vcfR(system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"))
##' data <- onemap_read_vcfR(vcfR.object=vcfR.object,
##'                  cross="outcross",
##'                  parent1=c("P1"),
##'                  parent2=c("P2"))
##' }
##'                 
##'@export                  
onemap_read_vcfR <- function(vcf=NULL,
                             cross = c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"),
                             parent1 =NULL,
                             parent2 =NULL,
                             f1=NULL,
                             only_biallelic = TRUE){
  
  if (is.null(vcf)) {
    stop("You must specify one vcf file.")
  }
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify samples as parents 1 and 2.")
  }
  
  vcfR.obj <- read.vcfR(vcf)
  n.mk <- dim(vcfR.obj@gt)[1]
  n.ind <- dim(vcfR.obj@gt)[2]-1
  INDS <- dimnames(vcfR.obj@gt)[[2]][-1]
  CHROM <- vcfR.obj@fix[,1]
  POS <- as.numeric(vcfR.obj@fix[,2])
  
  if(is.vector(vcfR.obj@gt)){
    jump <- 1
  } else if(dim(vcfR.obj@gt)[1] == 0){
    jump <- 1
  } else jump <- 0
  
  if(jump == 1){
    warning("Input vcfR.objR object do not have markers. An empty object onemap will be generated.")
    onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcfR.obj@gt)[[2]]==parent1) -1 
  P2 <- which(dimnames(vcfR.obj@gt)[[2]]==parent2) -1
  
  MKS <- vcfR.obj@fix[,3]
  if (any(MKS == "." | is.na(MKS))) {
    MKS <- paste0(vcfR.obj@fix[,1],"_", vcfR.obj@fix[,2])
    # Add tag if is duplicated positions (split form of mnps)
    for(i in 2:length(MKS)) {
      if(MKS[i] == paste0(strsplit(MKS[i-1], "_")[[1]][1:2], collapse = "_")) {
        z <- z + 1
        MKS[i] <- paste0(MKS[i], "_",z)
      } else {
        z <- 0
      }
    }
  }
  
  # Geno matrix
  GT_matrix <- extract.gt(vcfR.obj)
  
  if(length(P1)==0 | length(P2)==0) stop("One or both parents names could not be found in your data")
  
  # This function do not consider phased genotypes
  GT_matrix[grep("[.]", GT_matrix)] <- "./."
  GT_names <- names(table(GT_matrix))
  
  phased <- any(grepl("[|]", GT_names))
  if(phased)
    GT_matrix <- gsub("[|]", "/", as.matrix(GT_matrix))
  
  GT_names <- names(table(GT_matrix))
  GT_names_up <- strsplit(GT_names, "/")
  max.alleles <- max(as.numeric(do.call(c, GT_names_up[-1])))
  
  if(phased){
    GT_names_up[[1]] <- 0 # avoiding warning
    GT_names_up <- sapply(GT_names_up, function(x) paste(sort(as.numeric(x)), collapse = "/"))
    GT_names_up[1] <- "./."
    only_diff <- which(GT_names_up != GT_names)
    repl <- GT_names_up[only_diff]
    sear <- GT_names[only_diff]
    for(i in 1:length(sear)){
      GT_matrix[which(GT_matrix == sear[i])] <- repl[i]
    }
  }
  
  # keep only biallelic
  if(only_biallelic | cross != "outcross"){
    rx <- number_range(2, max.alleles)
    rm_multi <- which(apply(GT_matrix, 1, function(x) any(grepl(rx, x))))
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
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
      return(onemap.obj)
    }
    GT_matrix <- as.matrix(GT_matrix)
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
    cat <- paste0(P1_1[idx], "/", P2_1[idx]) # 18
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx]) # 18
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx]) # 18
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx]) # 18
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    cat <- paste0(P1_2[idx], "/", P2_2[idx]) # 18
    cat.rev <- paste0(P2_2[idx], "/", P1_2[idx]) # 18
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 3
    
    idx <- which(mk.type=="D1.10")
    idx.sub <- which(P1_1[idx] == P2_1[idx])
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub]) # 6
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub]) # 6
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_1[idx][idx.sub]) # 6
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_2[idx][idx.sub]) # 6 
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
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.H.B.2"
    
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
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
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
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
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
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
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
    warning("Input vcfR.objR object do not have markers. An empty object onemap will be generated.")
    
    onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Removing parents
  if(is.null(f1)){
    GT_matrix <- apply(GT_matrix[,-c(P1,P2), drop=F],2,as.numeric)
    if(!is(GT_matrix, "matrix")) GT_matrix <- t(as.matrix(GT_matrix)) # If there is only one marker
    colnames(GT_matrix)  <-  INDS[-c(P1,P2)] 
  } else{
    F1 <- which(dimnames(vcfR.obj@gt)[[2]]==f1) - 1
    GT_matrix <- apply(GT_matrix[,-c(P1,P2,F1), drop=F],2,as.numeric)
    if(!is(GT_matrix, "matrix")) GT_matrix <- t(as.matrix(GT_matrix))
    colnames(GT_matrix)  <-  INDS[-c(P1,P2,F1)] 
  }
  rownames(GT_matrix)  <- MKS
  
  legacy_crosses <- setNames(c("outcross", "f2", "backcross", "riself", "risib"), 
                             c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"))
  
  onemap.obj <- structure(list(geno= t(GT_matrix),
                               n.ind = if(!is(GT_matrix, "matrix")) length(GT_matrix) else dim(GT_matrix)[2],
                               n.mar = n.mk,
                               segr.type = mk.type,
                               segr.type.num = as.numeric(mk.type.num),
                               n.phe = 0,
                               pheno = NULL,
                               CHROM = CHROM,
                               POS = POS,
                               input = "vcf"),
                          class=c("onemap",legacy_crosses[cross]))
  
  onemap.obj  <- rm_dupli_mks(onemap.obj)
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
##' write_onemap_raw(onemap_example_out, file.name = "onemap_example_out.raw")
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
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "b"
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "d"
    
    idx <- which(onemap.obj$segr.type == "C.A")
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "c"
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


