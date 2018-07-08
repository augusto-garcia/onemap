#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: onemap_read_vcfR.R                                            #
# Contains: onemap_read_vcfR                                          #
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
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' @examples
##'
##' vcfR.object <- vcfR::read.vcfR(system.file("extdata/vcf_example_out.vcf", package = "onemap"))
##' data <- onemap_read_vcfR(vcfR.object=vcfR.object,
##'                  cross="outcross",
##'                  parent1=c("P1"),
##'                  parent2=c("P2"))
##'@export                  

onemap_read_vcfR <- function(vcfR.object=NULL,
         cross = c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"),
         parent1 =NULL,
         parent2 =NULL){
  
  if (is.null(vcfR.object)) {
    stop("You must specify one vcfR object.")
  }
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify samples as parents 1 and 2.")
  }
  if(class(vcfR.object)!="vcfR"){
    stop("You must specify one vcfR object.")
  }
  
  vcf <- vcfR.object
  n.mk <- dim(vcf@gt)[1]
  n.ind <- dim(vcf@gt)[2]-1
  MKS <- vcf@fix[,3]
  INDS <- dimnames(vcf@gt)[[2]][-1]
  
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcf@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcf@gt[,i], split=":"), "[[", GT))
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcf@gt)[[2]]==parent1) - 1
  P2 <- which(dimnames(vcf@gt)[[2]]==parent2) - 1
  
  mk.type <- rep(NA, n.mk)
  if (cross == "outcross"){
    # Marker types
    mk.type[which(GT_matrix[,P1] == "0/1" & GT_matrix[,P2] == "0/0")] <- "D1.10.1"
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/1")] <- "D2.15.1"
    mk.type[which(GT_matrix[,P1] == "0/1" & GT_matrix[,P2] == "1/1")] <- "D1.10.2"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/1")] <- "D2.15.2"
    mk.type[which(GT_matrix[,P1] == "0/1" & GT_matrix[,P2] == "0/1")] <- "B3.7"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed of the dataset because one or both of parents have no informed genotypes (are missing data)")
     
    idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0")|
                   (GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1") |
                   (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1")|
                   (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0"))
      
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed from the dataset because both of parents are homozygotes, these markers are considered non-informative in outcrossing populations.")
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- vcf@fix[,1][-rm_mk]
      POS <- as.numeric(vcf@fix[,2][-rm_mk])
      mk.type <- mk.type[-rm_mk]
    } else {CHROM <- vcf@fix[,1]
    POS <- as.numeric(vcf@fix[,2])}
    
    # Codification for OneMap
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2

    idx <- which(mk.type=="B3.7")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 3
    
    idx <- which(mk.type=="D1.10.1" | mk.type=="D2.15.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 0
    
    idx <- which(mk.type=="D1.10.2" | mk.type=="D2.15.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 0
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    mk.type.num <- mk.type
    mk.type.num[mk.type=="D1.10.1" | mk.type=="D1.10.2"] <- 6
    mk.type.num[mk.type=="D2.15.1" | mk.type=="D2.15.2"] <- 7
    mk.type.num[mk.type=="B3.7"] <- 4
    mk.type[mk.type.num==6] <- "D1.10"
    mk.type[mk.type.num==7] <- "D2.15"
    
  } else if(cross== "f2 intercross"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.B.1"
    mk.type[which(GT_matrix[,P1][i] == "1/1" & GT_matrix[,P2][i] == "0/0")] <- "A.H.B.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed of the dataset because one or both of parents have no informed genotypes (are missing data)")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.") 
    
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- vcf@fix[,1][-rm_mk]
      POS <- as.numeric(vcf@fix[,2][-rm_mk])
      mk.type <- mk.type[-rm_mk]
    } else {CHROM <- vcf@fix[,1]
    POS <- as.numeric(vcf@fix[,2])}
    
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
    mk.type.num[mk.type=="A.H.B"] <- 1
    
  } else if(cross=="f2 backcross"){
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.H.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed of the dataset because one or both of parents have no informed genotypes (are missing data)")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.") 
    
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- vcf@fix[,1][-rm_mk]
      POS <- as.numeric(vcf@fix[,2][-rm_mk])
      mk.type <- mk.type[-rm_mk]
    } else {CHROM <- vcf@fix[,1]
    POS <- as.numeric(vcf@fix[,2])}
    
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2
    
    idx <- which(mk.type=="A.H.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 0
        
    idx <- which(mk.type=="A.H.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 0
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    
    mk.type <- mk.type.num <- rep("A.H", n.mk)
    mk.type.num[mk.type=="A.H"] <- NA
    
  } else if(cross=="ri self" || cross=="ri sib"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.B.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.B.2"
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed of the dataset because one or both of parents have no informed genotypes (are missing data)")
    
    idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
    if (length(idx) > 0)
      cat("Markers", MKS[idx], "were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in RILs populations.") 
    
    
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- vcf@fix[,1][-rm_mk]
      POS <- as.numeric(vcf@fix[,2][-rm_mk])
      mk.type <- mk.type[-rm_mk]
    } else {CHROM <- vcf@fix[,1]
    POS <- as.numeric(vcf@fix[,2])}
    
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
    mk.type.num[mk.type=="A.B"] <- NA
  }
  
  # Removing parents
  
  GT_matrix <- apply(GT_matrix[,-c(P1,P2)],2,as.numeric)
  rownames(GT_matrix) <- MKS
  colnames(GT_matrix) <- INDS[-c(P1,P2)] 
  
  legacy_crosses <- setNames(c("outcross", "f2", "backcross", "riself", "risib"), 
                             c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"))
  structure(list(geno= t(GT_matrix),
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
}
