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
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' @examples
##' 
##' \dontrun{
##' vcfR.object <- vcfR::read.vcfR(system.file("extdata/vcf_example_out.vcf", package = "onemap"))
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
	 f1=NULL){
  
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
  INDS <- dimnames(vcf@gt)[[2]][-1]
 
  MKS <- vcf@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcf@fix[,1],"_", vcf@fix[,2])
    
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcf@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcf@gt[,i], split=":"), "[[", GT))
  
  if(any(grepl("|", GT_matrix))){
    GT_matrix[GT_matrix=="1|0"] <- GT_matrix[GT_matrix=="0|1"] <- "0/1"
    GT_matrix[GT_matrix=="0|0"] <- "0/0"
    GT_matrix[GT_matrix=="1|1"] <- "1/1"
  }
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcf@gt)[[2]]==parent1) - 1
  P2 <- which(dimnames(vcf@gt)[[2]]==parent2) - 1
  
  if(length(P1)==0 | length(P2)==0) stop("One or both parents names could not be found in your data")

  F1 <- which(dimnames(vcf@gt)[[2]]==f1) - 1

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
      cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
     
    idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0")|
                   (GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1") |
                   (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1")|
                   (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0"))
      
    if (length(idx) > 0)
      cat( length(MKS[idx]), "Markers were removed from the dataset because both of parents are homozygotes, these markers are considered non-informative in outcrossing populations.\n")
    
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

  if(is.null(f1)){
        GT_matrix <- apply(GT_matrix[,-c(P1,P2)],2,as.numeric)
        rownames(GT_matrix)  <- MKS
        colnames(GT_matrix)  <-  INDS[-c(P1,P2)] 
    } else{
        GT_matrix <- apply(GT_matrix[,-c(P1,P2,F1)],2,as.numeric)
        rownames(GT_matrix)  <- MKS
        colnames(GT_matrix)  <-  INDS[-c(P1,P2,F1)] 
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
  
  new.onemap.obj <- create_probs(onemap.obj)
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
##' @param cross a character describing the cross type. It can be "outcrossing", 
##' "f2 intercross", "f2 backcross", "ri self" or "ri sib"
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
                             file.name = "out.raw", 
                             cross = c("outcross", "f2 backcross", "f2 intercross", "ri self", "ri sib")){

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
  
  if(class(onemap.obj)[2] == "outcross"){
    
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
    if(class(onemap.obj)[2] == "f2" | class(onemap.obj)[2] == "backcross"){

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
    if(class(onemap.obj)[2] == "riself" || class(onemap.obj)[2] == "risib"){

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


