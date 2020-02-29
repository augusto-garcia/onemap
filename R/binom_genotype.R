#' Estimates genotype error probabilities using binomial distribution for allele counts
#' 
#' Uses allele counts for each genotype from vcf, generated with NGS genotyping, to estimate 
#' genotype error probability. These probabilities can be used in multipoint approach to 
#' weighted the markers according with their sequecing quality
#' 
#' @param depth_matrix list containing ref allele counts for parents and progeny, output 
#' from extract_depth function
#' @param mean_phred mean quality phred value of the sequencing 
#' @return onemap object with error updated 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' 
#' @references 
#' Bilton,T. P. , Schofield,  M. R. , Black,  M. A. , Chagné, D. , Wilcox, P. L. 
#' and Dodds, K. G. (2018) Accounting for Errors in Low Coverage High-Throughput
#' Sequencing Data When Constructing Genetic Maps Using Biparental Outcrossed Populations.
#' \emph{Genetics}. ISSN 0016-6731.
#' @export
binom_genotype <- function(vcfR.object=NULL,
                        onemap.object= NULL,
                        vcf.par = c("AD", "DPR"),
                        parent1="P1",
                        parent2="P2",
                        recovering = FALSE,
                        mean_phred=20,
                        depths = NULL){
  if(is.null(depths)){
    depth_matrix <- extract_depth(vcfR.object, 
                                  onemap.object, 
                                  vcf.par,
                                  parent1,
                                  parent2,
                                  recovering)
  } else {
    depth_matrix <- depths
  }
  
 # Extract list objects
    for(i in 1:length(depth_matrix)){
        tempobj=depth_matrix[[i]]
        eval(parse(text=paste(names(depth_matrix)[[i]],"= tempobj")))
    }
    
    # Select the number of less frequent allele
    minor_matrix <- matrix(rep(NA, n.mks*(n.ind+2)), nrow = n.mks, ncol = n.ind+2)
    alt <- cbind(palt, oalt)
    ref <- cbind(pref, oref)
    size <- cbind(psize, osize)
    minor_matrix[which(ref >= alt)] <- alt[which(ref >= alt)]
    minor_matrix[which(ref < alt)] <- ref[which(ref < alt)]
    
    # Probabilities calculation by binomial distribution
    size[which(size==0)] <- NA
    het_matrix <- choose(size, minor_matrix)*(0.5^minor_matrix)*(0.5^(size-minor_matrix))
    homo_matrix <- choose(size, minor_matrix)*((10^(-mean_phred/10))^minor_matrix)*((1-(10^((-mean_phred/10))))^(size-minor_matrix))
    homo.alt <- choose(size, minor_matrix)*((10^(-mean_phred/10))^(size-minor_matrix))*((1-(10^((-mean_phred/10))))^minor_matrix)

    # Reviewed matrix 
    check_matrix <- matrix(rep(NA, n.mks*(n.ind+2)), nrow = n.mks, ncol = n.ind+2)
    check_matrix[which(het_matrix >= homo_matrix)] <- 2
    check_matrix[which(het_matrix == Inf)] <- 2
    idx <- which(het_matrix < homo_matrix) 
    check_matrix[idx][ref[idx] > alt[idx]] <- 1
    check_matrix[idx][ref[idx] < alt[idx]] <- 3
    check_matrix[is.na(check_matrix)] <- 0
    
    # Correcting genotypes according to segregation
    # Removing missing and non-informative
    rm.mis <- which(check_matrix[,1]==0 | check_matrix[,2]==0)
    rm.mis1 <- vector()
    for(i in 1:dim(check_matrix)[1]){
      if(length(table(check_matrix[i,]))==1)
        rm.mis1 < c(rm.mis1,i)
    }
    rm.mis <- sort(unique(c(rm.mis, rm.mis1)))
    
    rm.mk <- which(check_matrix[,1] == 1 & check_matrix[,2] == 1 | check_matrix[,1] == 1 & check_matrix[,2] == 3 | 
                     check_matrix[,1] == 3 & check_matrix[,2] == 1 | 
                     check_matrix[,1] == 3 & check_matrix[,2] == 3)
    
    if(length(rm.mis) > 0 | length(rm.mk) > 0){
      rm.tot <- unique(sort(c(rm.mis,rm.mk)))
    } else { rm.tot <- vector()}
    
    if (length(rm.tot) > 0){
      check_matrix <- check_matrix[-rm.tot,]
      het_matrix <- het_matrix[-rm.tot,]
      homo_matrix <- homo_matrix[-rm.tot,]
      homo.alt <- homo.alt[-rm.tot,]
      n.mks <- n.mks - length(rm.tot)
      mks <- mks[-rm.tot]
    }
    
    segr.type <- segr.type.num <- rep(NA, n.mks)
    
    pcheck <- check_matrix[,c(1,2)]
    ocheck <- check_matrix[,-c(1,2)]
    het_matrix <- het_matrix[,-c(1,2)]
    homo_matrix <- homo_matrix[,-c(1,2)]
    homo.alt <- homo.alt[,-c(1,2)]
    
    #D2.15
    idx <- which(pcheck[,1] == 1 & pcheck[,2] == 2)
    ocheck[idx,][ocheck[idx,]==3] <- 0
    segr.type[idx] <- "D2.15"
    segr.type.num[idx] <- 7
    idx <- which(pcheck[,1] == 3 & pcheck[,2] == 2)
    ocheck[idx,][ocheck[idx,]==1] <- 0
    ocheck[idx,][ocheck[idx,]==3] <- 1
    segr.type[idx] <- "D2.15"
    segr.type.num[idx] <- 7
    
    #D1.10
    idx <- which(pcheck[,1] == 2 & pcheck[,2] == 1)
    ocheck[idx,][ocheck[idx,]==3] <- 0
    segr.type[idx] <- "D1.10"
    segr.type.num[idx] <- 6
    idx <- which(pcheck[,1] == 2 & pcheck[,2] == 3)
    ocheck[idx,][ocheck[idx,]==1] <- 0
    ocheck[idx,][ocheck[idx,]==3] <- 1
    segr.type[idx] <- "D1.10"
    segr.type.num[idx] <- 6
    
    #B3.7
    idx <- which(pcheck[,1] == 2 & pcheck[,2] == 2)
    segr.type[idx] <- "B3.7"
    segr.type.num[idx] <- 4
    
    comp <- which(colnames(onemap.object$geno) %in% mks)
    
    # If some marker in onemap object are now non-informative and didn't recover any marker
    if(length(colnames(onemap.object$geno)) - length(comp) > 0 & length(mks) == length(comp)){
      cat(length(colnames(onemap.object$geno)) - length(comp),
          "markers of original onemap object were considered non-informative after new SNP calling and were removed from analysis \n")
      cat("This approach changed", (sum(t(ocheck) != onemap.object$geno[comp])/length(onemap.object$geno[comp]))*100, "% of the genotypes\n")  
      
      # If some marker in onemap object are now non-informative and some marker were recoved from vcf
    } else if(length(colnames(onemap.object$geno)) - length(comp) > 0 & length(mks) > length(comp)){
      cat(length(mks) - length(colnames(onemap.object$geno)[comp]), "markers were recovered from vcf file and added to onemap object \n")
      comp1 <- which(mks %in% colnames(onemap.object$geno))
      cat("This approach changed", (sum(t(ocheck)[comp1] != onemap.object$geno[comp])/length(onemap.object$geno[comp]))*100, "% of the genotypes\n")  
      
      # If any marker in onemap object were considered non-informative and some marker were recovered from vcf
    } else if(length(colnames(onemap.object$geno)) - length(comp) == 0 & length(mks) > length(comp)){
      comp1 <- which(mks %in% colnames(onemap.object$geno))
      cat("This approach changed", (sum(t(ocheck)[comp1] != onemap.object$geno)/length(onemap.object$geno))*100, "% of the genotypes\n") 
    } else{
      cat("This approach changed", (sum(t(ocheck) != onemap.object$geno)/length(onemap.object$geno))*100, "% of the genotypes\n")
    }
    
    cat("New onemap object contains", length(mks), "markers")
    
   # error prob is the probability of being homozygote when heterozigote genotype is inferred
    error_matrix <- matrix(rep(NA, n.mks*n.ind), nrow = n.mks)
    error_matrix[which(ocheck==2)] <- homo_matrix[which(ocheck==2)] + homo.alt[which(ocheck==2)] 
    error_matrix[which(ocheck==1 | ocheck==3)] <- het_matrix[which(ocheck==1 | ocheck==3)] + homo.alt[which(ocheck==1 | ocheck==3)]
    error_matrix[which(is.na(error_matrix))] <- 10^((-mean_phred/10))
    error_matrix[which(ocheck==0)] <- 10^-5

    geno <- t(ocheck)
    error <- t(error_matrix)
    colnames(geno) <- colnames(error) <- mks
    rownames(geno) <- rownames(geno) <- inds

    
    onemap_binom <- onemap.object
    onemap_binom$geno  <- geno 
    onemap_binom$error <- error
    onemap_binom$n.ind <- length(inds)
    onemap_binom$n.mar <- length(mks)
    onemap_binom$segr.type.num <- segr.type.num
    onemap_binom$segr.type <- segr.type
    if(length(rm.tot) > 0){
      onemap_binom$CHROM <- CHROM[-rm.tot]
      onemap_binom$POS <- POS[-rm.tot]
    }
    structure(onemap_binom)
}
    
