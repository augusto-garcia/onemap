#' OneMap interface with polyRAD package
#'
#' Put here a description
#' 
#' @param vcf path to vcf file
#' @param onemap.obj object of onemap class
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR objetc
#' @param f1 f1 individual identification if F2 cross type
#' @param crosstype string defining the cross type, by now it supports only 
#' outcross and f2 intercross
#' @param global_error number from 0 to 1 defining the global error to be considered together 
#' with the genotype errors or the genotype probabilities or NULL to not considered any global error
#' @param use_genotypes_errors if \code{TRUE} the error probability of each genotype will be considered in emission function of HMM
#' @param use_genotype_probs if \code{TRUE} the probability of each possible genotype will be considered in emission function of HMM
#' 
#' @return onemap object with error updated 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' @seealso \code{\link[onemap]{extract_depth}} 
#'     \code{\link[onemap]{binom_error}} and 
#'     \url{https://github.com/dcgerard/updog}.
#'
#' @references 
#'
#' Clark LV, Lipka AE, and Sacks EJ (2019) Improvements to 
#' Genotype Calling in Polyploids Using the polyRAD R Package. 
#' Plant and Animal Genome Conference XXVII, January 12-16, 
#' San Diego, California, USA. doi:10.13140/RG.2.2.18358.75847
#' 
#' Clark LV, Lipka AE, and Sacks EJ (2018) polyRAD: Genotype 
#' Calling with Uncertainty from Sequencing Data in Polyploids 
#' and Diploids. Plant and Animal Genome Conference XXVI, 
#' January 13-17, San Diego, California, USA. doi:10.13140/RG.2.2.27134.08001
#'
#' @import polyRAD 
#'   
#' @export
polyRAD_genotype <- function(vcf=NULL, 
                          onemap.obj = NULL,
                          parent1=NULL,
                          parent2=NULL,
                          f1=NULL,
                          crosstype=NULL,
                          global_error = NULL,
                          use_genotypes_errors = TRUE,
                          use_genotypes_probs = FALSE,
                          rm_multiallelic = TRUE){
  # Do the checks
   poly.test <- VCF2RADdata(vcf, phaseSNPs = FALSE, 
                           min.ind.with.reads = 0,
                           min.ind.with.minor.allele = 0)
  if(crosstype=="f2 intercross"){
    poly.test <- SetDonorParent(poly.test, f1)
    poly.test <- SetRecurrentParent(poly.test, f1)
  } else if(crosstype=="outcross"){
    poly.test <- SetDonorParent(poly.test, parent1)
    poly.test <- SetRecurrentParent(poly.test, parent2)
  }
  
  mydata2 <- PipelineMapping2Parents(poly.test, 
                                     freqAllowedDeviation = 0.06,
                                     useLinkage = FALSE)
  
  
  seed.uniq <- sample(100000, 1)
  Export_MAPpoly(mydata2, paste0("temp.file.", seed.uniq)) 
  
  genotypes <- read.table(paste0("temp.file.",seed.uniq), skip=12)
  
  file.remove(paste0("temp.file.",seed.uniq))
  
  # this will change according to the vcf - bug!! Need attention!
  if(any(grepl(":", as.character(genotypes$V1)))){
    temp_list <- strsplit(as.character(genotypes$V1), split = "_")
    temp <- sapply(temp_list, "[", 1)
    pos <- gsub(":", "_", temp)
  } else {
    temp_list <- strsplit(as.character(genotypes$V1), split = "_")
    pos <- sapply(temp_list, function(x) if(length(x) > 2) paste0(x[1:2], collapse = "_"))
  }
  
  # Remove multiallelic markers
  multi <- names(which(table(pos) > onemap.obj$n.ind))
  if(length(multi) != 0){
    ## From vcf
    multi.idx <- which(pos %in% multi)
    genotypes <- droplevels(genotypes[-multi.idx,])
    pos <- pos[-multi.idx]
    
    # removing the multiallelics from the onemap object
    split.onemap <- split_onemap(onemap.obj, mks= multi)
    mult.obj <- split.onemap[[2]]
    onemap.obj <- split.onemap[[1]]
  }
  
  pos.onemap <- colnames(onemap.obj$geno)
  genotypes <- genotypes[which(pos%in%pos.onemap),]
  keep.mks <- which(pos.onemap%in%pos)
  
  # Remove parents
  if(crosstype=="f2 intercross"){
    genotypes <- genotypes[-which(genotypes[,2]%in%parent1),]
    genotypes <- genotypes[-which(genotypes[,2] %in%parent2),]
  }
  
  # Updating geno matrix
  onemap.obj$geno <- onemap.obj$geno[,keep.mks]
  
  new.geno <- maxpostprob <- vector()
  for(i in 1:dim(genotypes)[1]){
    if(which.max(genotypes[i,3:5]) == 3){
      new.geno[i] <- 3
    }else if(which.max(genotypes[i,3:5]) == 2){
      new.geno[i] <- 2
    } else if(which.max(genotypes[i,3:5]) == 1){
      new.geno[i] <- 1
    }
    maxpostprob[i] <- unlist(genotypes[i,3:5][which.max(genotypes[i,3:5])])
  }
  
  new.geno <- matrix(new.geno,nrow = onemap.obj$n.ind, ncol = length(keep.mks))
  maxpostprob <- matrix(maxpostprob,nrow = onemap.obj$n.ind, ncol = length(keep.mks))
  colnames(new.geno) <- colnames(maxpostprob) <- colnames(onemap.obj$geno)
  rownames(new.geno) <- rownames(maxpostprob) <- rownames(onemap.obj$geno)
  
  # Same order priority: individuals, markers
  genotypes$V2 <- factor(genotypes$V2, levels = genotypes$V2[1:onemap.obj$n.ind])
  
  genotypes <- genotypes[order(genotypes$V2),]
  
  # Print how many genotypes changed
  cat("This approach changed", (1- sum(new.geno == onemap.obj$geno)/length(new.geno))*100,"% of the genotypes\n")
  
  onemap.obj$geno <- new.geno
  
  # Removing markers
  onemap.obj$n.mar <- length(keep.mks)
  onemap.obj$segr.type <- onemap.obj$segr.type[keep.mks]
  onemap.obj$segr.type.num <- onemap.obj$segr.type.num[keep.mks]
  onemap.obj$CHROM <- onemap.obj$CHROM[keep.mks]
  onemap.obj$POS <- onemap.obj$POS[keep.mks]
  
  probs <- as.matrix(genotypes[,3:5])
  
  if(use_genotypes_probs){
    onemap.obj.new <- create_probs(onemap.obj = onemap.obj,
                                     genotypes_probs = probs,
                                     global_error = global_error)
  } else if(use_genotypes_errors){
    onemap.obj.new <- create_probs(onemap.obj = onemap.obj,
                                     genotypes_errors = 1- maxpostprob,
                                     global_error = global_error)
  } else if(!is.null(global_error)){
    onemap.obj.new <- create_probs(onemap.obj = onemap.obj,
                                     global_error = global_error)
  }
  
  if(!rm_multiallelic){
    if(length(multi) > 0)
      onemap.obj.new <- combine_onemap(onemap.obj.new, mult.obj)
  }
  
  return(onemap.obj.new)
}
