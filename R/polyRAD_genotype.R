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
#' @param tech.issue need to be fixed
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
                          tech.issue=TRUE,
                          global_error = 1){
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
  temp <- gsub(":", "_", as.character(genotypes$V1))
  temp_list <- strsplit(temp, split = "_")
  pos <- sapply(temp_list, function(x) if(length(x) > 2) paste0(x[1:2], collapse = "_"))
  
  
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
  
  new.geno <- vector()
  for(i in 1:dim(genotypes)[1]){
    if(which.max(genotypes[i,3:5]) == 3){
      new.geno[i] <- 3
    }else if(which.max(genotypes[i,3:5]) == 2){
      new.geno[i] <- 2
    } else if(which.max(genotypes[i,3:5]) == 1){
      new.geno[i] <- 1
    }
  }
  
  new.geno <- matrix(new.geno,nrow = onemap.obj$n.ind, ncol = length(keep.mks))
  colnames(new.geno) <- colnames(onemap.obj$geno)
  rownames(new.geno) <- rownames(onemap.obj$geno)
  
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
  
  probs <- as.matrix(genotypes[,3:5]*(1- global_error))
  probs[which(probs == 0)] <- global_error
  
  polyrad.one <- create_probs(onemap.obj = onemap.obj, genotypes_probs =  probs)
  
  return(polyrad.one)
}
