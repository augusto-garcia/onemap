#' OneMap interface with polyRAD package
#'
#' Put here a description
#' 
#' @param 
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR objetc
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
polyRAD_error <- function(vcf=NULL, 
                          onemap.obj = NULL,
                          parent1=NULL,
                          parent2=NULL,
                          f1=NULL,
                          crosstype=NULL,
                          tech.issue=TRUE,
                          depths= NULL){
  # Do the checks
  
  if(!is.null(depths)){
    # The input od polyRAD need to be a VCF, then this part takes the allele depth from "depths" and put at AD field of input vcf
    idx <- system(paste0("grep -in 'CHROM' ", vcf), intern = T) # This part only works in linux OS
    idx.i <- strsplit(idx, split = ":")[[1]][1]
    seed <- sample(1:10000, 1)
    system(paste0("head -n ", idx.i," ", vcf, " > head.",seed))
    
    vcf.tab <- read.table(vcf, stringsAsFactors = F)
    vcf.init <- vcf.tab[,1:8]
    AD.colum <- rep("AD", dim(vcf.init)[1])
    vcf.init <- cbind(vcf.init, AD.colum)
    
    rs <- rownames(depths[[1]])
    vcf.init[,3] <- rs
    
    ind.n <- colnames(depths[[1]]) # The names came in different order
    
    header <- strsplit(idx, split = "\t")[[1]]
    ind.vcf <- header[10:length(header)]
    ind.n <- factor(ind.n, levels = ind.vcf)
    
    depths[[1]] <- depths[[1]][,order(ind.n)]
    depths[[2]] <- depths[[2]][,order(ind.n)]
    
    comb.depth <- matrix(paste0(as.matrix(depths[[1]]), ",", as.matrix(depths[[2]])), ncol = ncol(depths[[2]]))
    colnames(comb.depth) <- ind.vcf
    #hmc.file <- cbind(rs, comb.depth)
    
    vcf.body <- cbind(vcf.init, comb.depth)
    
    write.table(vcf.body, file = paste0("temp.body.", seed), quote = FALSE, sep = "\t", row.names = FALSE, col.names = F) 
    
    system(paste0("cat head.",seed," temp.body.",seed," > temp.",seed,".vcf"))
    poly.test <- VCF2RADdata(paste0("temp.",seed,".vcf"), phaseSNPs = FALSE, 
                             min.ind.with.reads = 0,
                             min.ind.with.minor.allele = 0)
    file.remove(paste0("head.",seed), paste0("temp.body.",seed) ,paste0("temp.",seed,".vcf"))
  } else {
  
    poly.test <- VCF2RADdata(vcf, phaseSNPs = FALSE, 
                           min.ind.with.reads = 0,
                           min.ind.with.minor.allele = 0)
  }
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
  
  # this will change according to the vcf - bug!!
  pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1)
  if(length(unique(pos)) ==1){
    pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",2)
    pos <- paste0(sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1), "_", pos)
  }  else {
  if(tech.issue) # Muda conforme o software de chamada, tem q ver como deixar universal
    pos <- gsub(":", "_", pos)
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
  cat("This approach changed", (1- sum(new.geno == onemap.obj$geno)/length(new.geno))*100,"% of the genotypes")
  
  onemap.obj$geno <- new.geno
  
  # Removing markers
  onemap.obj$n.mar <- length(keep.mks)
  onemap.obj$segr.type <- onemap.obj$segr.type[keep.mks]
  onemap.obj$segr.type.num <- onemap.obj$segr.type.num[keep.mks]
  onemap.obj$CHROM <- onemap.obj$CHROM[keep.mks]
  onemap.obj$POS <- onemap.obj$POS[keep.mks]
  
  polyrad.one <- create_probs(onemap.obj = onemap.obj, genotypes_probs =  genotypes[,3:5])
  
  return(polyrad.one)
}
