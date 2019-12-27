#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: create_depths_profile.R                                       #
# Contains: create_depths_profile                                     #
#                                                                     #
# Written by Cristiane Taniguti                                       #
# copyright (c) 2019, Cristiane Taniguti                              #
#                                                                     #
# First version: 11/2019                                              #
# License: GNU General Public License version 3                       #
#                                                                     #
#######################################################################

#' Create database and ggplot graphic of allele reads depths
#'
#'
#' @param onemap.obj an object of class \code{onemap}.
#' @param vcf.file an object of class \code{vcfR}.
#' @param parent1 a character specifying the first parent ID
#' @param parent2 a character specifying the second parent ID
#' @param vcf.par the vcf parameter that store the allele depth information. 
#' @param recovering logical. If TRUE, all markers in vcf are considere, if FALSE only those in onemap.obj
#' @param mks a vector of characters specifying the markers names to be considered or NULL to consider all markers
#' @param inds a vector of characters specifying the individual names to be considered or NULL to consider all individuals
#' @param GTfrom the graphic should contain the genotypes from onemap.obj or from the vcf? Specify using "onemap" or "vcf".
#' @param alpha define the transparency of the dots in the graphic
#'
#' @return an rds file and a ggplot graphic.
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
#' @seealso \code{\link[onemap]{onemap_read_vcfR}}
#' @keywords depth alleles 
#'   
#'@import tidyr ggplot2
#'@export
create_depths_profile <- function(onemap.obj = NULL, 
                                  vcfR.object = NULL, 
                                  parent1 = NULL,
                                  parent2 = NULL,
                                  vcf.par = "AD",
                                  recovering=TRUE,
                                  mks = NULL,
                                  inds = NULL,
                                  GTfrom= "onemap",
                                  alpha=1,
                                  rds.file = "data.rds"){
  # do the checks
  depths <- extract_depth(vcfR.object = vcfR.object, onemap.object = onemap.obj, vcf.par, parent1, parent2, recovering = recovering)
  
  # parents depth
  alt <- depths$palt %>% data.frame(mks=depths$mks) %>% gather("ind", "alt", -"mks")
  ref <- depths$pref %>% data.frame(mks=depths$mks) %>% gather("ind", "ref", -"mks")
  parents <- merge(alt,ref)
  
  # parents onemap genotypes
  ## By now, only for biallelic <<<<<<<<
  p1 <- p2 <- vector()
  p1[which(onemap.obj$segr.type == "D1.10")] <- 2
  p1[which(onemap.obj$segr.type == "D2.15")] <- 1
  p1[which(onemap.obj$segr.type == "B3.7")] <- 2
  
  p2[which(onemap.obj$segr.type == "D1.10")] <- 1
  p2[which(onemap.obj$segr.type == "D2.15")] <- 2
  p2[which(onemap.obj$segr.type == "B3.7")] <- 2
  
  p.gt <- data.frame(mks=colnames(onemap.obj$geno), p1, p2)
  colnames(p.gt) <- c("mks", parent1, parent2)
  p.gt <- gather(p.gt, "ind", "gt.onemap", -"mks")
  parents <- merge(parents, p.gt)
  
  # parents vcf genotypes
  id.parents <- which(colnames(vcfR.object@gt[,-1]) %in% c(parent1, parent2))
  gts <- vcfR.object@gt[,-1] %>% strsplit(":") %>% sapply("[",1) %>% matrix(ncol = dim(vcfR.object@gt)[2] -1)
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcfR.object@fix[,1],"_", vcfR.object@fix[,2])
  
  p.gt <- data.frame(mks = MKS, gts[,id.parents])
  colnames(p.gt) <- c("mks", parent1, parent2)
  p.gt <- gather(p.gt, "ind", "gt.vcf", -"mks")
  parents <- merge(parents, p.gt)
  
  # progeny depth
  alt <- depths$oalt %>% data.frame(mks=depths$mks) %>% gather("ind", "alt", -"mks")
  ref <- depths$oref %>% data.frame(mks=depths$mks) %>% gather("ind", "ref", -"mks")
  progeny <- merge(alt,ref)
  
  # progeny onemap genotypes
  gt <- data.frame(ind = rownames(onemap.obj$geno), onemap.obj$geno)
  gt <- gather(gt,  "mks","gt.onemap", -"ind")
  progeny <- merge(progeny, gt)
  
  # progeny vcf genotypes
  pro.gt <- data.frame(mks = MKS, gts[,-id.parents])
  colnames(pro.gt) <- c("mks", colnames(vcfR.object@gt)[-1][-id.parents])
  pro.gt <- gather(pro.gt, "ind", "gt.vcf", -"mks")
  progeny <- merge(progeny, pro.gt)
  data <- rbind(parents, progeny)
  data$gt.onemap[which(data$gt.onemap==0)] <- "missing"
  data$gt.onemap[which(data$gt.onemap==1)] <- "homozygous"
  data$gt.onemap[which(data$gt.onemap==2)] <- "heterozygote"
  data$gt.onemap[which(data$gt.onemap==3)] <- "homozygous"
  
  
  # removing phased
  data$gt.vcf <- gsub(pattern = "[|]", replacement = "/", data$gt.vcf)
  
  data$ind <- as.factor(data$ind)
  data$gt.onemap <- as.factor(data$gt.onemap)
  data$gt.vcf <- as.factor(data$gt.vcf)
  
  if(class(mks) == "character"){
    data <- data[which(data$mks %in% mks),]
  }
  
  if(class(inds) == "character"){
    data <- data[which(data$ind %in% inds),]
  }
  
  if(GTfrom == "onemap"){
    colors <- c("#58355e", "#4D9DE0", "#ADE25D")
    names(colors) <-  c("missing", "homozygous", "heterozygote")
    
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.onemap)) + 
      geom_point(alpha = alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      scale_colour_manual(name="Genotypes", values = colors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
      xlim(0,5*summary(data$ref[-which(data$ref == 0)])[5]) +
      ylim(0,5*summary(data$alt[-which(data$alt == 0)])[5])
  } else if(GTfrom == "vcf"){
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.vcf)) + 
      geom_point(alpha=alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+ 
      xlim(0,5*summary(data$ref[-which(data$ref == 0)])[5]) +
      ylim(0,5*summary(data$alt[-which(data$alt == 0)])[5])
  }
  
  saveRDS(data, file = rds.file)
  
  cat("Summary of reference counts: \n")
  print(summary(data$ref[-which(data$ref == 0)]))
  
  cat("Summary of alternative counts: \n")
  print(summary(data$alt[-which(data$alt == 0)]))
  
  return(p)
  
}
