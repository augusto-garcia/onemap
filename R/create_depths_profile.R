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

globalVariables(c("gt.onemap", "gt.vcf"))

#' Create database and ggplot graphic of allele reads depths
#'
#'
#' @param onemap.obj an object of class \code{onemap}.
#' @param vcfR.object an object of class \code{vcfR}.
#' @param parent1 a character specifying the first parent ID
#' @param parent2 a character specifying the second parent ID
#' @param f1 if your cross type is f2, you must define the F1 individual
#' @param vcf.par the vcf parameter that store the allele depth information. 
#' @param recovering logical. If TRUE, all markers in vcf are considere, if FALSE only those in onemap.obj
#' @param mks a vector of characters specifying the markers names to be considered or NULL to consider all markers
#' @param inds a vector of characters specifying the individual names to be considered or NULL to consider all individuals
#' @param GTfrom the graphic should contain the genotypes from onemap.obj or from the vcf? Specify using "onemap", "vcf" or "prob".
#' @param alpha define the transparency of the dots in the graphic
#' @param rds.file rds file name to store the data frame with values used to build the graphic
#' @param x_lim set scale limit for x axis
#' @param y_lim set scale limit for y axis
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
                                  f1 = NULL,
                                  vcf.par = "AD",
                                  recovering=FALSE,
                                  mks = NULL,
                                  inds = NULL,
                                  GTfrom= "onemap",
                                  alpha=1,
                                  rds.file = "data.rds",
                                  y_lim = NULL,
                                  x_lim = NULL){
  
  # Checks
  if(is(onemap.obj, c("f2 intercross", "f2 backcross")) & is.null(f1)) 
    stop("You must define f1 argument for this cross type \n")
  
  # Exclude multiallelic markers
  if(is(onemap.obj, "outcross")){
    idx.mks <- colnames(onemap.obj$geno)[which(!(onemap.obj$segr.type %in% c("B3.7", "D1.10", "D2.15")))]
    if(length(idx.mks) > 0){
      warning("Only biallelic codominant markers are supported. The multiallelic markers present in onemap object will not be plotted.\n") 
      onemaps <- split_onemap(onemap.obj, mks= idx.mks)
      onemap.obj <- onemaps[[1]]
    }
  } else if(is(onemap.obj,"f2")){
    idx.mks <- colnames(onemap.obj$geno)[which(!(onemap.obj$segr.type %in% c("A.H.B")))]
    if(length(idx.mks) > 0){
      warning("Only codominant markers are supported. The dominant markers present in onemap object will not be plotted.\n") 
      onemaps <- split_onemap(onemap.obj, mks= idx.mks)
      onemap.obj <- onemaps[[1]]
    }
  } else{
    stop("By now, this function is only available for outcrossing and f2 intercross populations\n")
  }
  
  if(is.null(parent1) | is.null(parent2)) stop("Parents ID must be defined.")
  
  # do the checks
  depths <- extract_depth(vcfR.object = vcfR.object, onemap.object = onemap.obj, vcf.par, parent1, parent2, f1= f1,recovering = recovering)
  
  # parents onemap genotypes
  ## Only for biallelic codominant markers
  p1 <- p2 <- vector()
  if(is(onemap.obj, "outcross")){
    # parents depth
    alt <- depths$palt %>% data.frame(mks=depths$mks) %>% gather("ind", "alt", -"mks")
    ref <- depths$pref %>% data.frame(mks=depths$mks) %>% gather("ind", "ref", -"mks")
    parents <- merge(alt,ref)
    
    p1[which(onemap.obj$segr.type == "D1.10")] <- 2
    p1[which(onemap.obj$segr.type == "D2.15")] <- 1
    p1[which(onemap.obj$segr.type == "B3.7")] <- 2
    
    p2[which(onemap.obj$segr.type == "D1.10")] <- 1
    p2[which(onemap.obj$segr.type == "D2.15")] <- 2
    p2[which(onemap.obj$segr.type == "B3.7")] <- 2
    id.parents <- c(parent1, parent2)
    p.gt <- data.frame(mks=colnames(onemap.obj$geno), p1, p2)
  } else if(is(onemap.obj,"f2")){
    # parents depth
    alt <- depths$palt %>% data.frame(mks=depths$mks)
    alt <- cbind(alt, f1)
    colnames(alt) <- c("alt", "mks", "ind")
    ref <- depths$pref %>% data.frame(mks=depths$mks)
    ref <- cbind(ref, f1)
    colnames(ref) <- c("ref", "mks", "ind")
    parents <- merge(alt,ref)
    p1 <- 2
    id.parents <- f1
    p.gt <- data.frame(mks=colnames(onemap.obj$geno), p1)
  } else{
    stop("By now, this function don't support this cross type\n")
  }
  
  colnames(p.gt) <- c("mks", id.parents)
  p.gt <- gather(p.gt, "ind", "gt.onemap", -"mks")
  parents <- merge(parents, p.gt)
  
  # parents vcf genotypes
  idx.parents <- which(colnames(vcfR.object@gt[,-1]) %in% id.parents)
  gts <- vcfR.object@gt[,-1] %>% strsplit(":") %>% sapply("[",1) %>% matrix(ncol = dim(vcfR.object@gt)[2] -1)
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcfR.object@fix[,1],"_", vcfR.object@fix[,2])
  
  p.gt <- data.frame(mks = MKS, gts[,idx.parents], stringsAsFactors = F)
  colnames(p.gt) <- c("mks", id.parents)
  p.gt <- gather(p.gt, "ind", "gt.vcf", -"mks")
  parents <- merge(parents, p.gt)
  
  parents <- data.frame(parents, A=NA, AB=NA, BA=NA, B=NA)
  
  # progeny depth
  alt <- depths$oalt %>% data.frame(mks=depths$mks) %>% gather("ind", "alt", -"mks")
  ref <- depths$oref %>% data.frame(mks=depths$mks) %>% gather("ind", "ref", -"mks")
  progeny <- merge(alt,ref)
  
  # progeny onemap genotypes
  gt <- data.frame(ind = rownames(onemap.obj$geno), onemap.obj$geno)
  gt <- gather(gt,  "mks","gt.onemap", -"ind")
  temp <- match(paste0(gt$mks, "_", gt$ind), rownames(onemap.obj$error))
  gt <- data.frame(gt, onemap.obj$error[temp,])
  colnames(gt) <- c("ind", "mks", "gt.onemap", "A", "AB", "BA", "B")
  progeny <- merge(progeny, gt)
  
  # progeny vcf genotypes
  pro.gt <- data.frame(mks = MKS, gts[,-idx.parents], stringsAsFactors = F)
  colnames(pro.gt) <- c("mks", colnames(vcfR.object@gt)[-1][-idx.parents])
  pro.gt <- gather(pro.gt, "ind", "gt.vcf", -"mks")
  progeny <- merge(progeny, pro.gt)
  data <- rbind(parents, progeny)
  data$gt.onemap[which(data$gt.onemap==0)] <- "missing"
  data$gt.onemap[which(data$gt.onemap==1)] <- "homozygous"
  data$gt.onemap[which(data$gt.onemap==2)] <- "heterozygote"
  data$gt.onemap[which(data$gt.onemap==3)] <- "homozygous"
  
  # removing phased
  data$gt.vcf <- gsub(pattern = "[|]", replacement = "/", data$gt.vcf)
  
  if(length(grep(pattern = ",", data$gt.vcf)) > 0){ # If the vcf do not have GT field
    data$gt.vcf[grep(pattern = ",", data$gt.vcf)] <- NA
  }
  
  data$ind <- as.factor(data$ind)
  data$gt.onemap <- as.factor(data$gt.onemap)
  data$gt.vcf <- as.factor(data$gt.vcf)
  
  if(class(mks) == "character"){
    data <- data[which(data$mks %in% mks),]
  }
  
  if(class(inds) == "character"){
    data <- data[which(data$ind %in% inds),]
  }
  
  if(is.null(y_lim))
    y_lim <- max(data$ref) 
  if(is.null(x_lim))
    x_lim <- max(data$alt)
  
  if(GTfrom == "onemap"){
    colors <- c("#58355e", "#4D9DE0", "#ADE25D")
    names(colors) <-  c("missing", "homozygous", "heterozygote")
    
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.onemap)) + 
      geom_point(alpha = alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      scale_colour_manual(name="Genotypes", values = colors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
      xlim(0, x_lim) +
      ylim(0, y_lim)
  } else if(GTfrom == "vcf"){
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.vcf)) + 
      geom_point(alpha=alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+ 
      xlim(0, x_lim) +
      ylim(0, y_lim)
  } else if(GTfrom == "prob"){
    errors <- apply(data[,7:10], 1, function(x) {
      if(all(is.na(x)) | all(x == 1)) {
        return(NA) 
      } else { 
        z <- 1 - x[which.max(x)]
        return(z)
      }
    })
    p <- data %>% ggplot(aes(x=ref, y=alt, color=errors)) + 
      geom_point(alpha=alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      scale_colour_gradient(low = "#70ED57", high = "#F62A2C")+
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+ 
      xlim(0, x_lim) +
      ylim(0, y_lim)
  }
  
  saveRDS(data, file = rds.file)
  
  cat("Summary of reference counts: \n")
  print(summary(data$ref[-which(data$ref == 0)]))
  
  cat("Summary of alternative counts: \n")
  print(summary(data$alt[-which(data$alt == 0)]))
  
  return(p)
  
}
