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
#' @param vcf path to VCF file.
#' @param vcfR.object object of class vcfR;
#' @param parent1 a character specifying the first parent ID
#' @param parent2 a character specifying the second parent ID
#' @param vcf.par the vcf parameter that store the allele depth information. 
#' @param recovering logical. If TRUE, all markers in vcf are consider, if FALSE only those in onemap.obj
#' @param mks a vector of characters specifying the markers names to be considered or NULL to consider all markers
#' @param inds a vector of characters specifying the individual names to be considered or NULL to consider all individuals
#' @param GTfrom the graphic should contain the genotypes from onemap.obj or from the vcf? Specify using "onemap", "vcf" or "prob".
#' @param alpha define the transparency of the dots in the graphic
#' @param rds.file rds file name to store the data frame with values used to build the graphic
#' @param x_lim set scale limit for x axis
#' @param y_lim set scale limit for y axis
#' @param verbose If \code{TRUE}, print tracing information.
#' 
#' @return an rds file and a ggplot graphic.
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#' 
#' @seealso \code{\link[onemap]{onemap_read_vcfR}}
#' 
#' @keywords depth alleles 
#'
#'
#'@import tidyr ggplot2
#'@importFrom vcfR read.vcfR extract.gt masplit
#'
#'@export
create_depths_profile <- function(onemap.obj = NULL, 
                                  vcfR.object = NULL,
                                  vcf = NULL, 
                                  parent1 = NULL,
                                  parent2 = NULL,
                                  vcf.par = "AD",
                                  recovering=FALSE,
                                  mks = NULL,
                                  inds = NULL,
                                  GTfrom= "onemap",
                                  alpha=1,
                                  rds.file = "data.rds",
                                  y_lim = NULL,
                                  x_lim = NULL, verbose=TRUE){
  
  if (is.null(vcf) & is.null(vcfR.object)) {
    stop("You must specify one vcf file.")
  }
  
  if(is.null(vcfR.object)){
    vcfR.object <- read.vcfR(vcf, verbose = F)
  } else vcfR.object <- vcfR.object
  
  # Exclude multiallelic markers
  if(inherits(onemap.obj, "outcross")){
    idx.mks <- colnames(onemap.obj$geno)[which(!(onemap.obj$segr.type %in% c("B3.7", "D1.10", "D2.15")))]
    if(length(idx.mks) > 0){
      idx.mks <- colnames(onemap.obj$geno)[which((onemap.obj$segr.type %in% c("B3.7", "D1.10", "D2.15")))]
      warning("Only biallelic codominant markers are supported. The multiallelic markers present in onemap object will not be plotted.\n") 
      onemap.obj <- split_onemap(onemap.obj, mks= idx.mks)
    }
  } else if(inherits(onemap.obj,"f2")){
    idx.mks <- colnames(onemap.obj$geno)[which(!(onemap.obj$segr.type %in% c("A.H.B")))]
    if(length(idx.mks) > 0){
      idx.mks <- colnames(onemap.obj$geno)[which((onemap.obj$segr.type %in% c("A.H.B")))]
      warning("Only codominant markers are supported. The dominant markers present in onemap object will not be plotted.\n") 
      onemap.obj <- split_onemap(onemap.obj, mks= idx.mks)
    }
  } 
  
  if(is.null(parent1) | is.null(parent2)) stop("Parents ID must be defined.")
  
  depths <- extract_depth(vcfR.object = vcfR.object, onemap.object = onemap.obj, vcf.par, parent1, parent2,recovering = recovering)
  
  # parents onemap genotypes
  ## Only for biallelic codominant markers
  p1 <- p2 <- vector()
  # parents depth
  alt <- depths$palt %>% cbind(mks=depths$mks) %>% as.data.frame %>% gather("ind", "alt", -"mks")
  ref <- depths$pref %>% cbind(mks=depths$mks) %>% as.data.frame %>% gather("ind", "ref", -"mks")
  parents <- merge(alt,ref)
  parents$mks <- gsub("[|]", ".", parents$mks)
  
  if(inherits(onemap.obj, "outcross")){
    p1[which(onemap.obj$segr.type == "D1.10")] <- 2
    p1[which(onemap.obj$segr.type == "D2.15")] <- 1
    p1[which(onemap.obj$segr.type == "B3.7")] <- 2
    
    p2[which(onemap.obj$segr.type == "D1.10")] <- 1
    p2[which(onemap.obj$segr.type == "D2.15")] <- 2
    p2[which(onemap.obj$segr.type == "B3.7")] <- 2
    
  } else  if(inherits(onemap.obj,c("riself", "risib", "f2"))){
    p1 <- 1
    p2 <- 3
  } else{
    p1 <- 1
    p2 <- 2
  }
  
  id.parents <- c(parent1, parent2)
  p.gt <- data.frame(mks=colnames(onemap.obj$geno), p1, p2)
  
  colnames(p.gt) <- c("mks", id.parents)
  p.gt <- gather(p.gt, "ind", "gt.onemap", -"mks")
  if(all(is.na(match(parents$mks, p.gt$mks))))
    parents$mks <- paste0(vcfR.object@fix[,1], "_", vcfR.object@fix[,2])[match(parents$mks, vcfR.object@fix[,3])]
  parents <- merge(parents, p.gt)
  
  # parents vcf genotypes
  idx.parents <- which(colnames(vcfR.object@gt[,-1]) %in% id.parents)
  gts <- extract.gt(vcfR.object)
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcfR.object@fix[,1], "_",vcfR.object@fix[,2])
  
  p.gt <- data.frame(mks = MKS, gts[,idx.parents], stringsAsFactors = F)
  colnames(p.gt) <- c("mks", colnames(gts)[idx.parents])
  p.gt <- gather(p.gt, "ind", "gt.vcf", -"mks")
  if(all(is.na(match(parents$mks, p.gt$mks))))
    parents$mks <- vcfR.object@fix[,3][match(parents$mks, paste0(vcfR.object@fix[,1], "_", vcfR.object@fix[,2]))]
  parents <- merge(parents, p.gt)
  
  if(inherits(onemap.obj, c("outcross", "f2"))){
    parents <- data.frame(parents, A=NA, AB=NA, BA=NA, B=NA)
  } else {
    parents <- data.frame(parents, A=NA, AB=NA)
  }
  # progeny depth
  alt <- depths$oalt %>% cbind(mks=depths$mks) %>% as.data.frame %>% gather("ind", "alt", -"mks")
  ref <- depths$oref %>% cbind(mks=depths$mks) %>% as.data.frame %>% gather("ind", "ref", -"mks")
  progeny <- merge(alt,ref)
  
  # progeny onemap genotypes
  gt <- as.data.frame(cbind(ind = rownames(onemap.obj$geno), onemap.obj$geno))
  gt <- gather(gt,  "mks","gt.onemap", -"ind")
  temp <- match(paste0(gt$mks, "_", gt$ind), rownames(onemap.obj$error))
  gt <- data.frame(gt, onemap.obj$error[temp,])
  
  if(inherits(onemap.obj, c("outcross", "f2"))){
    colnames(gt) <- c("ind", "mks", "gt.onemap", "A", "AB", "BA", "B")
  } else {
    colnames(gt) <- c("ind", "mks", "gt.onemap", "A", "AB")
  }
  progeny <- merge(progeny, gt)
  
  # progeny vcf genotypes
  pro.gt <- data.frame(mks = MKS, gts[, -idx.parents], stringsAsFactors = F)
  colnames(pro.gt) <- c("mks", colnames(gts)[-idx.parents])
  pro.gt <- gather(pro.gt, "ind", "gt.vcf", -"mks")
  progeny <- merge(progeny, pro.gt)
  
  data <- rbind(parents, progeny)
  
  # Add marker type
  data$mk.type <- onemap.obj$segr.type[match(data$mks, colnames(onemap.obj$geno))]
  
  data$gt.onemap.alt.ref <- data$gt.onemap
  data$gt.onemap.alt.ref[which(data$gt.onemap == 0)] <- "missing"
  data$gt.onemap.alt.ref[which(data$gt.onemap == 2)] <- "heterozygous"
  
  # Search the ref and alt alleles using parents
  data$alt <- as.numeric(data$alt)
  data$ref <- as.numeric(data$ref)
  
  # D2.15 = aa x ab
  idx <- data$gt.onemap == 1 & (data$mk.type == "D2.15")
  temp <- data[idx,]
  temp.ref <- temp[temp$ind==id.parents[1] & temp$alt < temp$ref,]$mks # if parent 1 have reference allele
  data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-ref"
  temp.ref <- temp[temp$ind==id.parents[1] & temp$alt > temp$ref,]$mks # if parent 1 have alternative allele
  data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-alt"
  # We do not expect heterozygous parents in these cases, if counts are the same, NA is inserted
  temp.ref <- temp[temp$ind==id.parents[1] & temp$alt == temp$ref,]$mks 
  if(length(temp.ref) > 0){
    data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-alt == ref"
  }
  
  # D1.10 = ab x aa
  idx <- data$gt.onemap == 1 & (data$mk.type == "D1.10")
  temp <- data[idx,]
  temp.ref <- temp[temp$ind==id.parents[2] & temp$alt < temp$ref,]$mks # if parent 2 have reference allele
  data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-ref"
  temp.ref <- temp[temp$ind==id.parents[2] & temp$alt > temp$ref,]$mks # if parent 2 have alternative allele
  data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-alt"
  # We do not expect heterozygous parents in these cases, if counts are the same, NA is inserted
  temp.ref <- temp[temp$ind==id.parents[2] & temp$alt == temp$ref,]$mks 
  if(length(temp.ref) > 0){
    data$gt.onemap.alt.ref[which(idx & data$mks %in% temp.ref)] <- "homozygous-alt == ref"
  }
  
  # B3.7 - can not take the alleles from parents, here we make by individual
  idx <- data$gt.onemap %in% c(1,3) & (data$mk.type == "B3.7" | data$mk.type == "A.H.B")
  data$gt.onemap.alt.ref[idx & data$alt > data$ref] <- "homozygous-alt"
  data$gt.onemap.alt.ref[idx & data$alt < data$ref] <- "homozygous-ref"
  data$gt.onemap.alt.ref[idx & data$alt == data$ref] <- "homozygous-alt == ref" # If counts are the same we can not recover the information
  
  # removing phased
  data$gt.vcf <- gsub(pattern = "[|]", replacement = "/", data$gt.vcf)
  
  if(length(grep(pattern = ",", data$gt.vcf)) > 0){ # If the vcf do not have GT field
    data$gt.vcf[grep(pattern = ",", data$gt.vcf)] <- NA
  }
  
  data$gt.vcf <- gsub(pattern = "[|]", replacement = "/", data$gt.vcf)
  data$gt.vcf.alt.ref <- data$gt.vcf
  data$gt.vcf.alt.ref[grep("[.]", data$gt.vcf.alt.ref)] <- "missing"
  data$gt.vcf.alt.ref[grep("0/0", data$gt.vcf.alt.ref)] <- "homozygous-ref"
  data$gt.vcf.alt.ref[grep("1/1", data$gt.vcf.alt.ref)] <- "homozygous-alt"
  data$gt.vcf.alt.ref[grepl("0/1", data$gt.vcf.alt.ref) | grepl("1/0", data$gt.vcf.alt.ref)] <- "heterozygous"
  
  data$ind <- as.factor(data$ind)
  data$gt.onemap.alt.ref <- as.factor(data$gt.onemap.alt.ref)
  data$gt.vcf.alt.ref <- as.factor(data$gt.vcf.alt.ref)
  
  # Remove multiallelic vcf
  rm.mks <- grep("[2-9]", data$gt.vcf)
  if(length(rm.mks) > 0)
    data <- data[-rm.mks,]
  
  if(inherits(mks,"character")){
    data <- data[which(data$mks %in% mks),]
  }
  
  if(inherits(inds, "character")){
    data <- data[which(data$ind %in% inds),]
  }
  
  data$ref <- as.numeric(data$ref)
  data$alt <- as.numeric(data$alt)
  
  if(is.null(y_lim))
    y_lim <- max(data$ref) 
  if(is.null(x_lim))
    x_lim <- max(data$alt)
  
  if(GTfrom == "onemap"){
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.onemap.alt.ref)) + 
      geom_point(alpha = alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      scale_colour_viridis_d() +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
      xlim(0, x_lim) +
      ylim(0, y_lim)
  } else if(GTfrom == "vcf"){
    p <- data %>% ggplot(aes(x=ref, y=alt, color=gt.vcf.alt.ref)) + 
      geom_point(alpha=alpha) +
      labs(title= "Depths",x="ref", y = "alt", color="Genotypes") +
      scale_colour_viridis_d() +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+ 
      xlim(0, x_lim) +
      ylim(0, y_lim)
  } else if(GTfrom == "prob"){
    if(inherits(onemap.obj, c("outcross", "f2"))){
      idx <- 7:10
    } else {
      idx <- 7:8
    }
    errors <- apply(data[,idx], 1, function(x) {
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
  
  if(verbose){
    cat("Summary of reference counts: \n")
    print(summary(data$ref[-which(data$ref == 0)]))
    
    cat("Summary of alternative counts: \n")
    print(summary(data$alt[-which(data$alt == 0)]))
  }
  return(p)
  
}
