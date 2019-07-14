## Bug in simulating F2

#' Function to simulate depths and to convert PedigreeSim file to vcf file
#' 
#' Given PedigreeSim .dat .map and .chrom files generate vcf file with depths count
#' in AD format field estimated by negative binomial. Works only for biallelic markers, 
#' the allels can be named as A (reference) and B, or by DNA bases A, T, C and G. In this last case, 
#' a vector indicanting the reference allele is needed.
#' 
#' @param inputfile file .dat output from PedigreeSim software
#' @param map.file file .map input in PedigreeSim software
#' @param chrom.file file.chrom input in PedigreeSIm software
#' @param out.file path to vcf output file
#' @param mean.depth mean of the negative binomial distribution to generate depth counts
#' @param disper.par dispersion parameter for negative binomial distribution
#' @param mean.phred Sequencing error parameter
#' @param chr.mb Chromossome size in mega base.
#' @param method Choose negative binomial ("neg.binom") distribution or updog ("updog") model to simulate counts values
#' @param miss.perc Percentage of missing data
#' @param pos Phisical map position of each marker
#' @param chr Chromosome where the marker is positioned
#' @param haplo.ref The homologue or haplotype which correspond to the reference genome
#' @param phase if TRUE the genotypes in VCF will be phased
#' 
#' @return vcf file located in out.file defined path
#' 
#' @seealso vcf file description 
#' <http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/>
#' 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' 
#' @references 
#' Voorrips, R. E. , Maliepaard, C. A. (2012) The simulation of meiosis in 
#' diploid and tetraploid organisms using various genetic models.
#' \emph{BMC Bioinformatics} 13  - ISSN 1471-2105 - 12 p.
#' 
#' Bilton,T. P. , Schofield,  M. R. , Black,  M. A. , Chagné, D. , Wilcox, P. L. 
#' and Dodds, K. G. (2018) Accounting for Errors in Low Coverage High-Throughput
#' Sequencing Data When Constructing Genetic Maps Using Biparental Outcrossed Populations.
#' \emph{Genetics}. ISSN 0016-6731.
#' @export
pedsim2vcf <- function(inputfile=NULL, 
                       map.file=NULL, 
                       chrom.file=NULL,
                       out.file="out.pedsim2vcf.txt", 
                       miss.perc = 0, 
                       counts=TRUE, 
                       mean.depth=20, 
                       p.mean.depth = 20, 
                       disper.par=2, 
                       chr.mb= 10, 
                       method = c("updog", "neg.binom", "poisson"), 
                       mean.phred=20, 
                       bias=0, 
                       od=0,
                       pos=NULL,
                       haplo.ref=NULL,
                       chr=NULL,
                       phase = FALSE){
  
  # Do the checks here
  
  data <- read.table(paste(inputfile), stringsAsFactors = FALSE, header = TRUE)
  
  # Infos
  rownames(data) <- data[,1]
  data <- data[,-1]
  n.ind <- dim(data)[2]/2
  n.mk <- dim(data)[1]
  
  # Reference haplotype - If there is no haploid reference, A will be the reference and B the alternative
  if(is.null(haplo.ref)){
    h.ref <- rep("a", n.mk)
  } else{
    h.ref <- data[,which(colnames(data) == haplo.ref)]
    alt <- vector()
    for(i in 1:n.mk)
      alt <- c(alt, levels(factor(unlist(data[i,])))[which(levels(factor(unlist(data[i,]))) != h.ref[i])])
  }
  
  # Genotypes matrix --- only for biallelic markers
  idx <- rep(1:(n.ind), each = 2)
  gt_matrix <- matrix(rep(NA,n.ind*n.mk), nrow = n.mk, ncol = n.ind)
  if(phase){
    for(j in 1:n.mk){
      for(i in 1:(length(idx)/2)){
        if(data[j,which(idx==i)[1]] == h.ref[j] & data[j,which(idx==i)[2]] == h.ref[j])
          gt_matrix[j,i] <- "0|0"
        if(data[j,which(idx==i)[1]] == h.ref[j] & data[j,which(idx==i)[2]] != h.ref[j])
          gt_matrix[j,i] <- "0|1"
        if(data[j,which(idx==i)[1]] != h.ref[j] & data[j,which(idx==i)[2]] == h.ref[j])
          gt_matrix[j,i] <- "1|0"
        if(data[j,which(idx==i)[1]] != h.ref[j] & data[j,which(idx==i)[2]] != h.ref[j])
          gt_matrix[j,i] <- "1|1"
      }
    }
  } else {
    for(j in 1:n.mk){
      for(i in 1:(length(idx)/2)){
        if(data[j,which(idx==i)[1]] == h.ref[j] & data[j,which(idx==i)[2]] == h.ref[j])
          gt_matrix[j,i] <- "0/0"
        if(data[j,which(idx==i)[1]] == h.ref[j] & data[j,which(idx==i)[2]] != h.ref[j])
          gt_matrix[j,i] <- "0/1"
        if(data[j,which(idx==i)[1]] != h.ref[j] & data[j,which(idx==i)[2]] == h.ref[j])
          gt_matrix[j,i] <- "0/1"
        if(data[j,which(idx==i)[1]] != h.ref[j] & data[j,which(idx==i)[2]] != h.ref[j])
          gt_matrix[j,i] <- "1/1"
      }
    }
  }
  
  if(counts==TRUE){
    if(method=="neg.binom"){
    # Negative binomial to estimate the depths (code adaptaded from Gusmap)
      depth <- prob.mat <- matrix(rep(NA, n.ind*n.mk),nrow=n.mk, ncol=n.ind)
      
      prob.mat[gt_matrix == "0/0"] <- 1
      prob.mat[gt_matrix == "0/1"] <- 0.5
      prob.mat[gt_matrix == "1/1"] <- (10^(-mean.phred/10))
      
      depth[which(!is.na(gt_matrix))] <- rnbinom(sum(!is.na(gt_matrix)),mu=mean.depth,size=2)
      
      # Avoiding missing
      idx <- which(depth==0)
      
      while(length(idx) >0){
        depth[which(depth==0)] <- rnbinom(length(idx),mu=mean.depth,size=2)
        idx <- which(depth==0)
      }
      
      countsA <- matrix(rbinom(n.ind*n.mk, depth, prob.mat), nrow = n.mk)
      countsB <- depth-countsA

      # AD matrix
      ad_matrix <- matrix(paste0(countsA,",",countsB), nrow = n.mk)
      
      info <- paste0("DP=", apply(depth,1,sum))
      
    } else if(method == "poisson"){
      depth <- prob.mat <- matrix(rep(NA, n.ind*n.mk),nrow=n.mk, ncol=n.ind)
      
      prob.mat[gt_matrix == "0/0"] <- 1
      prob.mat[gt_matrix == "0/1"] <- 0.5
      prob.mat[gt_matrix == "1/1"] <- (10^(-mean.phred/10))
      
      depth[which(!is.na(gt_matrix))] <- rpois(sum(!is.na(gt_matrix)),mean.depth)
      
      # Avoiding missing
      idx <- which(depth==0)
      
      while(length(idx) >0){
        depth[which(depth==0)] <- rpois(length(idx),mean.depth)
        idx <- which(depth==0)
      }
      
      countsA <- matrix(rbinom(n.ind*n.mk, depth, prob.mat), nrow = n.mk)
      countsB <- depth-countsA
      
      # AD matrix
      ad_matrix <- matrix(paste0(countsA,",",countsB), nrow = n.mk)
      
      info <- paste0("DP=", apply(depth,1,sum))
      
      } else if(method=="updog"){
        
        up_matrix <- size_matrix <- ref_matrix <- matrix(rep(NA, length(gt_matrix)), nrow = dim(gt_matrix)[1])
        
        idx <- which(gt_matrix == "0/1")
        up_matrix[idx] <- 1
        idx <- which(gt_matrix == "0/0")
        up_matrix[idx] <- 0
        idx <- which(gt_matrix == "1/1")
        up_matrix[idx] <- 2
        
        size_matrix <- matrix(rnbinom(length(gt_matrix),mu=mean.depth,size=disper.par), nrow = dim(gt_matrix)[1])

        # Parents with other depth
        if(!is.null(p.mean.depth)){
          size_matrix[,1:2] <- rnbinom(length(size_matrix[,1:2]),mu=p.mean.depth,size=disper.par)
        }
        
        mis <- which(size_matrix==0)
        
        while(length(mis) > 0){
          size_matrix[mis] <- rnbinom(length(mis),mu=mean.depth,size=disper.par)
          mis <- which(size_matrix==0)
        }
        
        for(i in 1:dim(up_matrix)[1]){
          ref_matrix[i,] <- updog::rflexdog(sizevec = size_matrix[i,],
                          geno=up_matrix[i,],
                          ploidy = 2,
                          seq=10^(-mean.phred/10),
                          bias=bias,
                          od = od)
        }
        
        ad_matrix <- matrix(paste0(size_matrix-ref_matrix, ",", ref_matrix), nrow = dim(up_matrix)[1])

        info <- apply(size_matrix, 1, sum)
      }
    
    # VCF format field
    format <- rep("GT:AD", n.mk)
    
    # Update geno matrix with counts
    ref_matrix <- matrix(as.numeric(unlist(lapply(strsplit(ad_matrix, split = ","), "[[",1))), nrow = n.mk)
    alt_matrix <- matrix(as.numeric(unlist(lapply(strsplit(ad_matrix, split = ","), "[[",2))), nrow = n.mk)
    
    # Select the number of less frequent allele
    minor_matrix <- matrix(rep(NA, n.mk*n.ind), nrow = n.mk, ncol = n.ind)
    minor_matrix[which(ref_matrix >= alt_matrix)] <- alt_matrix[which(ref_matrix >= alt_matrix)]
    minor_matrix[which(ref_matrix < alt_matrix)] <- ref_matrix[which(ref_matrix < alt_matrix)]
    
    # Probabilities calculation by binomial distribution
    tot_matrix <- ref_matrix + alt_matrix
    tot_matrix[which(tot_matrix==0)] <- NA
    het_matrix <- choose(tot_matrix, minor_matrix)*(0.5^minor_matrix)*(0.5^(tot_matrix-minor_matrix))
    homo_matrix <- choose(tot_matrix, minor_matrix)*((10^(-mean.phred/10))^minor_matrix)*((1-(10^((-mean.phred/10))))^(tot_matrix-minor_matrix))
    homo.oalt <- choose(tot_matrix, minor_matrix)*((10^(-mean.phred/10))^(tot_matrix-minor_matrix))*((1-(10^((-mean.phred/10))))^minor_matrix)
    
    # Reviewed matrix 
    check_matrix <- matrix(rep(NA, n.mk*n.ind), nrow = n.mk, ncol = n.ind)
    check_matrix[which(het_matrix >= homo_matrix)] <- "0/1"
    check_matrix[which(het_matrix == Inf)] <- "0/1"
    idx <- which(het_matrix < homo_matrix) 
    check_matrix[idx][ref_matrix[idx] > alt_matrix[idx]] <- "0/0"
    check_matrix[idx][ref_matrix[idx] < alt_matrix[idx]] <- "1/1"
    check_matrix[is.na(check_matrix)] <- "./."
    
    chang <- table((gt_matrix == check_matrix))
    
    if(dim(chang) ==2){
      cat("Counts simulation changed", (chang[1]/(chang[1] + chang[2]))*100,
          "% of the given genotypes\n")
    } else if(names(chang) == "FALSE"){
      cat("All genotypes were changed")
    } else{
      cat("None genotypes were changed")
    }
    
    # Reference allele is the most frequent
    # for(i in 1:dim(check_matrix)[[1]]){
    #   z <- table(check_matrix[i,])
    #   w <- strsplit(ad_matrix[i,], split = ",")
    #   if(length(z) == 2){
    #     if(all(names(z) == c("0/1", "1/1"))){
    #       idx1 <- which(check_matrix[i,] == "1/1")
    #       idx2 <- which(check_matrix[i,] == "0/0")
    #       check_matrix[i,][idx1] <- "0/0"
    #       w[idx1] <- lapply(w[idx1], "rev")
    #       check_matrix[i,][idx2] <- "1/1"
    #       w[idx2] <- lapply(w[idx2], "rev")
    #       ad_matrix[i,] <- unlist(lapply(w, function(x) paste0(x, collapse = ",")))
    #     }
    #   } else if(length(z) == 3){
    #     if(all(names(z) == c("0/0", "0/1", "1/1"))){
    #       if(z[[1]] < z[[3]]){
    #         idx1 <- which(check_matrix[i,] == "1/1")
    #         idx2 <- which(check_matrix[i,] == "0/0")
    #         check_matrix[i,][idx1] <- "0/0"
    #         w[idx1] <- lapply(w[idx1], "rev")
    #         check_matrix[i,][idx2] <- "1/1"
    #         w[idx2] <- lapply(w[idx2], "rev")
    #         ad_matrix[i,] <- unlist(lapply(w, function(x) paste0(x, collapse = ",")))
    #       }
    #     }
    #   }
    # }
    
    vcf_format <- matrix(paste0(check_matrix, ":", ad_matrix), nrow = n.mk)
    
    } else {
        vcf_format <- gt_matrix
        format <- rep("GT", n.mk)
        info <- "."
    }
    
   # Adding missing data
    miss <- sample(1:length(vcf_format), length(vcf_format)*(miss.perc/100))
    if(length(miss)>0){
      if(counts){
        vcf_format[miss] <- "./.:0,0"  
      } else{ vcf_format[miss] <- "./." }
    }
      
      
    names1 <- lapply(strsplit(colnames(data)[-1], "_"), function (x) x[-length(x)])
    names1 <- unique(sapply(names1, function(x) if(length(x)>1) paste0(x[1], "_", x[2]) else x)) 
    
    colnames(vcf_format) <- names1
  
    if(is.null(chr)){
    chr.info <- read.table(map.file, header = TRUE, stringsAsFactors = FALSE)
    chr <- chr.info$chromosome
    }
    
    if(is.null(pos)){
    pos.info <- read.table(chrom.file, header = TRUE, stringsAsFactors = FALSE)
    pos <- chr.info$position*((chr.mb*1000)/mean(pos.info$length))
    } 
    
    id <- rownames(data)
    if(is.null(haplo.ref)){
      ref <- sample(c("A","T", "C", "G"), n.mk, replace = TRUE)
      alt <- sample(c("A","T", "C", "G"), n.mk, replace = TRUE)
      
      while(!all(!ref==alt)){
        alt[which(ref==alt)] <- sample(c("A","T", "C", "G"), length(which(ref==alt)), replace = TRUE)
      }
    } else{
      ref <- h.ref 
      alt <- alt
    }
    
    qual <- filter <- rep(".", n.mk)
    # transformar em vetor
    vcf_file_mks <- data.frame("CHROM"=chr, "POS"=pos, "ID"= id,"REF"=ref, "ALT"=alt, 
                               "QUAL"=qual, "FILTER"=filter,"INFO"=info,"FORMAT"=format,vcf_format, stringsAsFactors = FALSE)
  
    vcf_vector <- apply(vcf_file_mks, 1, function(x) paste(x, collapse = "\t"))
  
    header1 <- paste0(colnames(vcf_file_mks), collapse = "\t")
    header <- paste0("##fileformat=VCFv4.1", "\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", '\n',
                   "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">",'\n',
                   "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", "\n", "#", header1)
  
    vcf <- c(header,vcf_vector)
    write.table(vcf, file = paste(out.file), quote = FALSE, row.names = FALSE,  col.names = FALSE)
}
