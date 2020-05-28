#' Function to simulate depths and to convert PedigreeSim file to vcf file
#' 
#' Given PedigreeSim .dat .map and .chrom files generate vcf file with depths count
#' in AD format field estimated by negative binomial. The function receives genotypes codification in
#' genotypes.dat according with Wu et. al 2002a table I observed bands. The a allele will be considered
#' the reference allele. Null alleles are not supported.
#' 
#' @param inputfile file .dat output from PedigreeSim software
#' @param map.file file .map input in PedigreeSim software
#' @param chrom.file file.chrom input in PedigreeSIm software
#' @param out.file path to vcf output file
#' @param mean.depth mean of the negative binomial distribution to generate depth counts
#' @param disper.par dispersion parameter for negative binomial distribution
#' @param mean.phred Sequencing error parameter
#' @param chr.mb Chromossome size in mega base.
#' @param method Choose negative binomial ("neg.binom"), poisson ("poisson") distributions or updog ("updog") model to simulate counts values
#' @param miss.perc Percentage of missing data
#' @param pos Phisical map position of each marker
#' @param chr Chromosome where the marker is positioned
#' @param phase if TRUE the genotypes in VCF will be phased
#' @param bias The bias parameter for updog model. Pr(a read after selected) / Pr(A read after selected).
#' @param od The overdispersion parameter for updog model. See the Details of the rho variable in betabinom.
#' @param disper.par Dispertion parameter for negative binomial 
#' @param haplo.ref character indicating the reference haplotype of genotypes.dat file
#' @param use.as.alleles if \code{TRUE} uses codification in genotypes dat to define the reference and alternative alleles fields in VCF
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
                       chr.mb= 10, 
                       method = c("updog", "neg.binom"), 
                       mean.phred=20, 
                       bias=1, 
                       od=0.001,
                       disper.par=2,
                       pos=NULL,
                       chr=NULL,
                       phase = FALSE,
                       haplo.ref=NULL,
                       use.as.alleles=FALSE){
  
  # Do the checks here
  data <- read.table(paste(inputfile), stringsAsFactors = FALSE, header = TRUE)
  
  # Infos
  rownames(data) <- data[,1]
  data <- data[,-1]
  n.ind <- dim(data)[2]/2
  n.mk <- dim(data)[1]
  
  # Genotypes matrix 
  idx <- rep(1:(n.ind), each = 2)
  gt_matrix <- gt_ref <- gt_alt <-  het_matrix <- matrix(rep(NA,n.ind*n.mk), nrow = n.mk, ncol = n.ind)
  
  data <- as.matrix(data)
  keep.alleles <- data
  
  if(any(data == "o"))
    stop("Null alleles are not supported.")
  
  if(is.null(haplo.ref)){
    h.ref <- data[,1]
  } else{
    h.ref <- data[,which(colnames(data) == haplo.ref)]
  }
  
  alt <- list()
  for(i in 1:n.mk){
    alt <- levels(factor(unlist(data[i,])))[which(levels(factor(unlist(data[i,]))) != h.ref[i])]
    data[i,][which(data[i,]== h.ref[i])] <- "0"
    for(w in 1:length(alt)){
      data[i,][which(data[i,]== alt[w])] <- w
    }
  }
  
  for(i in 1:(length(idx)/2)){
    gt_matrix[,i] <- paste0(data[,which(idx == i)[1]], "|", data[,which(idx == i)[2]])
    het_matrix[,i] <- data[,which(idx == i)[1]] != data[,which(idx == i)[2]]
    gt_ref[,i] <- data[,which(idx == i)[1]]
    gt_alt[,i] <- data[,which(idx == i)[2]]
  }
  
  if(!phase){
    gt_matrix <- gsub("[|]", "/", gt_matrix)
    gt_matrix[which(gt_matrix == "1/0")] <- "0/1"
    gt_matrix[which(gt_matrix == "2/0")] <- "0/2"
    gt_matrix[which(gt_matrix == "3/0")] <- "0/3"
    gt_matrix[which(gt_matrix == "2/1")] <- "1/2"
    gt_matrix[which(gt_matrix == "3/1")] <- "1/3"
    gt_matrix[which(gt_matrix == "3/2")] <- "2/3"
  }
  
  if(counts==TRUE){
    if(method=="neg.binom" ){
      # Negative binomial to estimate the depths (code adaptaded from Gusmap)
      depth <- prob.mat <- matrix(rep(NA, n.ind*n.mk),nrow=n.mk, ncol=n.ind)
      
      prob.mat[which(gt_matrix == "0/0" | gt_matrix == "0|0")] <- 1
      prob.mat[het_matrix] <- 0.5
      prob.mat[is.na(prob.mat)] <- (10^(-mean.phred/10))
      
      depth[which(!is.na(gt_matrix))] <- rnbinom(sum(!is.na(gt_matrix)),mu=mean.depth,size = disper.par)
      
      # Avoiding missing
      idx <- which(depth==0)
      
      while(length(idx) >0){
        depth[which(depth==0)] <- rnbinom(length(idx),mu=mean.depth,size=disper.par)
        idx <- which(depth==0)
      }
      ref_matrix <- matrix(rbinom(n.ind*n.mk, depth, prob.mat), nrow = n.mk)
      
      
      alt_matrix <- depth-ref_matrix      
      
      info <- paste0("DP=", apply(depth,1,sum))
      
    } else if(method=="updog"){
      
      up_matrix <- size_matrix <- ref_matrix <- matrix(rep(NA, length(gt_matrix)), nrow = dim(gt_matrix)[1])
      
      idx <- which(het_matrix)
      up_matrix[idx] <- 1
      idx <- which(gt_matrix == "0/0" | gt_matrix == "0|0")
      up_matrix[idx] <- 0
      idx <- which(is.na(up_matrix))
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
      
      alt_matrix <- size_matrix-ref_matrix
      
      info <- apply(size_matrix, 1, sum)
    }
    
    # VCF format field
    format <- rep("GT:AD", n.mk)
    
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
    
    # search in initial file the alleles from heterozygotes
    idx <- which(het_matrix >= homo_matrix | het_matrix == Inf)
    check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])]
    
    if(phase){ 
      check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- paste0("0|",sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])], "[|]"), unique))
    } else { 
      check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- paste0("0/",sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])], "/"), unique))
    }
    
    # search in initial file the alleles from homozigotes
    idx <- which(het_matrix < homo_matrix)
    check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])]
    
    if(phase){ 
      allele <- sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])], "[|]"), "[", 1)
      check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- paste0(allele, "|", allele) 
    } else {
      allele <- sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])], "/"), "[", 1)
      check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- paste0(allele, "/", allele) 
    }
    
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
    
    
    ad_matrix <- matrix(NA, nrow = nrow(check_matrix), ncol = ncol(check_matrix))
    for(j in 1:dim(check_matrix)[1]){
      for(i in c(3,2,1)){
        if(any(grepl(i, check_matrix[j,]))){
          ncol <- i + 1
          break
        }
      }
      ad_matrix_temp <- matrix(0, nrow= dim(check_matrix)[2], ncol = ncol)
      for(i in 0:(ncol-1)){
        idx <- which(gt_ref[j,] == i & gt_alt[j,] == i)
        idx.sub <- which(ref_matrix[j,idx] != 0)
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- which(alt_matrix[j,idx] != 0)
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
        
        idx <- which(gt_ref[j,] == i & gt_alt[j,] != i)  
        idx.sub <- gt_ref[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- gt_alt[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
        
        idx <- which(gt_ref[j,] != i & gt_alt[j,] == i)  
        idx.sub <- gt_ref[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- gt_alt[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
      }
      ad_matrix[j,] <- apply(ad_matrix_temp, 1, function(x) paste0(x, collapse = ","))
    }
    
    vcf_format <- matrix(paste0(check_matrix, ":", ad_matrix), nrow = n.mk)
    
  } else {
    vcf_format <- check_matrix <- gt_matrix
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
    pos <- round(pos,0)
  } 
  
  id <- rownames(data)
  
  # Defining alleles in field REF and ALT
  ## REF
  if(use.as.alleles){
    ref <- h.ref
  } else{
    ref <- sample(c("A","T", "C", "G"), n.mk, replace = TRUE)
  }
  alt <- rep(NA, length(ref))
  ## ALT
  done <- vector() # vector to store markers already avaliated
  guide <- 1:dim(check_matrix)[1]
  for(j in c(3,2,1)){
    if(length(done) != 0){
      alleles <- apply(check_matrix[-done,], 1, function(x) any(grepl(j, x)))
    } else {
      alleles <- apply(check_matrix, 1, function(x) any(grepl(j, x)))
    }
    if(sum(alleles) != 0){
      alt_temp <- vector()
      if(length(done) != 0)
        ref_temp <- ref[-done][alleles]
      else ref_temp <- ref[alleles]
      
      for(i in 1:length(alleles)){
        if(use.as.alleles){
          if(length(done) == 0)
            temp <- levels(factor(keep.alleles[which(alleles)[i],]))
          else temp <- levels(factor(keep.alleles[-done,][which(alleles)[i],]))
        } else  temp <- sample(c("A","T", "C", "G"), 4, replace = F)
        alt_temp <- rbind(alt_temp, temp[-which(temp == ref_temp[i])])
      }
      rm.colum <- 3-j
      if(rm.colum != 0 & !use.as.alleles){
        alt_temp <- alt_temp[,-c(1:rm.colum)]
      }
      if(length(done) != 0){ 
        if(is(alt_temp, "matrix"))
          alt[-done][alleles] <- apply(alt_temp, 1, function(x) paste0(x, collapse = ","))
        else alt[-done][alleles] <- alt_temp
        done <- c(done,guide[-done][alleles])
      } else {
        if(is(alt_temp, "matrix"))
          alt[alleles] <- apply(alt_temp, 1, function(x) paste0(x, collapse = ","))
        else alt[alleles] <- alt_temp
        done <- c(done,guide[alleles])
      }
    }
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
