#' Extract allele counts of progeny and parents of vcf file
#' 
#' Uses vcfR package and onemap object to generates list of vectors with
#' reference alelle count and total counts for each marker and genotypes 
#' included in onemap object
#' 
#' @param vcfR.object object output from vcfR package
#' @param onemap.object onemap object output from read_onemap, read_mapmaker or onemap_read_vcf function
#' @param vcf.par vcf format field that contain alelle counts informations, usually AD and DPR
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR objetc
#' @param recovering TRUE/FALSE, if TRUE avaliate all markers from vcf file, if FALSE avaliate only markers in onemap object
#' @return list containing the following components: \item{palt}{a \code{matrix} with parent 1 and 2 
#' alternative alelle counts.} \item{pref}{a \code{matrix} with parent 1 and 2 
#' reference alelle counts.} \item{psize}{a \code{matrix} with parent 1 and 2 
#' total alelle counts.}\item{oalt}{a \code{matrix} with progeny 
#' alternative alelle counts.}\item{oref}{a \code{matrix} with progeny 
#' reference alelle counts.}\item{osize}{a \code{matrix} with progeny 
#' total alelle counts.}\item{n.mks}{total number of markers.} 
#' \item{n.ind}{total number of individuals in progeny.} \item{inds}{progeny individuals identification.}
#' \item{mks}{markers identification.} \item{onemap.object}{same onemap.object inputed}
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' @export

extract_depth <- function(vcfR.object=NULL,
                         onemap.object= NULL,
                         vcf.par = c("GQ","AD", "DPR"),
                         parent1="P1",
                         parent2="P2",
                         recovering = FALSE){
  if(is.null(vcfR.object)) 
   stop("You must specify one vcfR object.")

  if(class(vcfR.object)!="vcfR")
    stop("You must specify one vcfR object.")

  if(is.null(onemap.object)) 
    stop("You must specify one onemap object.")

  if(class(onemap.object)[1]!="onemap")
    stop("You must specify one onemap object.")
  
  # if(deparse(substitute(vcfR.object)) != onemap.object$input)
  #   stop("The onemap object declared is not compatible with the vcfR object")
  
  # Infos
  ind <- rownames(onemap.object$geno)
  IND <- colnames(vcfR.object@gt)[-1]
  mks <- colnames(onemap.object$geno)
  MKS <- vcfR.object@fix[,3]
  n.mks <- length(mks)
  n.ind <- length(ind)
  N.MKs <- dim(vcfR.object@gt)[1]
  N.IND <- dim(vcfR.object@gt)[2]-1
  
  # If there are no marker names
  if(all(is.na(MKS)))
    MKS <- paste0("M", 1:length(MKS))

  parents <- c(which(IND == parent1),which(IND == parent2))
  if(recovering==FALSE){
  rm.mks <- which(MKS %in% mks == FALSE)
  rm.ind <- which(IND[-parents] %in% ind == FALSE)
  CHROM <- onemap.object$CHROM
  POS <- onemap.object$POS
  } else {
    CHROM <- vcfR.object@fix[,1]
    POS <- vcfR.object@fix[,2]
    rm.mks <- NULL 
    rm.ind <- NULL
  }
  
  n.par <- length(strsplit(vcfR.object@gt[1,1], split=":")[[1]])
  
  if(vcf.par=="GQ")
    par <- which(strsplit(vcfR.object@gt[1,1], split=":")[[1]]=="GQ")
  if(vcf.par=="AD")
    par <- which(strsplit(vcfR.object@gt[1,1], split=":")[[1]]=="AD")
  if(vcf.par=="DPR")
    par <-  which(strsplit(vcfR.object@gt[1,1], split=":")[[1]]=="DPR")
  if(length(par)==0)
    stop("There is no ", vcf.par, " field in the vcfR.object file. Error probabilities can't be generated")
  
  # Spliting fields
  split.gt <- strsplit(vcfR.object@gt, split=":")
  
  # Checking format of missing data
  miss <- unlist(lapply(split.gt,  length))
  miss.num <- which(miss < par)
  
  # Replacing missing data with compatible format
  if(vcf.par=="AD" | vcf.par=="DPR"){
    vcfR.object@gt[miss.num] <- paste0(rep("0,0",par),":", collapse = "")
    split.gt <- strsplit(vcfR.object@gt, split=":")
  } else {vcfR.object@gt[miss.num] <- paste0(rep("0",par),":", collapse = "")}
  
  # Extracting choosed parameter matrix
  par_matrix <- matrix(unlist(lapply(split.gt,  "[[", par)), nrow = N.MKs, ncol = N.IND+1)[,-1]
  
  # Replacing missing data with compatible format
  if(length(which(par_matrix == ".")) > 0)
    par_matrix[which(par_matrix == ".")] <- "0,0"
  
  if(length(rm.mks)>0 & length(rm.ind)>0){
      par_matrix <- par_matrix[-rm.mks, -rm.ind]
      IND <- IND[-rm.ind]
      MKS <- MKS[-rm.mks]
  } else if(length(rm.mks)>0){
      par_matrix <- par_matrix[-rm.mks,]
      MKS <- MKS[-rm.mks]
  } else if(length(rm.ind)>0){
      par_matrix <- par_matrix[,-rm.ind]
      IND <- IND[-rm.ind]
  }

    n.ind <- N.IND - length(rm.ind)
    n.mks <- N.MKs - length(rm.mks)
    
  # The probabilities must be calculated if AD or DPR parameters were choosed
  if(vcf.par=="AD" | vcf.par=="DPR"){
    ref_matrix <- matrix(as.numeric(unlist(lapply(strsplit(par_matrix, split = ","), "[[",1))), nrow = n.mks, ncol = n.ind)
    alt_matrix <- matrix(as.numeric(unlist(lapply(strsplit(par_matrix, split = ","), "[[",2))), nrow = n.mks, ncol = n.ind)
  }

    colnames(alt_matrix) <- colnames(ref_matrix) <- IND
    rownames(ref_matrix) <- rownames(alt_matrix) <- MKS
    
    if(vcf.par=="DPR"){
      size_matrix <- ref_matrix
      ref_matrix <- size_matrix - alt_matrix 
      alt_matrix <- size_matrix - ref_matrix
    } else {
    size_matrix <- ref_matrix + alt_matrix
    }
    
    # Separing offspring and parents
    idx <- c(which(IND==parent1),which(IND==parent2))
    palt <- alt_matrix[,idx]
    pref <- ref_matrix[,idx]
    psize <- size_matrix[,idx]
    oalt <- alt_matrix[,-idx]
    oref <- ref_matrix[,-idx]
    osize <- size_matrix[,-idx]
    IND <- IND[-idx]

    n.ind <- dim(oref)[2]
    n.mks <- dim(oref)[1]
    
    structure(list(palt=palt,
                   pref=pref,
                   psize=psize,
                   oalt=oalt,
                   oref=oref,
                   osize=osize,
                   n.mks=n.mks,
                   n.ind=n.ind,
                   inds = IND,
                   mks = MKS,
                   CHROM = CHROM,
                   POS = POS,
                   onemap.object = onemap.object))   
}
