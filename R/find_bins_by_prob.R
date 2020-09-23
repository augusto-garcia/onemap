##' Not ready yet - use carefully 
##' 
##' @export
find_bins_by_probs <- function(onemap_obj, threshold.probs = 0.0001, threshold.count = 0.08){
  mks.types <- onemap_obj$segr.type.num
  MKS <- colnames(onemap_obj$geno)
  n.mks <- length(MKS)
  n.ind <- onemap_obj$n.ind
  types <- unique(mks.types)
  threshold.num <- n.ind*4*threshold.count
  
  mars <- rep(colnames(onemap_obj$geno), onemap_obj$n.ind)
  temp <- cbind(mars, onemap_obj$error)
  error_list <- split(onemap_obj$error, mars)
  
  bins.all.types <- list()
  for(z in 1:length(types)){
    idx.mks <- which(mks.types == types[z])
    
    MKS <- colnames(onemap_obj$geno)[idx.mks]
    
    diff_geno <- cbind(t(combn(MKS, 2)),mk.type=types[z],NA)
    diff_geno <- data.frame(diff_geno, stringsAsFactors = F)
    colnames(diff_geno) <- c("Marker1", "Marker2", "Marker type", "total")
    
    total <- vector()
    for(i in 1:(dim(diff_geno)[1])){
      total[i] <- sum(sqrt((error_list[[diff_geno[i,1]]] - error_list[[diff_geno[i,2]]])^2) > threshold.probs)
    }
    
    diff_geno[,4] <- total
    bins <- diff_geno[diff_geno$total <= threshold.num,] # How many genotypes are different? Less than the threshold?
    
    if(dim(bins)[1] == 0) {
      warning("No bins were found") 
      return(0)
    } else {
      idx <- match(bins$Marker1,bins$Marker2)
      idx <- idx[-which(is.na(idx))]
      # Some markers can have the differences with some but not with others
      bins[which(bins$Marker1%in%bins$Marker2),1] <- bins$Marker1[idx] 
      
      bins[,1:2] <- lapply(bins[,1:2],as.factor)
      split.bins <- split(bins, bins$Marker1, drop = T)
      bins.all.types <- c(bins.all.types, split.bins)
    }
  }
  return(bins.all.types)
}

##' Not ready yet - use carefully 
##' 
##' @export 
create_data_bins_by_prob <- function(onemap_obj, bins){
  if(!is(bins,"list")) return(onemap_obj)
  
  redun <- unique(as.character(unlist(sapply(bins, "[",2))))
  
  idx <- which(colnames(onemap_obj$geno) %in% redun)
  new.onemap <- onemap_obj
  new.onemap$geno <- onemap_obj$geno[,-idx]
  new.onemap$n.mar <- new.onemap$n.mar - length(idx)
  new.onemap$segr.type <- new.onemap$segr.type[-idx]
  new.onemap$segr.type.num <- new.onemap$segr.type.num[-idx]
  new.onemap$CHROM <- new.onemap$CHROM[-idx]
  new.onemap$POS <- new.onemap$POS[-idx]
  new.onemap$error <- onemap_obj$error[-(idx + rep(c(0:(new.onemap$n.ind-1))*onemap_obj$n.mar, each=length(idx))),]
  
  return(new.onemap)
}
