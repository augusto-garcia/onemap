#' Export OneMap maps to be visualized in VIEWpoly
#' 
#' @param seqs.list a list with `sequence` objects
#' 
#' @return object of class viewmap
#' 
#' @export
export_viewpoly <- function(seqs.list){
  ph.p1 <- ph.p2 <- maps <- list()
  for(i in 1:length(seqs.list)){
    parents <- parents_haplotypes(seqs.list[[i]])
    ph.p1[[i]] <- parents[,c(5,6)]
    ph.p2[[i]] <- parents[,c(7,8)]
    chr <- seqs.list[[i]]$data.name$CHROM[seqs.list[[i]]$seq.num]
    pos <- seqs.list[[i]]$data.name$POS[seqs.list[[i]]$seq.num]
    
    maps[[i]] <- data.frame(mk.names = colnames(seqs.list[[i]]$data.name$geno)[seqs.list[[i]]$seq.num],
                            l.dist = c(0,cumsum(kosambi(seqs.list[[i]]$seq.rf))),
                            g.chr = if(is.null(chr)) rep(NA, length(seqs.list[[i]]$seq.num)) else chr,
                            g.dist = if(is.null(pos)) rep(NA, length(seqs.list[[i]]$seq.num)) else pos,
                            alt = rep(NA, length(seqs.list[[i]]$seq.num)),
                            ref = rep(NA, length(seqs.list[[i]]$seq.num)))  
  }
  
  structure(list(d.p1 = NULL,
                 d.p2 = NULL,
                 ph.p1,
                 ph.p2,
                 maps,
                 software = "onemap"),
            class = "viewmap")
}
