#' Export OneMap maps to be visualized in VIEWpoly
#' 
#' @param seqs.list a list with `sequence` objects
#' 
#' @return object of class viewmap
#' 
#' @export
export_viewpoly <- function(seqs.list){
  if(!(inherits(seqs.list,c("list", "sequence")))) stop(deparse(substitute(seqs.list))," is not an object of class 'list' or 'sequence'")
  
  ph.p1 <- ph.p2 <- maps <- d.p1 <- d.p2 <- list()
  for(i in 1:length(seqs.list)){
    
    # only for biallelic markers
    types <- seqs.list[[i]]$data.name$segr.type[seqs.list[[i]]$seq.num]
    if(inherits(seqs.list[[i]]$data.name, "outcross")){
      labs.p1 <- c("B3.7"=1, "D1.10" = 1, "D2.15" = 0)
      labs.p2 <- c("B3.7"=1, "D1.10" = 0, "D2.15" = 1)
      d.p1[[i]] <- labs.p1[match(types, names(labs.p1))]
      d.p2[[i]] <- labs.p2[match(types, names(labs.p2))]
    } else if(inherits(seqs.list[[i]]$data.name, "intercross")){
      labs <- c("A.H.B" = 1)
      d.p1[[i]] <- labs[match(types, names(labs))]
      d.p2[[i]] <- labs[match(types, names(labs))]
    } else if(inherits(seqs.list[[i]]$data.name, "backcross")){
      labs <- c("A.H" = 1)
      d.p1[[i]] <- labs[match(types, names(labs))]
      labs <- c("A.H" = 0)
      d.p2[[i]] <- labs[match(types, names(labs))]
    } else if(inherits(seqs.list[[i]]$data.name, "ri")){
      labs <- c("A.B" = 1)
      d.p1[[i]] <- labs[match(types, names(labs))]
      labs <- c("A.B" = 2)
      d.p2[[i]] <- labs[match(types, names(labs))]
    }
    
    parents <- parents_haplotypes(seqs.list[[i]])
    ph.p1[[i]] <- parents[,c(5,6)]
    ph.p2[[i]] <- parents[,c(7,8)]
    rownames(ph.p1[[i]]) <- rownames(ph.p2[[i]]) <- colnames(seqs.list[[i]]$data.name$geno)[seqs.list[[i]]$seq.num]
    names(ph.p1[[i]]) <- c("a", "b")
    names(ph.p2[[i]]) <- c("c", "d")
    
    chr <- seqs.list[[i]]$data.name$CHROM[seqs.list[[i]]$seq.num]
    pos <- seqs.list[[i]]$data.name$POS[seqs.list[[i]]$seq.num]
    
    maps[[i]] <- data.frame(mk.names = colnames(seqs.list[[i]]$data.name$geno)[seqs.list[[i]]$seq.num],
                            l.dist = c(0,cumsum(kosambi(seqs.list[[i]]$seq.rf))),
                            g.chr = if(is.null(chr)) rep(NA, length(seqs.list[[i]]$seq.num)) else chr,
                            g.dist = if(is.null(pos)) rep(NA, length(seqs.list[[i]]$seq.num)) else pos,
                            alt = rep(NA, length(seqs.list[[i]]$seq.num)),
                            ref = rep(NA, length(seqs.list[[i]]$seq.num)))  
  }
  
  structure(list(d.p1 = d.p1,
                 d.p2 = d.p2,
                 ph.p1 = ph.p1,
                 ph.p2= ph.p2,
                 maps = maps,
                 software = "onemap"),
            class = "viewmap")
}

#' Export genotype probabilities in MAPpoly format (input for QTLpoly)
#' 
#' @param input.map object of class `sequence`
#' 
#' @return object of class `mappoly.genoprob`
#' 
#' @export
export_mappoly_genoprob <- function(input.map){
  if(!(inherits(input.map,c("sequence")))) stop(deparse(substitute(seqs.list))," is not an object of class 'sequence'")
  
  probs <- cbind(ind = rep(1:input.map$data.name$n.ind, each = length(input.map$seq.num)),
                 marker = rep(colnames(input.map$data.name$geno)[input.map$seq.num], input.map$data.name$n.ind),
                 pos = c(0,cumsum(kosambi(input.map$seq.rf))),
                 as.data.frame(t(input.map$probs)))
  
  if(inherits(input.map$data.name, "outcross") | inherits(input.map$data.name, "f2")){
    phase <- list('1' = c(1,2,3,4),
                  '2' = c(2,1,4,3),
                  '3' = c(3,4,1,2),
                  "4" = c(4,3,2,1))
    
    seq.phase <- rep(c(1,input.map$seq.phases), input.map$data.name$n.ind)
    
    # Adjusting phases
    for(i in 1:length(seq.phase))
      probs[i,4:7] <- probs[i,phase[[seq.phase[i]]]+3]
  }
  
  colnames(probs)[4:7] <- c("a:c", "a:d", "b:c", "b:d")
  
  genoprob <- array(unlist(t(probs[,4:7])), 
                    dim = c(4, length(input.map$seq.num), input.map$data.name$n.ind),
                    dimnames = list(c("a:c", "a:d", "b:c", "b:d"),
                                    colnames(input.map$data.name$geno)[input.map$seq.num],
                                    rownames(input.map$data.name$geno)))
  
  map <- cumsum(c(0,kosambi(input.map$seq.rf)))
  names(map) <- colnames(input.map$data.name$geno)[input.map$seq.num]
  structure(list(probs = genoprob,
                 map = map), 
            class = "mappoly.genoprob")
}
