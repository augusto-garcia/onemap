#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: write_haplotypes.R                                            ##
## Contains: parents_haplotypes, progeny_haplotypes,                   ##
## plot.onemap_progeny_haplotypes,                                     ##
## plot.onemap_progeny_haplotypes_counts                               ##
##                                                                     ##
## Written by  Getulio Caixeta Ferreira and Cristiane Taniguti         ##
##                                                                     ##
## First version: 2020/05/26                                           ##
## Last update: 2020/05/26                                             ##
## License: GNU General Public License version 3 or later              ##
##                                                                     ##
#######################################################################

globalVariables(c("grp", "for.split", ".", "pos", "prob", "pos2", "homologs"))
globalVariables(c("V1", "V2", "V3", "V4", "H1_P1", "H1_P2", "H2_P1", "H2_P2"))


#' Generates data.frame with parents estimated haplotypes 
#'
#' @param ... objects of class sequence
#' @param group_names vector of characters defining the group names
#'
#' @return data.frame with group ID (group), marker number (mk.number) 
#' and names (mk.names), position in centimorgan (dist) and parents haplotypes 
#' (P1_1, P1_2, P2_1, P2_2)
#' 
#' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
#' @export
parents_haplotypes <- function(..., group_names=NULL){
  input <- list(...)
  if(length(input) == 0) stop("argument '...' missing, with no default")
  # Accept list of sequences or list of list of sequences
  if(is(input[[1]], "sequence")) input.map <- input else input.map <- unlist(input, recursive = FALSE)
  if(!all(sapply(input.map, function(x) is(x, "sequence")))) stop(paste("Input objects must be of 'sequence' class. \n"))
  if(is.null(group_names)) group_names <- paste("Group",seq(input.map), sep = " - ")
  
  if(all(sapply(input, function(x) is(x, "sequence")))){
    n <- length(sapply(input, function(x) is(x, "sequence")))
  } else n <- 1
  
  input_temp <- input
  out_dat <- data.frame()
  for(z in 1:n){
    if(all(sapply(input_temp, function(x) is(x, "sequence")))) input <- input_temp[[z]]
    marnames <- colnames(input$data.name$geno)[input$seq.num]
    if(length(input$seq.rf) == 1 && input$seq.rf == -1) {
      # no information available for the order
      cat("\nParameters not estimated.\n\n")
    }
    else {
      # convert numerical linkage phases to strings
      link.phases <- matrix(NA,length(input$seq.num),2)
      link.phases[1,] <- rep(1,2)
      for (i in 1:length(input$seq.phases)) {
        switch(EXPR=input$seq.phases[i],
               link.phases[i+1,] <- link.phases[i,]*c(1,1),
               link.phases[i+1,] <- link.phases[i,]*c(1,-1),
               link.phases[i+1,] <- link.phases[i,]*c(-1,1),
               link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
        )
      }
      
      ## display results
      marnumbers <- input$seq.num
      distances <- c(0,cumsum(get(get(".map.fun", envir=.onemapEnv))(input$seq.rf)))
      ## whith diplotypes for class 'outcross'
      if(is(input$data.name, c("outcross", "f2"))){
        ## create diplotypes from segregation types and linkage phases
        link.phases <- apply(link.phases,1,function(x) paste(as.character(x),collapse="."))
        parents <- matrix("",length(input$seq.num),4)
        for (i in 1:length(input$seq.num))
          parents[i,] <- return_geno(input$data.name$segr.type[input$seq.num[i]],link.phases[i])
        out_dat_temp <- data.frame(group= group_names[z], mk.number = marnumbers, mk.names = marnames, dist = as.numeric(distances), 
                                   P1_1 = parents[,1],
                                   P1_2 = parents[,2],
                                   P2_1 = parents[,3],
                                   P2_2 = parents[,4])
        out_dat <- rbind(out_dat, out_dat_temp)
      }
      ## whithout diplotypes for other classes
      else if(is(input$data.name, c("backcross", "riself", "risib"))){
        cat("There is only a possible phase for this cross type\n")
      }
      else warning("invalid cross type")
    }
  }
  return(out_dat)
}


#' Generate data.frame with genotypes estimated by HMM and its probabilities
#'
#' @param ... Map(s) or list(s) of maps. Object(s) of class sequence.
#' @param ind vector with individual index to be evaluated or "all" to include all individuals
#' @param most_likely logical; if  \code{TRUE}, the most likely genotype receive 1 and all the rest 0. 
#' If there are more than one most likely both receive 0.5.
#' if FALSE (default) the genotype probability is plotted.
#' @param group_names Names of the groups.
#' 
#' @return a data.frame information: individual (ind) and group (grp) ID, position in centimorgan (pos), 
#' position adjusted for build graphic (pos2), f1 homologs (f1), progeny homologs (prog), 
#' parents homologs (parents) IDs and genotypes probabilities (prob)
#' 
#' @import dplyr
#' @import tidyr
#' 
#' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
#' @export
progeny_haplotypes <- function(...,
                               ind = 1, 
                               group_names=NULL,
                               most_likely = FALSE){
  #input map
  input <- list(...)
  if(length(input) == 0) stop("argument '...' missing, with no default")
  # Accept list of sequences or list of list of sequences
  if(is(input[[1]], "sequence")) input.map <- input else input.map <- unlist(input, recursive = FALSE)
  if(!all(sapply(input.map, function(x) is(x, "sequence")))) stop(paste("Input objects must be of 'sequence' class. \n"))
  if(is.null(group_names)) group_names <- paste("Group",seq(input.map), sep = " - ")
  n.mar <- sapply(input.map, function(x) length(x$seq.num))
  n.ind <- sapply(input.map, function(x) ncol(x$probs))/n.mar
  ind.names <- lapply(input.map, function(x) rownames(x$data.name$geno))
  ind.names <- unique(unlist(ind.names)) 
  if(length(unique(n.ind)) != 1) stop("At least one of the sequences have different number of individuals in dataset.")
  n.ind <- unique(n.ind)
  if(ind[1] == "all"){
    ind <- 1:n.ind
  } 
  
  probs <- lapply(1:length(input.map), function(x) cbind(ind = rep(1:n.ind, each = n.mar[x]),
                                                         grp = group_names[x],
                                                         marker = input.map[[x]]$seq.num,
                                                         pos = c(0,cumsum(get(get(".map.fun", envir=.onemapEnv))(input.map[[x]]$seq.rf))),
                                                         as.data.frame(t(input.map[[x]]$probs))))
  probs <- lapply(probs, function(x) split.data.frame(x, x$ind)[ind])
  
  if(is(input.map[[1]]$data.name, "outcross") | is(input.map[[1]]$data.name, "f2")){
    phase <- list('1' = c(1,2,3,4),
                  '2' = c(2,1,4,3),
                  '3' = c(3,4,1,2),
                  "4" = c(4,3,2,1))
    
    seq.phase <- lapply(input.map, function(x) c(1,x$seq.phases))
    
    # Adjusting phases
    for(g in seq_along(input.map)){
      for(i in seq_along(ind)){
        for(m in seq(n.mar[g])){
          probs[[g]][[i]][m:n.mar[g],5:8] <- probs[[g]][[i]][m:n.mar[g],phase[[seq.phase[[g]][m]]]+4]
        }
      }
    }
    
    # If there are two most probably genotypes both will receive 0.5
    probs <- do.call(rbind,lapply(probs, function(x) do.call(rbind, x)))
    
    if(most_likely){
      probs[,5:8] <- t(apply(probs[,5:8], 1, function(x) as.numeric(x == max(x))/sum(x == max(x))))
    }
    
    probs <- probs %>%
      mutate(H1_P1 = V1 + V2,
             H1_P2 = V3 + V4,
             H2_P1 = V1 + V3,
             H2_P2 = V2 + V4) 
    
    if(is(input.map[[1]]$data.name, "outcross")){
      cross <- "outcross"
    } else {
      cross <- "f2"
    }
    
  } else {
    #probs <- probs1
    probs <- do.call(rbind,lapply(probs, function(x) do.call(rbind, x)))
    if(most_likely){
      probs[,5:6] <- t(apply(probs[,5:6], 1, function(x) as.numeric(x == max(x))/sum(x == max(x))))
    }
    
    if (is(input.map[[1]]$data.name, c("backcross"))){
      cross <- "backcross"
      probs <- probs %>% 
        mutate(H1_P1 = V1 + V2, # homozigote parent
               H1_P2 = 0,
               H2_P1 = V2, 
               H2_P2 = V1) 
      
    } else if (is(input.map[[1]]$data.name, c("riself", "risib"))){
      cross <- "rils"
      probs <- probs %>%
        mutate(H1_P1 = V1,
               H1_P2 = V2,
               H2_P1 = V1,
               H2_P2 = V2)
    }
  }
  
  probs <- probs %>% 
    select(ind, grp, pos, H1_P1, H1_P2, H2_P1, H2_P2) %>% 
    gather(homologs, prob, H1_P1, H1_P2, H2_P1, H2_P2) 
  
  new.col <- t(sapply(strsplit(probs$homologs, "_"), "[", 1:2))
  colnames(new.col) <- c("homologs", "parents")
  
  probs <- cbind(probs, new.col)
  probs <- probs[,-4]
  probs <- as.data.frame(probs)
  
  probs$ind <- ind.names[probs$ind]
  
  if(most_likely) flag <- "most.likely" else flag <- "by.probs"
  
  class(probs) <- c("onemap_progeny_haplotypes", cross, "data.frame", flag)
  return(probs)
}

##' Plots progeny haplotypes
##' 
##' @param x object of class onemap_progeny_haplotypes
##' @param col Color of parentes' homologous.
##' @param position "split" or "stack"; if "split" (default) the parents' homologous are plotted separately. if "stack" the parents' homologous are plotted together.
##' @param show_markers logical; if  \code{TRUE}, the markers (default) are plotted.
##' @param main An overall title for the plot; default is \code{NULL}.
##' @param ncol number of columns of the facet_wrap
##' @param ... currently ignored
##' 
##' @method plot onemap_progeny_haplotypes
##' @import ggplot2
#' @import dplyr
#' @import tidyr
#' 
##' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' 
##' @export
plot.onemap_progeny_haplotypes <- function(x,
                                           col = NULL, 
                                           position = "stack",
                                           show_markers = TRUE, 
                                           main = "Genotypes", ncol=4, ...){
  
  colors <- ifelse(is(x,"outcross"), "for.split", "parents")  
  
  probs <- cbind(x, for.split= paste0(x$homologs, "_", x$parents))
  
  probs <- probs %>% group_by(ind, grp, for.split) %>%
    do(rbind(.,.[nrow(.),])) %>%
    do(mutate(.,
              pos2 = c(0,pos[-1]-diff(pos)/2), # Because we don't know exactly where 
              # the recombination occurs, we ilustrate it in the mean point between 
              # markers
              pos = c(pos[-nrow(.)], NA)))
  
  p <- ggplot(probs, aes(x = pos, col=get(colors), alpha = prob)) + ggtitle(main) +
    facet_wrap(~ ind + grp , ncol = ncol) +
    scale_alpha_continuous(range = c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(alpha = "Prob", col = "Allele", x = "position (cM)")
  
  if(is.null(col)) p <- p + scale_color_brewer(palette="Set1")
  else p <- p + scale_color_manual(values = rev(col))
  
  if(position == "stack"){
    p <- p + geom_line(aes(x = pos2, y = homologs), size = ifelse(show_markers, 4, 5)) + labs(y = "homologs")
    if(show_markers) p <- p + geom_point(aes(y = homologs), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  if(position == "split"){
    p <- p + geom_line(aes(x = pos2, y = for.split), size = ifelse(show_markers, 4, 5))
    if(show_markers) p <- p + geom_point(aes(y = for.split), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  return(p)
}


#' Convert phased vcf to onemap_progeny_haplotypes object. Turns possible to use
#' the plot.onemap_progeny_haplotypes (By now only for outcrossing population).
#' 
#' @param vcfR.object vcfR object
#' @param ind.id vector of characters with progeny individuals to be evaluated
#' @param group_names string with chromosomes or group names to be evaluated
#' @param parent1 character with parent 1 ID
#' @param parent2 character with parent 2 ID
#' @param crosstype character defining crosstype (outcross, f2 intercross, f2 backcross, ril sib, ril self)
#' 
#' @export
vcf2progeny_haplotypes <- function(vcfR.object, 
                                   ind.id=NULL, 
                                   group_names = NULL, 
                                   parent1, 
                                   parent2,
                                   crosstype){
  if(is.null(ind.id)){
    stop("You should define one individual.")
  } 
  n.mk <- dim(vcfR.object@gt)[1]
  n.ind <- dim(vcfR.object@gt)[2]-1
  INDS <- dimnames(vcfR.object@gt)[[2]][-1]
  
  if(length(which(INDS %in% ind.id)) != length(ind.id)){
    stop("At least one of the individuals in ind.id was not found in vcfR object.")
  }
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcfR.object@fix[,1],"_", vcfR.object@fix[,2])
  
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcfR.object@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcfR.object@gt[,i], split=":"), "[[", GT))
  
  CHROM <- vcfR.object@fix[,1]
  
  if(is.null(group_names)) group_names <- CHROM[1]
  
  if(length(which(unique(CHROM) %in% group_names)) != length(group_names)){
    stop("At least one of the groups in group_names was not found in vcfR object.")
  }
  
  progeny_haplotypes_obj_chr <- data.frame()
  for(chr in 1:length(group_names)){ ### Need optimization
    CHROM.now <- which(CHROM %in% group_names[chr])
    POS <- as.numeric(vcfR.object@fix[,2])[CHROM.now]
    
    colnames(GT_matrix) <- INDS
    rownames(GT_matrix) <- MKS
    
    P1.idx <- grep(parent1, INDS)
    P2.idx <- grep(parent2, INDS)
    
    P1_1 <- sapply(strsplit(GT_matrix[CHROM.now,P1.idx], "[|]"), "[",1)
    P1_2 <- sapply(strsplit(GT_matrix[CHROM.now,P1.idx], "[|]"), "[",2)
    P2_1 <- sapply(strsplit(GT_matrix[CHROM.now,P2.idx], "[|]"), "[",1)
    P2_2 <- sapply(strsplit(GT_matrix[CHROM.now,P2.idx], "[|]"), "[",2)
    
    progeny_haplotypes_obj_ind <- data.frame()
    for(ind in 1:length(ind.id)){
      ind.idx <- grep(ind.id[ind], INDS)
      ind.number <- grep(ind.id[ind], INDS[-c(P1.idx, P2.idx)])
      ind_1 <- sapply(strsplit(GT_matrix[CHROM.now,ind.idx], "[|]"), "[",1)
      ind_2 <- sapply(strsplit(GT_matrix[CHROM.now,ind.idx], "[|]"), "[",2)
      
      p.names <- c("P1", "P2", "P1", "P2")
      Hs <-  rep(list(rep(NA, length(ind_1))),2)
      progeny_haplotypes_obj <- data.frame()
      for(w in 1:2){
        ref.frags<-rep(1, length(ind_1))
        comp <- list(cbind(P1_1, P1_2, P2_1,P2_2,H1=ind_1, H2=ind_2))
        idx.cum <- 1
        while(any(is.na(Hs[[w]]))){
          # count how many equal characters are consecutive
          frags <- rep(list(list(0,0,0,0)),length(comp))
          for(z in 1:length(comp)){
            for(j in 1:4){
              if(is.vector(comp[[z]]))
                comp[[z]] <- t(as.matrix(comp[[z]]))
              idx.comp <- comp[[z]][,j] == comp[[z]][,4+w]
              frags[[z]][[j]] <- sequence(rle(as.character(idx.comp))$length)
            }
          }
          # Find the higher fragment in ind1
          new.zs <- list()
          inter.cum <-max(ref.frags)
          for(z in 1:length(frags)){
            max.ind1 <- unlist(lapply(frags[[z]], max))
            # I could't adapt to inbred because it found two possible match
            idx.ind1 <- which.max(max.ind1)
            which.max.ind1 <- which.max(frags[[z]][[idx.ind1]])
            frag <- (which.max.ind1 - max.ind1[idx.ind1]+1):which.max.ind1 
            Hs[[w]][ref.frags==idx.cum][frag] <- p.names[idx.ind1]
            # The fragment should be removed and the split the remaining in two
            ref.frags.new <- ref.frags
            if(frag[1] == 1)
              frag1 <- NULL else  frag1 <- comp[[z]][1:(frag[1]-1),]
            if(frag[length(frag)] == dim(comp[[z]])[1])
              frag2 <- NULL else  frag2 <- comp[[z]][(frag[length(frag)]+1):dim(comp[[z]])[1],]
            if(!is.null(frag1)){
              new.zs <- c(new.zs,list(frag1))
              inter.cum <- inter.cum + 1
              ref.frags.new[ref.frags==idx.cum][1:(frag[1]-1)] <- inter.cum
            }
            
            if(!is.null(frag2)){
              new.zs <- c(new.zs,list(frag2))
              inter.cum <- inter.cum + 1
              ref.frags.new[ref.frags==idx.cum][(frag[length(frag)]+1):length(ref.frags.new[ref.frags==idx.cum])] <- inter.cum
            }
            ref.frags <- ref.frags.new
            idx.cum <- idx.cum + 1
          }
          comp <- new.zs
        }
        
        num.mk <- length(CHROM.now)
        df.H <- data.frame(ind = rep(ind.id[ind], 2*num.mk),
                           grp = rep(CHROM[CHROM.now], 2),
                           pos = rep(POS, 2),
                           prob = rep(0, 2*num.mk),
                           homologs = rep(c("H1", "H2")[w], num.mk),
                           parents = rep(rep(c("P1","P2"),each =num.mk)))
        
        for(i in 1:length(Hs[[w]])){
          df.H$prob[which(df.H$pos == POS[i] & df.H$parents == Hs[[w]][i])] <- 1
        }
        # bind homologs
        progeny_haplotypes_obj <- rbind(progeny_haplotypes_obj, df.H)
      }
      # bind individuals
      progeny_haplotypes_obj_ind <- rbind(progeny_haplotypes_obj_ind, progeny_haplotypes_obj)
    }
    # bind chromosomes
    progeny_haplotypes_obj_chr <- rbind(progeny_haplotypes_obj_chr,progeny_haplotypes_obj_ind)
  }
  
  crosstype <- switch(crosstype, "outcross" = "outcross", "f2 intercross"="f2",
                      "f2 backcross"="backcross", "ril sib"="rils", "ril self"="rils")
  
  flag <- "most.likely"
  class(progeny_haplotypes_obj_chr) <- c("onemap_progeny_haplotypes", crosstype, "data.frame", flag)
  return(progeny_haplotypes_obj_chr)
}


#' By now only for outcrossing
#' 
#' Genotypes with same probability for two genotypes will be removed
#' 
#' @param x object of class onemap_progeny_haplotypes
#' 
#' @import dplyr
#' @import tidyr
#'@export
progeny_haplotypes_counts <- function(x){
  if(!is(x, "onemap_progeny_haplotypes")) stop("Input need is not of class onemap_progeny_haplotyes")
  if(!is(x, "most.likely")) stop("The most likely genotypes must receive maximum probability (1)")
  cross <- class(x)[2]
  
  # Some genotypes receveis prob of 0.5, here we need to make a decision about them
  doubt <- x[which(x$prob == 0.5),]
  if(dim(doubt)[1] > 0){
    nondupli <- which(!duplicated(paste0(doubt$ind, doubt$grp, doubt$pos)))
    x[which(x$prob == 0.5),][nondupli,]$prob <- 1
  }
  
  x <- x[which(x$prob == 1),]
  x <- x[order(x$ind, x$grp, x$prob, x$homologs,x$pos),]
  
  x <- x %>% group_by(ind, grp, homologs) %>%
    mutate(seq = sequence(rle(as.character(parents))$length) == 1) %>%
    summarise(counts = sum(seq) -1) %>% ungroup()
  
  class(x) <- c("onemap_progeny_haplotypes_counts", cross, "data.frame")
  return(x)
}


globalVariables(c("counts", "colorRampPalette", "parents"))


##' Plot recombination breakpoints counts for each individual
##'
##' @param x object of class onemap_progeny_haplotypes_counts
##' @param by_homolog logical, if TRUE plots counts by homolog (two for each individuals), if FALSE plots total counts by individual
##' @param n.graphics integer defining the number of graphics to be plotted, they separate the individuals in different plots 
##' @param ncol integer defining the number of columns in plot
##' @param ... currently ignored
##' 
##' @method plot onemap_progeny_haplotypes_counts
##' @import ggplot2
##' @importFrom ggpubr ggarrange
##' @import dplyr
##' @import tidyr
##' @importFrom RColorBrewer brewer.pal
##' @importFrom grDevices colorRamp colorRampPalette
##' 
##' @export
plot.onemap_progeny_haplotypes_counts <- function(x, 
                                                  by_homolog = FALSE, # Do not use TRUE yet
                                                  n.graphics =NULL, 
                                                  ncol=NULL, ...){
  if(!is(x, "onemap_progeny_haplotypes_counts")) stop("Input need is not of class onemap_progeny_haplotyes_counts")
  p <- list()
  if(by_homolog){ ## Bugfix! 
    if(is.null(n.graphics) & is.null(ncol)){
      n.ind <- dim(x)[1]
      if(n.ind/25 <= 1) {
        n.graphics = 1
        ncol=1 
      }else { n.graphics = round(n.ind/25,0)
      ncol=round(n.ind/25,0)
      }
    }
    size <- dim(x)[1]
    if(size%%n.graphics == 0){
      div.n.graphics <- rep(1:n.graphics, each= size/n.graphics) 
    } else {           
      div.n.graphics <- c(rep(1:n.graphics, each = round(size/n.graphics,0)), rep(n.graphics, size%%n.graphics))
    }
    
    y_lim_counts <- max(x$counts)
    div.n.graphics <- div.n.graphics[1:size]
    p <- x %>% mutate(div.n.graphics = div.n.graphics) %>%
      split(., .$div.n.graphics) %>%
      lapply(., function(x) ggplot(x, aes(x=homologs, y=counts)) +
               geom_bar(stat="identity", aes(fill=grp)) + theme_minimal() + 
               coord_flip() + 
               scale_fill_brewer(palette="Set1") +
               facet_grid(ind~., switch = "y") +
               theme(axis.title.y = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     strip.text.y.left = element_text(angle = 0)) + 
               labs(fill="groups") +
               ylim(0,y_lim_counts)
      )      
  } else {
    x <- x %>% group_by(ind, grp) %>%
      summarise(counts = sum(counts))
    
    n.ind <- length(unique(x$ind))
    if(is.null(n.graphics) & is.null(ncol)){
      if(n.ind/25 <= 1) {
        n.graphics = 1
        ncol=1 
      }else { 
        n.graphics = round(n.ind/25,0)
        ncol=round(n.ind/25,0)
      }
    }
    
    size <-n.ind
    if(size%%n.graphics == 0){
      div.n.graphics <- rep(1:n.graphics, each= size/n.graphics) 
    } else {           
      div.n.graphics <-   c(rep(1:n.graphics, each = round(size/n.graphics,0)),rep(n.graphics, size%%n.graphics))
    }
    div.n.graphics <- div.n.graphics[1:n.ind]
    div.n.graphics <- rep(div.n.graphics, each = length(unique(x$grp)))
    
    x$ind <- factor(as.character(x$ind), levels = sort(as.character(unique(x$ind))))
    
    temp <- x %>% ungroup() %>% group_by(ind) %>%
      summarise(total = sum(counts))
      
    y_lim_counts <- max(temp$total)
    nb.cols <- n.ind
    mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
    set.seed(20)
    mycolors <- sample(mycolors)
    p <- x %>% ungroup() %>%  mutate(div.n.graphics = div.n.graphics) %>%
      split(., .$div.n.graphics) %>%
      lapply(., function(x) ggplot(x, aes(x=ind, y=counts, fill=grp)) +
               geom_bar(stat="identity") + coord_flip() + 
               scale_fill_manual(values=mycolors) +
               theme(axis.title.y = element_blank(), 
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
               labs(fill="groups") +
               ylim(0,y_lim_counts)
      ) 
  }
  p <- ggarrange(plotlist = p, common.legend = T, label.x = 1, ncol = ncol, nrow = round(n.graphics/ncol,0))
  return(p)
}
