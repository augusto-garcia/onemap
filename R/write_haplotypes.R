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
globalVariables(c("V1", "V2", "V3", "V4", "P1_H1", "P1_H2", "P2_H1", "P2_H2"))


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
#' group (grp), genotypes probabilities (prob), parent (homologs), and the parents homologs (homologs) 
#' in that homolog and position. 
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
  if(is.null(ind.names)) ind.names <- 1:n.ind
  if(ind[1] == "all"){
    ind <- 1:n.ind
  } 
  
  probs <- lapply(1:length(input.map), function(x) cbind(ind = rep(1:n.ind, each = n.mar[x]),
                                                         grp = group_names[x],
                                                         marker = input.map[[x]]$seq.num,
                                                         pos = c(0,cumsum(kosambi(input.map[[x]]$seq.rf))),
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
    
    # When P1_H1 it means the allele in homolog 1 of the parent 1
    # when H1_P1 is the parent 1 allele in the progeny homolog 1
    if(is(input.map[[1]]$data.name, "outcross")){
      probs <- probs %>%
        mutate(P1_H1 = V1 + V2,
               P1_H2 = V3 + V4,
               P2_H1 = V1 + V3,
               P2_H2 = V2 + V4)  %>%
        select(ind, marker, grp, pos, P1_H1, P1_H2, P2_H1, P2_H2) %>% 
        gather(parents, prob, P1_H1, P1_H2, P2_H1, P2_H2) 
      
      new.col <- t(sapply(strsplit(probs$parents, "_"), "[", 1:2))
      colnames(new.col) <- c("parents", "parents.homologs")
      cross <- "outcross"
    } else {
      probs <- probs %>%
        mutate(H1_P1 = V1 + V2,
               H1_P2 = V3 + V4,
               H2_P1 = V1 + V3,
               H2_P2 = V2 + V4)  %>%
        select(ind, marker, grp, pos, H1_P1, H1_P2, H2_P1, H2_P2) %>% 
        gather(parents, prob, H1_P1, H1_P2, H2_P1, H2_P2) 
      
      new.col <- t(sapply(strsplit(probs$parents, "_"), "[", 1:2))
      colnames(new.col) <- c("progeny.homologs", "parents")
      cross <- "f2"
    }
  } else {
    probs <- do.call(rbind,lapply(probs, function(x) do.call(rbind, x)))
    if(most_likely){
      probs[,5:6] <- t(apply(probs[,5:6], 1, function(x) as.numeric(x == max(x))/sum(x == max(x))))
    }
    
    if (is(input.map[[1]]$data.name, c("backcross")) | is(input.map[[1]]$data.name, c("riself", "risib"))){
      if(is(input.map[[1]]$data.name, c("backcross"))){
        cross <- "backcross"
        probs <- probs %>% 
          mutate(H1_P1 = V1 + V2, # homozygote parent
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
      probs <- probs %>%  select(ind, marker, grp, pos, H1_P1, H1_P2, H2_P1, H2_P2) %>% 
        gather(parents, prob, H1_P1, H1_P2, H2_P1, H2_P2) 
      new.col <- t(sapply(strsplit(probs$parents, "_"), "[", 1:2))
      colnames(new.col) <- c("progeny.homologs", "parents")
    }
  }
  
  probs <- cbind(probs, new.col)
  probs <- cbind(probs[,-5], allele = probs[,5])
  probs <- as.data.frame(probs)
  probs$marker = colnames(input.map[[1]]$data.name$geno)[probs$marker]
  probs$ind <- ind.names[probs$ind]
  
  if(most_likely) flag <- "most.likely" else flag <- "by.probs"
  
  class(probs) <- c("onemap_progeny_haplotypes", cross, "data.frame", flag)
  return(probs)
}

##' Plots progeny haplotypes
##' 
##' Figure is generated with the haplotypes for each selected individual. As a representation, the recombination breakpoints are here considered 
##' to be in the mean point of the distance between two markers.  It is important to highlight that it did not reflects the exact breakpoint position, 
##' specially if the genetic map have low resolution. 
##' 
##' @param x object of class onemap_progeny_haplotypes
##' @param col Color of parentes' homologous.
##' @param position "split" or "stack"; if "split" (default) the alleles' are plotted separately. if "stack" the parents' alleles are plotted together.
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
  
  if(is(x, "outcross")){
    n <- c("H1" = "P1", "H2" = "P2")
    progeny.homologs <- names(n)[match(x$parents, n)]
    probs <- cbind(x, progeny.homologs)
  } else {
    probs <- x
  }
  
  probs <- probs %>% group_by(ind, grp, allele) %>%
    do(rbind(.,.[nrow(.),])) %>%
    do(mutate(.,
              pos2 = c(0,pos[-1]-diff(pos)/2), # Because we don't know exactly where 
              # the recombination occurs, we ilustrate it in the mean point between 
              # markers
              pos = c(pos[-nrow(.)], NA)))
  
  if(is(x, "outcross")){
    p <- ggplot(probs, aes(x = pos, col=allele, alpha = prob)) + ggtitle(main) 
  } else {
    p <- ggplot(probs, aes(x = pos, col=parents, alpha = prob)) + ggtitle(main) 
  }
  
  p <- p +  facet_wrap(~ ind + grp , ncol = ncol) +
    scale_alpha_continuous(range = c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(alpha = "Prob", col = "Allele", x = "position (cM)")
  
  if(is.null(col)) p <- p + scale_color_brewer(palette="Set1")
  else p <- p + scale_color_manual(values = rev(col))
  
  if(position == "stack"){
    p <- p + geom_line(aes(x = pos2, y = progeny.homologs), size = ifelse(show_markers, 4, 5)) + labs(y = "progeny homologs")
    if(show_markers) p <- p + geom_point(aes(y = progeny.homologs), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  if(position == "split"){
    p <- p + geom_line(aes(x = pos2, y = allele), size = ifelse(show_markers, 4, 5))
    if(show_markers) p <- p + geom_point(aes(y = allele), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  return(p)
}

#' Generate graphic with the number of break points for each individual considering the most likely genotypes estimated by the HMM.
#' Genotypes with same probability for two genotypes are removed. By now, only available for outcrossing and f2 intercross. 
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
  # Here we keep the genotype of the marker before it
  doubt <- x[which(x$prob == 0.5),]
  if(dim(doubt)[1] > 0){
    repl <- x[which(x$prob == 0.5)-1,"prob"]
    x[which(x$prob == 0.5),"prob"] <- repl
  }
  
  x <- x[which(x$prob == 1),]
  x <- x[order(x$ind, x$grp, x$prob, x$parents,x$pos),]
  
  x <- x %>% group_by(ind, grp, parents) %>%
    mutate(seq = sequence(rle(as.character(homologs))$length) == 1) %>%
    summarise(counts = sum(seq) -1) %>% ungroup()
  
  class(x) <- c("onemap_progeny_haplotypes_counts", cross, "data.frame")
  return(x)
}


globalVariables(c("counts", "colorRampPalette", "alleles"))


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
                                                  by_homolog = FALSE, 
                                                  n.graphics =NULL, 
                                                  ncol=NULL, ...){
  if(!is(x, "onemap_progeny_haplotypes_counts")) stop("Input need is not of class onemap_progeny_haplotyes_counts")
  
  p <- list()
  n.ind <- length(unique(x$ind))
  nb.cols <- n.ind
  mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
  set.seed(20)
  mycolors <- sample(mycolors)
  
  if(by_homolog){ 
    if(is.null(n.graphics) & is.null(ncol)){
      n.ind <- dim(x)[1]/2
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
      lapply(., function(x) ggplot(x, aes(x=parents, y=counts)) +
               geom_bar(stat="identity", aes(fill=grp)) + theme_minimal() + 
               coord_flip() + 
               scale_fill_manual(values=mycolors) +
               facet_grid(ind ~ ., switch = "y") +
               theme(axis.title.y = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     strip.text.y.left = element_text(angle = 0)) + 
               labs(fill="groups") +
               ylim(0,y_lim_counts)
      )      
  } else {
    x <- x %>% ungroup %>% group_by(ind, grp) %>%
      summarise(counts = sum(counts))
    
    
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
