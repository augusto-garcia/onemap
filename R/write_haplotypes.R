#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: write_haplotypes.R                                            ##
## Contains: parents_haplotypes, progeny_haplotypes,                   ##
## plot.onemap_progeny_haplotypes                                      ##
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
#' @param ind vector with individual index to be evaluated
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
  ind.select <- ind
  
  probs <- lapply(1:length(input.map), function(x) cbind(ind = rep(1:n.ind[x], each = n.mar[x]),
                                                         grp = group_names[x],
                                                         marker = input.map[[x]]$seq.num,
                                                         #pos = c(0,cumsum(get(get(".map.fun", envir=.onemapEnv))(input.map[[x]]$seq.rf))),
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
  class(probs) <- c("onemap_progeny_haplotypes", cross, "data.frame")
  return(probs)
}

##' Plots progeny haplotypes
##' 
##' @param x object of class onemap_progeny_haplotypes
##' @param col Color of parentes' homologous.
##' @param position "split" or "stack"; if "split" (default) the parents' homologous are plotted separately. if "stack" the parents' homologous are plotted together.
##' @param show_markers logical; if  \code{TRUE}, the markers (default) are plotted.
##' @param main An overall title for the plot; default is \code{NULL}.
##' @param ... currently ignored
##' 
##' @method plot onemap_progeny_haplotypes
##' @import ggplot2
##' 
##' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' 
##' @export
plot.onemap_progeny_haplotypes <- function(x,
                                           col = NULL, 
                                           position = "stack",
                                           show_markers = TRUE, 
                                           main = "Genotypes", ...){
  
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
    facet_grid(grp ~ ind,switch = "y") +
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
