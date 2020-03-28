##' Draw Genotypes
##'
##' Provides a draw of the individuals genotypes.
##'
##'@import ggplot2
##'@import dplyr
##'@import tidyr
##'
##'
##' @param ... Map(s) or list(s) of maps. Object(s) of class sequence.
##' @param ind The individual number to be ploted.
##' @param most.likely logical; if  \code{TRUE}, the most likely genotype is plotted; if FALSE (default) the genotype probability is plotted.
##' @param col Color of parentes' homologous.
##' @param position "split" or "stack"; if "split" (default) the parents' homologous are plotted separately. if "stack" the parents' homologous are plotted together.
##' @param show.markers logical; if  \code{TRUE}, the markers (default) are plotted.
##' @param group.names Names of the groups.
##' @param main An overall title for the plot; default is \code{NULL}.
##'
##' @author Getulio Caixeta Ferreira, \email{getulio.caifer@@gmail.com}
##' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'  data("onemap_example_out")
##'  twopt <- rf_2pts(onemap_example_out)
##'  lg <- group(make_seq(twopt, "all"))
##'  map1 <- map(make_seq(order_seq(input.seq = make_seq(lg, 1), twopt.alg = "rcd"), "force"))
##'  map2 <- map(make_seq(order_seq(input.seq = make_seq(lg, 2), twopt.alg = "rcd"), "force"))
##'  map3 <- map(make_seq(order_seq(input.seq = make_seq(lg, 3), twopt.alg = "rcd"), "force"))
##'  draw_geno(map1, map2, map3, ind = 1)
##'  
##'  data("onemap_example_f2")
##'  twopt <- rf_2pts(onemap_example_f2)
##'  lg <- group(make_seq(twopt, "all"))
##'  map <- list(
##'  map(make_seq(order_seq(input.seq = make_seq(lg, 1),twopt.alg = "rcd"), "force")),
##'  map(make_seq(order_seq(input.seq = make_seq(lg, 2),twopt.alg = "rcd"), "force")),
##'  map(make_seq(order_seq(input.seq = make_seq(lg, 3),twopt.alg = "rcd"), "force"))
##'  )
##'  draw_geno(map, 1:2)
##'
##' }
##'@export

draw_geno <- function(..., 
                      ind = 1, 
                      most.likely = FALSE, 
                      col = NULL, 
                      position = "stack", 
                      show.markers = TRUE, 
                      group.names = NULL, 
                      main = "Genotypes"){
  #input map
  input <- list(...)
  if(length(input) == 0) stop("argument '...' missing, with no default")
  # Accept list of sequences or list of list of sequences
  if(is(input[[1]], "sequence")) input.map <- input else input.map <- unlist(input, recursive = FALSE)
  if(!all(sapply(input.map, function(x) is(x, "sequence")))) stop(paste("Input objects must be of 'sequence' class. \n"))
  if(is.null(group.names)) group.names <- paste("Group",seq(input.map), sep = " - ")
  n.mar <- sapply(input.map, function(x) length(x$seq.num))
  n.ind <- sapply(input.map, function(x) ncol(x$probs))/n.mar
  ind.select <- ind
  phase <- list('1' = c(1,2,3,4),
                '2' = c(2,1,4,3),
                '3' = c(3,4,1,2),
                "4" = c(4,3,2,1))
  
  seq.phase <- lapply(input.map, function(x) c(1,x$seq.phases))
  
  probs <- lapply(1:length(input.map), function(x) cbind(ind = rep(1:n.ind[x], each = n.mar[x]),
                                                         grp = group.names[x],
                                                         marker = input.map[[x]]$seq.num,
                                                         pos = c(0,cumsum(kosambi(input.map[[x]]$seq.rf))),
                                                         as.data.frame(t(input.map[[x]]$probs))))
  probs <- lapply(probs, function(x) split.data.frame(x, x$ind)[ind])
  
  for(g in seq_along(input.map)){
    for(i in seq_along(ind)){
      for(m in seq(n.mar[g])){
        probs[[g]][[i]][m:n.mar[g],5:8] <- probs[[g]][[i]][m:n.mar[g],phase[[seq.phase[[g]][m]]]+4]
      }
    }
  }
  probs <- do.call(rbind,lapply(probs, function(x) do.call(rbind, x)))
  if(most.likely){
    probs[,5:8] <- t(apply(probs[,5:8], 1, function(x) as.numeric(x == max(x))/sum(x == max(x))))
  }
  probs <- probs %>%
    group_by(ind, grp) %>%
    do(rbind(.,.[nrow(.),])) %>% 
    do(mutate(.,
              pos2 = c(0,pos[-1]-diff(pos)/2),
              pos = c(pos[-nrow(.)], NA))) %>% 
    mutate(P1_1 = V1 + V2,
           P1_2 = V3 + V4,
           P2_1 = V1 + V3,
           P2_2 = V2 + V4) %>% 
    select(ind, grp, pos, pos2, P1_1, P1_2, P2_1, P2_2) %>% 
    gather(allele, prob, P1_1, P1_2, P2_1, P2_2) %>% 
    mutate(homolog = factor(if_else(allele %in% c("P1_1","P1_2"), "P1", "P2"), levels = c("P2","P1")),
           allele = factor(allele, levels = c("P2_2", "P2_1", "P1_2", "P1_1")))
  
  p <- ggplot(probs, aes(x = pos, col = allele, alpha = prob)) + ggtitle(main) +
    facet_grid(grp ~ ind,switch = "y") +
    scale_alpha_continuous(range = c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE))
  
  if(is.null(col)) p <- p + scale_color_brewer(palette="Set1")
  else p <- p + scale_color_manual(values = rev(col))
  
  if(position == "stack"){
    p <- p + geom_line(aes(x = pos2, y = homolog), size = ifelse(show.markers, 4, 5))
    if(show.markers) p <- p + geom_point(aes(y = homolog), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  if(position == "split"){
    p <- p + geom_line(aes(x = pos2, y = allele), size = ifelse(show.markers, 4, 5))
    if(show.markers) p <- p + geom_point(aes(y = allele), size = 5, stroke = 2, na.rm = T, shape = "|")
  } 
  return(p)
}

#if(phase = "CC"){}

#if(phase = "CR"){
#  1 = 2
#  2 = 1
#  3 = 4
#  4 = 3
#}
#if(phase = "RC"){
#  1 = 3
#  2 = 4
#  3 = 1
#  4 = 2
#}
#if(phase = "RR"){
#  1 = 4
#  2 = 3
#  3 = 2
#  4 = 1
#}


