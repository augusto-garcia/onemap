##' Draw Genotypes
##'
##' Provides a draw of the individuals genotypes.
##'
##'@import ggplot2
##'@import dplyr
##'
##' @param input.map map(s). Object(s) of class \code{sequence}.
##' @param ind Individual number
##' @param most.likely TRUE or FALSE
##' @param position "split" or "stack"
##' @param main an overall title for the plot. Default is \code{NULL}.
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

draw_geno <- function(input.map, ind, most.likely = F, position = "stack", main = "Genotypes"){
  probs <- as.data.frame(t(input.map$probs))
  n.mar <- length(input.map$seq.num)
  n.ind <- nrow(probs)/n.mar
  ind.select <- ind
  if(most.likely){
    probs <- as.data.frame(t(apply(probs, 1, function(x) as.numeric(x == max(x))))  )
  }
  p <- cbind(expand.grid(mar = c(0,cumsum(kosambi(input.map$seq.rf))),ind = 1:n.ind), probs) %>% 
    filter(ind %in% ind.select) %>% 
    mutate(a = V1 + V2,
           b = V3 + V4,
           c = V1 + V3,
           d = V2 + V4) %>% 
    select(ind, mar, a, b, c, d) %>% 
    gather(allele, prob, a, b, c, d) %>% 
    mutate(homolog = factor(if_else(allele %in% c("a","b"), 1, 2), levels = c(2,1)),
           allele = factor(allele, levels = c("d","c","b","a")),) %>% 
    ggplot() +
    #scale_color_manual(values = c("red","darkgreen","darkblue","orange"))+
    facet_wrap(~ind, nrow = ceiling(length(ind.select)^.5)) + ggtitle(main)
  if(position == "stack") p <- p + geom_line(aes(x = mar, y = homolog, col = allele, alpha = prob),size = 5, show.legend = F)
  if(position == "split") p <- p + geom_line(aes(x = mar, y = allele, col = allele, alpha = prob),size = 5, show.legend = F)
  p
}
