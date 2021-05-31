#' Assign markers to linkage groups
#'
#' Identifies linkage groups of markers using the results of two-point
#' (pairwise) analysis and UPGMA method. Function adapted from MAPpoly package
#' written by Marcelo Mollinari.
#'
#' @param input.seq an object of class \code{mappoly.rf.matrix}
#'
#' @param expected.groups when available, inform the number of expected 
#' linkage groups (i.e. chromosomes) for the species
#'
#' @param inter if \code{TRUE} (default), plots a dendrogram highlighting the
#'    expected groups before continue
#'
#' @param comp.mat if \code{TRUE}, shows a comparison between the reference
#'     based and the linkage based grouping, if the sequence information is
#'     available (default = FALSE)
#'
#' @param verbose logical. If \code{TRUE} (default), current progress is shown;
#'     if \code{FALSE}, no output is produced
#'
#' @return Returns an object of class \code{group}, which is a list
#'     containing the following components:
#'     \item{data.name}{the referred dataset name}
#'     \item{hc.snp}{a list containing information related to 
#'     the UPGMA grouping method}
#'     \item{expected.groups}{the number of expected linkage groups}
#'     \item{groups.snp}{the groups to which each of the markers belong}
#'     \item{seq.vs.grouped.snp}{comparison between the genomic group information
#'     (when available) and the groups provided by \code{group_upgma}}
#'     \item{LOD}{minimum LOD Score to declare linkage.} 
#'     \item{max.rf}{maximum recombination fraction to declare linkage.}
#'     \item{twopt}{name of the object of class \code{rf.2ts}
#'     used as input, i.e., containing information used to assign
#'     markers to linkage groups.}
#'
#' @examples
#'
#' data("vcf_example_out")
#' twopts <- rf_2pts(vcf_example_out)
#' input.seq <- make_seq(twopts, "all")
#' lgs <- group_upgma(input.seq, expected.groups = 3, comp.mat=TRUE, inter = FALSE)
#' plot(lgs)
#'    
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @author Cristiane Taniguti \email{chtaniguti@usp.br}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @importFrom graphics abline pie
#' @importFrom stats as.dendrogram as.dist cutree hclust lm predict quantile rect.hclust
#' @importFrom dendextend color_branches
#' 
#' @export 
group_upgma <- function(input.seq, expected.groups = NULL,
                        inter = TRUE, comp.mat = FALSE, verbose = TRUE)
{
  ## checking for correct objects
  if(!any(is(input.seq,"sequence")))
    stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
  
  ## extracting data
  if(is(input.seq$data.name, "outcross") | is(input.seq$data.name, "f2"))
  {
    ## making a list with necessary information
    n.mrk <- length(input.seq$seq.num)
    LOD <- lapply(input.seq$twopt$analysis,
                  function(x, w){
                    m<-matrix(0, length(w), length(w))
                    for(i in 1:(length(w)-1)){
                      for(j in (i+1):length(w)){
                        z<-sort(c(w[i],w[j]))
                        m[j,i]<-m[i,j]<-x[z[1], z[2]]
                      }
                    }
                    return(m)
                  }, input.seq$seq.num
    )
    mat<-t(get_mat_rf_out(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
  } else {
    #if(input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like))
    #stop("You must estimate parameters before running 'rf_graph_table' ")
    ## making a list with necessary information
    n.mrk <- length(input.seq$seq.num) 
    LOD<-matrix(0, length(input.seq$seq.num), length(input.seq$seq.num))
    for(i in 1:(length(input.seq$seq.num)-1)){
      for(j in (i+1):length(input.seq$seq.num)){
        z<-sort(c(input.seq$seq.num[i],input.seq$seq.num[j]))
        LOD[j,i]<-LOD[i,j]<-input.seq$twopt$analysis[z[1], z[2]]
      }
    }
    mat<-t(get_mat_rf_in(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
  }
  
  diag(mat) <- 0
  mn<-input.seq$data.name$CHROM[input.seq$seq.num]
  mn[is.na(mn)]<-"NH"
  dimnames(mat)<-list(mn, mn)
  
  hc.snp<-hclust(as.dist(mat), method = "average")
  ANSWER <- "flag"
  if(interactive() && inter)
  {
    dend.snp <- as.dendrogram(hc.snp)
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
    {
      dend1 <- dendextend::color_branches(dend.snp, k = expected.groups)
      plot(dend1, leaflab = "none")
      if(is.null(expected.groups))
        expected.group <- as.numeric(readline("Enter the number of expected groups: "))
      z<-rect.hclust(hc.snp, k = expected.groups, border = "red")
      groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
      xy<-sapply(z, length)
      xt<-as.numeric(cumsum(xy)-ceiling(xy/2))
      yt<-.1
      points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
      text(x = xt, y = yt, labels = pmatch(xy, table(groups.snp, useNA = "ifany")), adj = .5)
      ANSWER <- readline("Enter 'Y/n' to proceed or update the number of expected groups: ")
      if(substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "no" | substr(ANSWER, 1, 1) == "N")
        stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
        expected.groups <- as.numeric(ANSWER)
    }
  }
  if(is.null(expected.groups))
    stop("Inform the 'expected.groups' or use 'inter = TRUE'")
  
  # Distribution of SNPs into linkage groups
  seq.vs.grouped.snp <- NULL
  if(all(unique(mn) == "NH") & comp.mat)
  {
    comp.mat <- FALSE
    seq.vs.grouped.snp <- NULL
    warning("There is no physical reference to generate a comparison matrix")
  }
  groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
  if(comp.mat){
    seq.vs.grouped.snp<-matrix(0, expected.groups, length(na.omit(unique(mn)))+1,
                               dimnames = list(1:expected.groups, c(na.omit(unique(mn)),"NH")))
    for(i in 1:expected.groups)
    {
      x<-table(names(which(groups.snp==i)))
      seq.vs.grouped.snp[i,names(x)]<-x
    }
    idtemp2 <- apply(seq.vs.grouped.snp, 1, which.max)
    seq.vs.grouped.snp<-cbind(seq.vs.grouped.snp[,unique(idtemp2)], seq.vs.grouped.snp[,"NH"])
    cnm<-colnames(seq.vs.grouped.snp)
    cnm[cnm==""]<-"NoChr"
    colnames(seq.vs.grouped.snp)<-cnm
  } else {
    seq.vs.grouped.snp <- NULL
  }
  names(groups.snp)<-input.seq$seq.num
  structure(list(data.name = input.seq$data.name, 
                 hc.snp = hc.snp, 
                 expected.groups = expected.groups,
                 n.groups = length(unique(groups.snp)),
                 groups = groups.snp, 
                 seq.num =input.seq$seq.num,
                 seq.vs.grouped.snp = seq.vs.grouped.snp, 
                 LOD = input.seq$LOD, 
                 max.rf = input.seq$max.rf,
                 twopt = input.seq$twopt),
            class = "group.upgma")
}

#' @export
print.group.upgma <- function(x, detailed = TRUE, ...) {
  cat("  This is an object of class 'group.upgma'\n  ------------------------------------------\n")
  ## criteria
  cat("  Criteria used to assign markers to groups:\n\n")
  cat("    - Number of markers:         ", length(x$groups), "\n")
  cat("    - Number of linkage groups:  ", length(unique(x$groups)), "\n")
  cat("    - Number of markers per linkage groups: \n")
  w<-data.frame(table(x$groups, useNA = "ifany"))
  colnames(w) = c("   group", "n.mrk")
  print (w, row.names = FALSE)
  cat("  ------------------------------------------\n")
  
  ## printing summary
  if(!is.null(x$seq.vs.grouped.snp)){
    print(x$seq.vs.grouped.snp)
    cat("  ------------------------------------------\n")
  }
}

#' @export
plot.group.upgma <- function(x, ...) {
  dend <- as.dendrogram(x$hc.snp)
  dend1 <- dendextend::color_branches(dend, k = x$expected.groups)
  plot(dend1, leaflab = "none")
  z<-rect.hclust(x$hc.snp, k = x$expected.groups, border = "red")
  xy<-sapply(z, length)
  xt<-as.numeric(cumsum(xy)-ceiling(xy/2))
  yt<-.1
  points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
  text(x = xt, y = yt, labels = pmatch(xy, table(x$groups, useNA = "ifany")), adj = .5)
}