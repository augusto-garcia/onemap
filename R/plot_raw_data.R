#########################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: plot_raw_data.R                                               ##
## Contains: plot.onemap, create_dataframe_for_plot_outcross,          ##
## plot_by_segreg_type                                                 ##
##                                                                     ##
## Written by Antonio Augusto Franco Garcia with minor modifications   ##
## by Marcelo Mollinari, Gabriel Rodrigues Alves Margarido and         ##
## Cristiane Taniguti                                                  ##
## copyright (c) 2015 Antonio Augusto Franco Garcia                    ##
##                                                                     ##
## First version: 2015/03/31                                           ##
## Last update: 2017/12/18                                             ##
## License: GNU General Public License version 3 or later              ##
##                                                                     ##
#########################################################################

globalVariables(c("ind", "variable", "value"))
globalVariables(c("Type", "segr.type"))
globalVariables(c("marker", "geno"))

##' Draw a graphic of raw data for any OneMap population
##'
##' Shows a heatmap (in ggplot2, a graphic of geom "tile") for raw data.
##' Lines correspond to markers and columns to individuals.
##' The function can plot a graph for all marker types, depending of the cross type (dominant/codominant markers, in all combinations).
##' The function receives a onemap object of class \code{onemap}, reads information
##' from genotypes from this object, converts it to a long dataframe format
##' using function melt() from package reshape2() or internal function create_dataframe_for_plot_outcross(), converts numbers from the object
##' to genetic notation (according to the cross type), then plots the graphic.
##' If there is more than 20 markers, removes y labels
##' For outcross populations, it can show all markers together, or it can split them according the segregation
##' pattern.
##'
##' @param x an object of class \code{onemap}, with data and additional information
##'
##' @param all a TRUE/FALSE option to indicate if results will be
##'     plotted together (if TRUE) or splitted based on their
##'     segregation pattern. Only used for outcross populations.
##'
##' @param ... currently ignored
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##' 
##' \dontrun{
##' data(example_bc) # Loads a fake backcross dataset installed with onemap
##' plot(example_bc) # This will show you the graph
##'
##' # You can store the graphic in an object, then save it with a number of properties
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(example_bc)
##' ggplot2::ggsave("MyRawData_bc.jpg", g, width=7, height=4, dpi=600)
##'
##' data(onemap_example_f2) # Loads a fake backcross dataset installed with onemap
##' plot(onemap_example_f2) # This will show you the graph
##'
##' # You can store the graphic in an object, then save it with a number of properties
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(onemap_example_f2)
##' ggplot2::ggsave("MyRawData_f2.jpg", g, width=7, height=4, dpi=600)
##'
##' data(example_out) # Loads a fake full-sib dataset installed with onemap
##' plot(example_out) # This will show you the graph for all markers
##' plot(example_out, all=FALSE) # This will show you the graph splitted for marker types
##'
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(example_out, all=FALSE)
##' ggplot2::ggsave("MyRawData_out.jpg", g, width=9, height=4, dpi=600)
##'}
##' @export
plot.onemap <- function(x, all=TRUE, ...) {
    # Creating the data frame
    if (is(x, "outcross")) {
        df.OM <- create_dataframe_for_plot_outcross(x)
        df.OM$geno <- as.numeric(as.character((df.OM$geno)))

        # Defining the label for genotypes
        temp <- which(df.OM[,3]==0)
        df.OM$geno[temp] <- "-"

        temp <- which(df.OM$Mk.type == "A.1")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ac"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ad"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "bc"
        df.OM$geno[temp][which(df.OM[temp,3]==4)] <- "bd"

        temp <- which(df.OM$Mk.type == "A.2")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ac"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "ba"
        df.OM$geno[temp][which(df.OM[temp,3]==4)] <- "bc"

        temp <- which(df.OM$Mk.type == "A.3")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ac"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "bc"
        df.OM$geno[temp][which(df.OM[temp,3]==4)] <- "b"

        temp <- which(df.OM$Mk.type == "A.4")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "b"
        df.OM$geno[temp][which(df.OM[temp,3]==4)] <- "o"

        temp <- which(df.OM$Mk.type == "B1.5")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "b"

        temp <- which(df.OM$Mk.type == "B2.6")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "b"

        temp <- which(df.OM$Mk.type == "B3.7")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==3)] <- "b"

        temp <- which(df.OM$Mk.type == "C.8")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "o"

        temp <- which(df.OM$Mk.type == "D1.9")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ac"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "bc"

        temp <- which(df.OM$Mk.type == "D1.10")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ab"

        temp <- which(df.OM$Mk.type == "D1.11")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "b"

        temp <- which(df.OM$Mk.type == "D1.12")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "a"

        temp <- which(df.OM$Mk.type == "D1.13")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "o"

        temp <- which(df.OM$Mk.type == "D2.14")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ac"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "bc"

        temp <- which(df.OM$Mk.type == "D2.15")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "ab"

        temp <- which(df.OM$Mk.type == "D2.16")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "b"

        temp <- which(df.OM$Mk.type == "D2.17")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "ab"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "a"

        temp <- which(df.OM$Mk.type == "D2.18")
        df.OM$geno[temp][which(df.OM[temp,3]==1)] <- "a"
        df.OM$geno[temp][which(df.OM[temp,3]==2)] <- "o"

        df.OM$geno <- as.factor(df.OM$geno)
    }
  else {
        df.OM <- as.matrix(x$geno)
        df.OM <- reshape2::melt(df.OM) # function from package reshape
        df.OM[is.na(df.OM)] <- 0 # To avoid problems with NAs
        df.OM <- data.frame("ind"=1:x$n.ind, "marker"= df.OM$Var2, "geno" = df.OM$value)
        df.OM$geno <- factor(df.OM$geno)
    }
    # Defining the label for genotypes
    if (is(x, "backcross")) {
        if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","2"))))
            labels.OM <- c("-","AA","AB")
        else if (all(levels(df.OM$geno)==c("1","2")))
            labels.OM <- c("AA","AB")
    } else if (is(x, "riself") || is(x, "risib")) {
        if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","3"))))
            labels.OM <- c("-","AA","BB")
        else if (all(levels(df.OM$geno)==c("1","3")))
            labels.OM <- c("AA","BB")
    } else if (is(x, "f2")) {
        if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","2","3","4","5"))))
            labels.OM <- c("-","AA","AB","BB","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("1","2","3","4","5"))))
            labels.OM <- c("AA","AB","BB","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","2","3"))))
            labels.OM <- c("-","AA","AB","BB")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("1","2","3"))))
            labels.OM <- c("AA","AB","BB")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","2","3","5"))))
            labels.OM <- c("-","AA","AB","BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("1","2","3","5"))))
            labels.OM <- c("AA","AB","BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("0","1","2","3","4"))))
            labels.OM <- c("-","AA","AB","BB","not BB")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("1","2","3","4"))))
            labels.OM <- c("AA","AB","BB","not BB")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("0","4","5"))))
            labels.OM <- c("-","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$geno)==c("4","5"))))
            labels.OM <- c("not BB","not AA")
    }
    # Plotting
    g <- ggplot(data=df.OM, aes(x=ind, y=marker, fill=factor(geno)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker")
    if (is(x, "outcross")) {
      if(length(which(df.OM$geno=="-")) != 0){
        g <- g + scale_fill_manual(name="Genotypes",
                                   values = c("black",'#e31a1c','#1f78b4','#6a3d9a','#33a02c','#ff7f00',
                                              '#b2df8a','#fb9a99','#fdbf6f','#a6cee3'))

      } else {
        g <- g + scale_fill_manual(name="Genotypes",
                                   values =c('#e31a1c','#1f78b4','#6a3d9a','#33a02c','#ff7f00',
                                             '#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#a6cee3'))
      }
      if (x$n.mar>20)
            g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

        if (all==TRUE) g
        else  g <- g + facet_wrap( ~Mk.type, ncol=1, scales = "free", strip.position = "right") +
            theme(axis.ticks.x =  element_blank(), axis.text.x = element_blank(),
                  strip.text.y = element_text(angle = 0))
    }
    else {
        if (is(x, "backcross") || is(x, "riself") || is(x, "risib")) {
            g <- g + scale_fill_manual(name="Genotype", labels=labels.OM,
                                       values=c("#F21A00","#3B9AB2","#EBCC2A"))
        } else if (is(x, "f2")) {
            g <- g + scale_fill_manual(name="Genotype", labels=labels.OM,
                                       values=c("#000000", "#ECCBAE", "#046C9A", "#D69C4E", "#85D4E3", "#74A089"))
        }
        if (x$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    return(g)
}
##'


##' Create a dataframe suitable for a ggplot2 graphic
##'
##' An internal function that prepares a dataframe suitable for
##' drawing a graphic of raw data using ggplot2, i. e., a data frame
##' with long format
##'
##' @param x an object of classes \code{onemap} and \code{outcross}, with data and additional information
##'
##' @return a dataframe
##'
create_dataframe_for_plot_outcross <- function(x) {
  # Markers of A type
  for (i in 1:4) {
    if (length(which(x$segr.type==paste("A.",i,sep="")))!=0) {
      F1.A <- as.matrix(data.frame(x$geno[,which(x$segr.type==paste("A.",i,sep=""))]))
      colnames(F1.A) <- colnames(x$geno)[which(x$segr.type==paste("A.",i,sep=""))]
      F1.A <- reshape2::melt(F1.A)
      F1.A <- data.frame("ind"=1:x$n.ind, "marker"= F1.A$Var2, "geno" = F1.A$value,
                         Mk.type=paste("A.",i,sep=""))
      F1.A$geno <- factor(F1.A$geno)
      if (!exists("F1.A.")) F1.A. <- data.frame(F1.A)
      else F1.A. <- rbind(F1.A.,F1.A)
    }
  }
  if (!exists("F1.A.")) F1.A. <- NULL
  # Markers of B type
  for (i in 1:3) {
    if (length(which(x$segr.type==paste("B",i,".",i+4,sep="")))!=0) {
      F1.B <- as.matrix(x$geno[,which(x$segr.type==paste("B",i,".",i+4,sep=""))])
      colnames(F1.B) <- colnames(x$geno)[which(x$segr.type==paste("B",i,".",i+4,sep=""))]
      F1.B <- reshape2::melt(F1.B)
      F1.B <- data.frame("ind"=1:x$n.ind, "marker" = F1.B$Var2, "geno" = F1.B$value,
                         Mk.type=paste("B",i,".",i+4,sep=""))
      F1.B$geno <- factor(F1.B$geno)
      if (!exists("F1.B.")) F1.B. <- F1.B
      else F1.B. <- rbind(F1.B.,F1.B)
    }
  }
  if (!exists("F1.B.")) F1.B. <- NULL
  # Markers of C type
  for (i in 8) {
    if (length(which(x$segr.type==paste("C.",i,sep="")))!=0) {
      F1.C <- as.matrix(x$geno[,which(x$segr.type==paste("C.",i,sep=""))])
      colnames(F1.C) <- colnames(x$geno)[which(x$segr.type==paste("C.",i,sep=""))]
      F1.C <- reshape2::melt(F1.C)
      F1.C <- data.frame("ind"=1:x$n.ind, "marker" = F1.C$Var2, "geno" = F1.C$value,
                         Mk.type=paste("C.",i,sep=""))
      F1.C$geno <- factor(F1.C$geno)
      if (!exists("F1.C.")) F1.C. <- data.frame(F1.C)
      else F1.C. <- rbind(F1.C.,F1.C)
    }
  }
  if (!exists("F1.C.")) F1.C. <- NULL
  # Markers of D1 type
  for (i in 9:13) {
    if (length(which(x$segr.type==paste("D1.",i,sep="")))!=0) {
      F1.D1 <- as.matrix(x$geno[,which(x$segr.type==paste("D1.",i,sep=""))])
      colnames(F1.D1) <- colnames(x$geno)[which(x$segr.type==paste("D1.",i,sep=""))]
      F1.D1 <- reshape2::melt(F1.D1)
      F1.D1 <- data.frame("ind"=1:x$n.ind, "marker" = F1.D1$Var2, "geno" = F1.D1$value,
                          Mk.type=paste("D1.",i,sep=""))
      F1.D1$geno <- factor(F1.D1$geno)
      if (!exists("F1.D1.")) F1.D1. <- data.frame(F1.D1)
      else F1.D1. <- rbind(F1.D1.,F1.D1)
    }
  }
  if (!exists("F1.D1.")) F1.D1. <- NULL

  # Markers of D2 type
  for (i in 14:18) {
    if (length(which(x$segr.type==paste("D2.",i,sep="")))!=0) {
      F1.D2 <- as.matrix(x$geno[,which(x$segr.type==paste("D2.",i,sep=""))])
      colnames(F1.D2) <- colnames(x$geno)[which(x$segr.type==paste("D2.",i,sep=""))]
      F1.D2 <- reshape2::melt(F1.D2)
      F1.D2 <- data.frame("ind"=1:x$n.ind, "marker" = F1.D2$Var2, "geno" = F1.D2$value,
                          Mk.type=paste("D2.",i,sep=""))
      F1.D2$geno <- factor(F1.D2$geno)
      if (!exists("F1.D2.")) F1.D2. <- data.frame(F1.D2)
      else F1.D2. <- rbind(F1.D2.,F1.D2)
    }
  }
  if (!exists("F1.D2.")) F1.D2. <- NULL
  # Defining classes and combining
    # It was not working if any class was missing; now it is
    if (!is.null(F1.A.)) class(F1.A.) <- "data.frame"
    if (!is.null(F1.B.)) class(F1.B.) <- "data.frame"
    if (!is.null(F1.C.)) class(F1.C.) <- "data.frame"
    if (!is.null(F1.D1.)) class(F1.D1.) <- "data.frame"
    if (!is.null(F1.D2.)) class(F1.D2.) <- "data.frame"
    return(rbind(F1.A.,F1.B.,F1.C.,F1.D1.,F1.D2.))
}
##'


##' Draw a graphic showing the number of markers of each segregation pattern.
##'
##' The function receives an object of class \code{onemap}.
##' For outcrossing populations, it can show detailed information (all 18 possible categories),
##' or a simplified version.
##'
##' @param x an object of class \code{onemap}
##' @param subcateg a TRUE/FALSE option to indicate if results will be plotted showing
##' all possible categories (only for outcrossing populations)
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##' data(example_out) #Outcrossing data
##' plot_by_segreg_type(example_out)
##' plot_by_segreg_type(example_out, subcateg=FALSE)
##'
##' data(example_bc)
##' plot_by_segreg_type(example_bc)
##'
##' data(mapmaker_example_f2)
##' plot_by_segreg_type(mapmaker_example_f2)
##'
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' # data(example_out) #Outcrossing data
##' # g <- plot_by_segreg_type(example_out)
##' # ggplot2::ggsave("SegregationTypes.jpg", g, width=7, height=4, dpi=600)
##'
##' @export
plot_by_segreg_type <- function(x, subcateg=TRUE) {
  # Create a dataframe, indicating the category and subcategory
  df <- data.frame(segr.type=factor(x$segr.type),Type=999)
  t <- c("A","B","C","D1","D2", "A.H", "A.H.B", "C.A", "D.B")
  for (i in t){
    if (length(grep(i, x$segr.type))>0)
      df[grep(i, x$segr.type),]$Type <- i
  }
  # Plot the graphic, with option for subcategory for outcross
  if (subcateg==TRUE){
    g <- ggplot(data=df, aes(x=Type, fill=segr.type))
    g <- g + geom_bar()
    g <- g + xlab("Segregation Type") + ylab("Count")
  }
  else {
    g <- ggplot(data=df, aes(x=Type, fill="orange"))
    g <- g + geom_bar()
    g <- g + xlab("Segregation Type") + ylab("Count") + theme(legend.position="none")
  }
  return(g)
}

#####
