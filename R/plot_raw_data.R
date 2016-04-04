#########################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: plot_raw_data.R                                               ##
## Contains: plot.onemap, create_dataframe_for_plot_outcross,          ##
## plot_by_segreg_type                                                 ##
##                                                                     ##
## Written by Antonio Augusto Franco Garcia with minor modifications   ##
## by Marcelo Mollinari and Gabriel Rodrigues Alves Margarido          ##
## copyright (c) 2015 Antonio Augusto Franco Garcia                    ##
##                                                                     ##
## First version: 2015/03/31                                           ##
## Last update: 2016/01/15                                             ##
## License: GNU General Public License version 3 or later              ##
##                                                                     ##
#########################################################################

globalVariables(c("ind", "variable", "value"))
globalVariables(c("Type", "segr.type"))

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
##' data(fake.bc.onemap) # Loads a fake backcross dataset installed with onemap
##' plot(fake.bc.onemap) # This will show you the graph
##'
##' # You can store the graphic in an object, then save it with a number of properties
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(fake.bc.onemap)
##' ggsave("MyRawData_bc.jpg", g, width=7, height=4, dpi=600)
##'
##' data(fake.f2.onemap) # Loads a fake backcross dataset installed with onemap
##' plot(fake.f2.onemap) # This will show you the graph
##'
##' # You can store the graphic in an object, then save it with a number of properties
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(fake.f2.onemap)
##' ggsave("MyRawData_f2.jpg", g, width=7, height=4, dpi=600)
##'
##' data(example.out) # Loads a fake full-sib dataset installed with onemap
##' plot(example.out) # This will show you the graph for all markers
##' plot(example.out, all=FALSE) # This will show you the graph splitted for marker types
##'
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(example.out, all=FALSE)
##' ggsave("MyRawData_out.jpg", g, width=9, height=4, dpi=600)
##'
##' @export
plot.onemap <- function(x, all=TRUE, ...) {
    # Creating the data frame
    if (is(x, "outcross")) {
        df.OM <- create_dataframe_for_plot_outcross(x)
    }
    else {
        df.OM <- data.frame(x$geno)
        df.OM <- reshape2::melt(df.OM) # function from package reshape
        df.OM[is.na(df.OM)] <- 0 # To avoid problems with NAs
        df.OM <- cbind(ind=1:x$n.ind, df.OM)
        df.OM$value <- factor(df.OM$value)
    }
    # Defining the label for genotypes
    if (is(x, "backcross")) {
        if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2"))))
            labels.OM <- c("-","AA","AB")
        else if (all(levels(df.OM$value)==c("1","2")))
            labels.OM <- c("AA","AB")
    } else if (is(x, "riself") || is(x, "risib")) {
        if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2"))))
            labels.OM <- c("-","AA","BB")
        else if (all(levels(df.OM$value)==c("1","2")))
            labels.OM <- c("AA","BB")
    } else if (is(x, "f2")) {
        if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2","3","4","5"))))
            labels.OM <- c("-","AA","AB","BB","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$value)==c("1","2","3","4","5"))))
            labels.OM <- c("AA","AB","BB","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2","3"))))
            labels.OM <- c("-","AA","AB","BB")
        else if (suppressWarnings(all(levels(df.OM$value)==c("1","2","3"))))
            labels.OM <- c("AA","AB","BB")
        else if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2","3","5"))))
            labels.OM <- c("-","AA","AB","BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$value)==c("1","2","3","5"))))
            labels.OM <- c("AA","AB","BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$value)==c("0","1","2","3","4"))))
            labels.OM <- c("-","AA","AB","BB","not BB")
        else if (suppressWarnings(all(levels(df.OM$value)==c("1","2","3","4"))))
            labels.OM <- c("AA","AB","BB","not BB")
        else if (suppressWarnings(all(levels(df.OM$value)==c("0","4","5"))))
            labels.OM <- c("-","not BB","not AA")
        else if (suppressWarnings(all(levels(df.OM$value)==c("4","5"))))
            labels.OM <- c("not BB","not AA")
    }
    # Plotting
    g <- ggplot(data=df.OM, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker")
    if (is(x, "outcross")) {
        g <- g + scale_fill_manual(name="Genotype",
                                   values=c("#00A08A", "#5BBCD6",  "#F2AD00", "#F98400", "#FF0000"))
        if (x$n.mar>20)
            g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        if (all==TRUE) g
        else g <- g + facet_grid( . ~ Mk.type) + theme(axis.text.x=element_text(size=5))
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
    if (length(which(x$segr.type==paste("A.",i,sep=""))!=0)) {
        F1.A <- data.frame(x$geno[,which(x$segr.type==paste("A.",i,sep=""))])
        colnames(F1.A) <- colnames(x$geno)[which(x$segr.type==paste("A.",i,sep=""))]
        F1.A <- reshape2::melt(F1.A)
        F1.A <- cbind(ind=1:x$n.ind, F1.A, Mk.type=paste("A.",i,sep=""))
        F1.A$value <- factor(F1.A$value)
        if (!exists("F1.A.")) F1.A. <- data.frame(F1.A)
        else F1.A. <- rbind(F1.A.,F1.A)
        } else F1.A. <- NULL
    }
    # Markers of B type
    for (i in 1:3) {
    if (length(which(x$segr.type==paste("B",i,".",i+4,sep=""))!=0)) {
        F1.B <- data.frame(x$geno[,which(x$segr.type==paste("B",i,".",i+4,sep=""))])
        colnames(F1.B) <- colnames(x$geno)[which(x$segr.type==paste("B",i,".",i+4,sep=""))]
        F1.B <- reshape2::melt(F1.B)
        F1.B <- cbind(ind=1:x$n.ind, F1.B, Mk.type=paste("B",i,".",i+4,sep=""))
        F1.B$value <- factor(F1.B$value)
        if (!exists("F1.B.")) F1.B. <- data.frame(F1.B)
        else F1.B. <- rbind(F1.B.,F1.B)
        } else F1.B. <- NULL
    }
    # Markers of C type
    for (i in 8) {
    if (length(which(x$segr.type==paste("C.",i,sep=""))!=0)) {
        F1.C <- data.frame(x$geno[,which(x$segr.type==paste("C.",i,sep=""))])
        colnames(F1.C) <- colnames(x$geno)[which(x$segr.type==paste("C.",i,sep=""))]
        F1.C <- reshape2::melt(F1.C)
        F1.C <- cbind(ind=1:x$n.ind, F1.C, Mk.type=paste("C.",i,sep=""))
        F1.C$value <- factor(F1.C$value)
        if (!exists("F1.C.")) F1.C. <- data.frame(F1.C)
        else F1.C. <- rbind(F1.C.,F1.C)
        } else F1.C. <- NULL
    }
    # Markers of D1 type
    for (i in 9:13) {
    if (length(which(x$segr.type==paste("D1.",i,sep=""))!=0)) {
        F1.D1 <- data.frame(x$geno[,which(x$segr.type==paste("D1.",i,sep=""))])
        colnames(F1.D1) <- colnames(x$geno)[which(x$segr.type==paste("D1.",i,sep=""))]
        F1.D1 <- reshape2::melt(F1.D1)
        F1.D1 <- cbind(ind=1:x$n.ind, F1.D1, Mk.type=paste("D1.",i,sep=""))
        F1.D1$value <- factor(F1.D1$value)
        if (!exists("F1.D1.")) F1.D1. <- data.frame(F1.D1)
        else F1.D1. <- rbind(F1.D1.,F1.D1)
        } else F1.D1. <- NULL
    }
    # Markers of D2 type
    for (i in 14:18) {
    if (length(which(x$segr.type==paste("D2.",i,sep=""))!=0)) {
        F1.D2 <- data.frame(x$geno[,which(x$segr.type==paste("D2.",i,sep=""))])
        colnames(F1.D2) <- colnames(x$geno)[which(x$segr.type==paste("D2.",i,sep=""))]
        F1.D2 <- reshape2::melt(F1.D2)
        F1.D2 <- cbind(ind=1:x$n.ind, F1.D2, Mk.type=paste("D2.",i,sep=""))
        F1.D2$value <- factor(F1.D2$value)
        if (!exists("F1.D2.")) F1.D2. <- data.frame(F1.D2)
        else F1.D2. <- rbind(F1.D2.,F1.D2)
        } else F1.D2. <- NULL
    }
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
##' data(example.out) #Outcrossing data
##' plot_by_segreg_type(example.out)
##' plot_by_segreg_type(example.out, subcateg=FALSE)
##'
##' data(fake.bc.onemap)
##' plot_by_segreg_type(fake.bc.onemap)
##'
##' data(fake.f2.onemap)
##' plot_by_segreg_type(fake.f2.onemap)
##'
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' data(example.out) #Outcrossing data
##' g <- plot_by_segreg_type(example.out)
##' ggsave("SegregationTypes.jpg", g, width=7, height=4, dpi=600)
##'
##' @export
plot_by_segreg_type <- function(x, subcateg=TRUE) {
    # Create a dataframe, indicating the category and subcategory
    df <- data.frame(segr.type=factor(x$segr.type),Type=999)
    t <- c("A","B","C","D1","D2")
    for (i in t){
        if (length(grep(i, x$segr.type))>0)
            df[grep(i, x$segr.type),]$Type <- i
    }
    # Plot the graphic, with option for subcategory for outcross
    if (is(x, "outcross") && subcateg==TRUE){
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
