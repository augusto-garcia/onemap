#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: plot_raw_data.R                                               #
# Contains: plot.bc.onemap, plot.riself.onemap, plot.risib.onemap,    #
# plot.f2.onemap, create_dataframe_for_plot_outcross, plot.outcross,  #
# plot_by_segreg_type                                                 #
#                                                                     #
# Written by Antonio Augusto Franco Garcia                            #
# copyright (c) 2015 Antonio Augusto Franco Garcia                    #
#                                                                     #
# First version: 2015/03/31                                           #
# Last update: 2015/07/25                                             #
# License: GNU General Public License version 3 or later              #
#                                                                     #
#######################################################################

##' Draw a graphic of raw data for a backcross population
##'
##' Shows a heatmap (in ggplot2, a graphic of geom "tile") for raw data.
##' Lines correspond to markers and columns for individuals.
##' The function receives a onemap object of class bc.onemap, reads information
##' from genotypes from this object, convert it to a long dataframe format
##' using function melt() from package reshape2(), converts numbers from the object
##' to genetic notation (AA, AB, -), then plot the graphic.
##' If there is more than 20 markers, removes y labels
##'
##' @param x an object of class bc.onemap, with data and additional information
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
##' ggsave("MyRawData.jpg", g, width=7, height=4, dpi=600)
##'
##' @export
plot.bc.onemap <- function(x) {
    # Creating the data frame
    df.BC <- data.frame(x$geno)
    df.BC <- reshape2::melt(df.BC) # function from package reshape
    df.BC[is.na(df.BC)] <- 0 # To avoid problems with NAs
    df.BC <- cbind(ind=1:x$n.ind, df.BC)
    df.BC$value <- factor(df.BC$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.BC$value)==c("0","1","2")))) labels.bc <- c("-","AA","AB")
    else if (all(levels(df.BC$value)==c("1","2"))) labels.bc <- c("AA","AB")
    # Plotting
    g <- ggplot(data=df.BC, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype",labels=labels.bc,
                          values=c("#F21A00","#3B9AB2","#EBCC2A"))
    if (x$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    return(g)
}
##'

##' Draw a heatmap graphic of raw data for a RIL population (made by selfing)
##'
##' Draw a graphic of raw data for a RIL population (made by selfing), using ggplot2.
##' Lines correspond to markers and columns for individuals.
##' The graphic is a "heatmap", whose geom in ggplot2 is "tile".
##' The function receives a onemap object of class riself.onemap, reads information
##' from genotypes from this object, convert it to a long dataframe format
##' using function melt() from package reshape2(), converts numbers from the object
##' to genetic notation (AA, BB, -), then plot the graphic.
##' If there is more than 20 markers, removes y labels
##'
##' @param x an object of class riself.onemap, with data and additional information
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##' # Please, see examples for plot.bc.onemap; the only difference is the dataset
##'
##' @export
plot.riself.onemap <- function(x) {
    # Creating the data frame
    df.RIL <- data.frame(x$geno)
    df.RIL <- reshape2::melt(df.RIL) # function from package reshape
    df.RIL[is.na(df.RIL)] <- 0 # To avoid problems with NAs
    df.RIL <- cbind(ind=1:x$n.ind, df.RIL)
    df.RIL$value <- factor(df.RIL$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.RIL$value)==c("0","1","2")))) labels.ril <- c("-","AA","BB")
    else if (all(levels(df.RIL$value)==c("1","2"))) labels.ril <- c("AA","BB")
    # Plotting
    g <- ggplot(data=df.RIL, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",labels=labels.ril, values=c("#F21A00","#3B9AB2","#EBCC2A"))
    if (x$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    return(g)
}
##'

##' Draw a heatmap graphic of raw data for a RIL population (made by sibing)
##'
##' Draw a graphic of raw data for a RIL population (made by sibing), using ggplot2.
##' In fact, the graphic for the raw data will have the same aspect of the ones
##' for RILs made by selfing. Therefore, this function will only call
##' plot.riself.onemap()
##'
##' @param x an object of class risib.onemap, with data and additional information
##'
##' @return a ggplot graphic
##'
##' @examples
##' # Please, see examples for plot.bc.onemap; the only difference is the dataset
##'
##' @export
plot.risib.onemap <- function(x) {
    plot.riself.onemap(x)
}
##'

##' Draw a graphic of raw data for an F2 population
##'
##' Draw a graphic of raw data for a f2 population, using ggplot2.
##' Lines correspond to markers and columns for individuals.
##' The function can plot graph for dominant/codominant markers, in all combinations.
##' The graphic is a "heatmap", whose geom in ggplot2 is "tile".
##' The function receives a onemap object of class f2.onemap, reads information
##' from genotypes from this object, convert it to a long dataframe format
##' using function melt() from package reshape2(), converts numbers from the object
##' to genetic notation (AA, AB, BB, not AA, not BB, -), then plot the graphic.
##' If there is more than 20 markers, removes y labels
##'
##' @param x an object of class f2.onemap, with data and additional information
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##' data(fake.f2.onemap) # Loads a fake backcross dataset installed with onemap
##' plot(fake.f2.onemap) # This will show you the graph
##'
##' # You can store the graphic in an object, then save it with a number of properties
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(fake.f2.onemap)
##' ggsave("MyRawData.jpg", g, width=7, height=4, dpi=600)
##'
##' @export
plot.f2.onemap <- function(x) {
    # Creating the data frame
    df.F2 <- data.frame(x$geno.mmk[[1]])
    df.F2 <- reshape2::melt(df.F2) # function from package reshape
    df.F2[is.na(df.F2)] <- 0 # To avoid problems with NAs
    df.F2 <- cbind(ind=1:x$n.ind, df.F2)
    df.F2$value <- factor(df.F2$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4","5")))) labels.f2 <-
        c("-","AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4","5")))) labels.f2 <-     c("AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3")))) labels.f2 <- c("-","AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3")))) labels.f2 <- c("AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","5")))) labels.f2 <- c("-","AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","5")))) labels.f2 <- c("AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4")))) labels.f2 <- c("-","AA","AB","BB","not BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4")))) labels.f2 <- c("AA","AB","BB","not BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","4","5")))) labels.f2 <-
        c("-","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("4","5")))) labels.f2 <-
        c("not BB","not AA")
    # Plotting
    g <- ggplot(data=df.F2, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype", labels=labels.f2, values=c("#000000", "#ECCBAE", "#046C9A", "#D69C4E", "#85D4E3", "#74A089"))
    if (x$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    return(g)
}
##'


##' Create a dataframe suitable for a ggplot2 graphic
##'
##' An internal function that prepares a dataframe suitable for
##' drawing a graphic of raw data using ggplot2, i. e., a data frame
##' with long format
##'
##' @param x an object of class outcross, with data and additional information
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
        }
    }
    class(F1.A.) <- "data.frame"
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
        }
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
        }
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
        }
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
        }
    }
    # Defining classes and combining
    class(F1.A.) <- "data.frame"
    class(F1.B.) <- "data.frame"
    class(F1.C.) <- "data.frame"
    class(F1.D1.) <- "data.frame"
    class(F1.D2.) <- "data.frame"
    return(rbind(F1.A.,F1.B.,F1.C.,F1.D1.,F1.D2.))
}
##'


##' Draw a heatmap graphic for raw data in an outcross population
##'
##' Draw a graphic of raw data for an outcross population, using ggplot2.
##' Lines correspond to markers and columns for individuals.
##' The function can plot a graph for all types of markers.
##' The graphic is a "heatmap", whose geom in ggplot2 is "tile".
##' The function receives a onemap object of class outcross, reads information
##' from genotypes from this object, convert it to a long dataframe format
##' using onemap internal function create_dataframe_for_plot_outcross(),
##' then plot the graphic.
##' If there is more than 20 markers, removes y labels
##' It can show all markers together, or it can split them according the segregation
##' pattern.
##'
##' @param x an object of class outcross, with data and additional information
##' @param all a TRUE/FALSE option to indicate if results will be plotted together (if TRUE)
##' or splitted based on their segregation pattern
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##' data(example.out) # Loads a fake backcross dataset installed with onemap
##' plot(example.out) # This will show you the graph for all markers
##' plot(example.out, all=FALSE) # This will show you the graph splitted for marker types
##'
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(example.out, all=FALSE)
##' ggsave("MyRawData.jpg", g, width=9, height=4, dpi=600)
##'
##' @export
plot.outcross <- function(x, all=TRUE) {
    df <- create_dataframe_for_plot_outcross(x)
    g <- ggplot(data=df, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype",
                          values=c("#00A08A", "#5BBCD6",  "#F2AD00", "#F98400", "#FF0000"))
    if (length(levels(x$variable))>20) g <- g +
                          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (all==TRUE) g
    else g <- g + facet_grid( . ~ Mk.type) + theme(axis.text.x=element_text(size=5))
    return(g)
}
##


##' Draw a graphic showing the number of markers of each segregation pattern.
##' 
##' The function receives a onemap object of class outcross, f2.onemap, bc.onemap,
##' risib.onemap or riself.onemap.
##' For outcrossing populations, it can show detailed information (all 18 possible categories),
##' or a simplified version.
##'
##' @param x an object of class outcross, f2.onemap, bc.onemap, risib.onemap or riself.onemap
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
    if (subcateg==TRUE & class(x)=="outcross"){
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
