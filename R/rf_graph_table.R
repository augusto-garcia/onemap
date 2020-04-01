#####################################################################################
##                                                                                 ##
## Package: onemap                                                                 ##
##                                                                                 ##
## File: rf_graph_table.R                                                          ##
## Contains: rf_graph_table                                                        ##
##                                                                                 ##
## Written by Marcelo Mollinari                                                    ##
## Upgraded to ggplot2/plotly by Rodrigo Amadeu and Cristiane Taniguti             ##
## copyright (c) 2018, Marcelo Mollinari and Rodrigo Amadeu and Cristiane Taniguti ##
##                                                                                 ##
## First version: 2009/05/03                                                       ##
## Last update: 2018/05/15                                                         ##
## Description was modified by Augusto Garcia on 2015/07/25                        ##
## Upgrade to ggplot2/plotly by Rodrigo Amadeu on 2018/02/15                       ##
## License: GNU General Public License version 2 (June, 1991) or later             ##
##                                                                                 ##
#####################################################################################

globalVariables(c("x", "y", "x.type", "rf"))
globalVariables(c("y.type", "x.missing"))
globalVariables(c("y.missing", "LOD.CC"))
globalVariables(c("LOD.CR", "LOD.RC", "LOD.RR"))

##' Plots pairwise recombination fractions and LOD Scores in a heatmap
##'
##' Plots a matrix of pairwise recombination fraction or
##' LOD Scores using a color scale. Any value of the
##' matrix can be easily accessed using an interactive plotly-html interface,
##' helping users to check for possible problems.
##'
##' The color scale varies from red (small distances or big LODs) to purple.
##' When hover on a cell, a dialog box is displayed with some information
##' about corresponding markers for that cell (line (y) \eqn{\times} column (x)). They are:
##' \eqn{i}) the name of the markers; \eqn{ii}) the number of
##' the markers on the data set; \eqn{iii}) the segregation types; \eqn{iv})
##' the recombination fraction between the markers and \eqn{v}) the LOD-Score
##' for each possible linkage phase calculated via two-point analysis. For
##' neighbor markers, the multipoint recombination fraction is printed;
##' otherwise, the two-point recombination fraction is printed. For markers of
##' type \code{D1} and \code{D2}, it is impossible to calculate recombination
##' fraction via two-point analysis and, therefore, the corresponding cell will
##' be empty (white color). For cells on the diagonal of the matrix, the name, the number and
##' the type of the marker are printed, as well as the percentage of missing
##' data for that marker.
##'
##'@import ggplot2 
##'
##' @param input.seq an object of class \code{sequence} with a predefined
##' order.
##' @param graph.LOD logical. If \code{TRUE}, displays the LOD heatmap, otherwise,
##' displays the recombination fraction heatmap.
##' @param main character. The title of the plot.
##' @param inter logical. If \code{TRUE}, an interactive HTML graphic is plotted.
##' Otherwise, a default graphic device is used.
##' @param html.file character naming the html file with iterative graphic.
##' @param mrk.axis character, "names" to display marker names in the axis, "numbers" to display
##' marker numbers and "none" to display axis free of labels.
##' @param lab.xy character vector with length 2, first component is the label of x axis and second of the y axis.
##' @param n.colors integer. Number of colors in the pallete.
##' 
##' @author Rodrigo Amadeu, \email{rramadeu@@gmail.com}
##' @keywords utilities
##' @examples
##'
##'\dontrun{
##' ##outcross example
##'   data(onemap_example_out)
##'   twopt <- rf_2pts(onemap_example_out)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'   LG1 <- make_seq(groups,1)
##'   LG1.rcd <- rcd(LG1)
##'   rf_graph_table(LG1.rcd, inter=FALSE)
##'
##'   ##Now, using interactive plotly
##'   rf_graph_table(LG1.rcd, inter=TRUE, html.file= "LG1.rcd.html")
##'
##'   ##F2 example
##'   data(onemap_example_f2)
##'   twopt <- rf_2pts(onemap_example_f2)
##'   all_mark <- make_seq(twopt,"all")
##'   groups <- group(all_mark)
##'
##'   ##"pre-allocate" an empty list of length groups$n.groups (3, in this case)
##'   maps.list<-vector("list", groups$n.groups)
##'
##'   for(i in 1:groups$n.groups){
##'     ##create linkage group i
##'     LG.cur <- make_seq(groups,i)
##'     ##ordering
##'     map.cur<-order_seq(LG.cur, subset.search = "sample")
##'     ##assign the map of the i-th group to the maps.list
##'     maps.list[[i]]<-make_seq(map.cur, "force")
##'   }
##'   ##Plot LOD/recombination fraction matrices for each group
##'   require(gridExtra)
##'   plot1 <- rf_graph_table(maps.list[[1]], main="Group 1",inter=FALSE)
##'   plot2 <- rf_graph_table(maps.list[[2]], main="Group 2",inter=FALSE)
##'   plot3 <- rf_graph_table(maps.list[[3]], main="Group 3",inter=FALSE)
##'   grid.arrange(plot1, plot2, plot3, nrow=3)
##' }
##'@export

rf_graph_table <- function(input.seq,
                             graph.LOD=FALSE,
                             main=NULL,
                             inter=FALSE,
                             html.file = NULL,
                             mrk.axis="numbers",
                             lab.xy=NULL,
                             n.colors=4){

    ## checking for correct objects
    if(!is(input.seq,"sequence"))
        stop(deparse(substitute(input.seq))," is not an object of class 'sequence'")
    if(!(mrk.axis=="names" | mrk.axis=="numbers" | mrk.axis=="none"))
      stop("This mrk.axis argument is not defined, choose 'names', 'numbers' or 'none'")

    ## extracting data
    if(is(get(input.seq$data.name), "outcross") | is(get(input.seq$data.name), "f2"))
    {
        if(input.seq$seq.phases[1] == -1 || input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like))
            stop("You must estimate parameters before running 'rf_graph_table' ")
        ## making a list with necessary information
        n.mrk <- length(input.seq$seq.num)
        LOD <- lapply(get(input.seq$twopt)$analysis,
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
    }else
    {
        if(input.seq$seq.rf[1] == -1 || is.null(input.seq$seq.like))
            stop("You must estimate parameters before running 'rf_graph_table' ")
        ## making a list with necessary information
        n.mrk <- length(input.seq$seq.num)

        LOD<-matrix(0, length(input.seq$seq.num), length(input.seq$seq.num))
        for(i in 1:(length(input.seq$seq.num)-1)){
            for(j in (i+1):length(input.seq$seq.num)){
                z<-sort(c(input.seq$seq.num[i],input.seq$seq.num[j]))
                LOD[j,i]<-LOD[i,j]<-get(input.seq$twopt)$analysis[z[1], z[2]]
            }
        }
        mat<-t(get_mat_rf_in(input.seq, LOD=TRUE,  max.rf = 0.501, min.LOD = -0.1))
    }

    ##Scaling the LODs to print them properly
    ## range.LOD<-range(as.dist(t(mat)), na.rm=TRUE)
    ## range.rf<-range(as.dist(mat), na.rm=TRUE)
    mat[row(mat) > col(mat) & mat > 0.5] <- 0.5 ## if there are recombinations greater than 0.5 (for numerical convergence problems), assuming 0.5
    mat[row(mat) < col(mat)][mat[row(mat) < col(mat)] < 10E-2]<-10E-2
    diag(mat)<-NA

    ##Write multipoint estimates
    for (i in 1:(n.mrk-1)){
        mat[i+1,i] <- input.seq$seq.rf[i]
    }
    
    colnames(mat) <- rownames(mat)<- colnames(get(input.seq$data.name, pos=1)$geno)[input.seq$seq.num]
    
    if (mrk.axis == "numbers")
      colnames(mat) <- rownames(mat)<- input.seq$seq.num
    
    ##Write NAs in two-point recombination fractions between markers of type D1 and D2
    types <- get(input.seq$data.name, pos=1)$segr.type.num[input.seq$seq.num]
    which.D1D2<-outer((substr(types, 1,2)== 6),(substr(types, 1,2)==7))
    which.D1D2<-which.D1D2+t(which.D1D2)
    which.D1D2[which.D1D2==1]<-NA
    which.D1D2[which.D1D2==0]<-1
    diag.si<-rbind(1:(ncol(which.D1D2)-1),2:ncol(which.D1D2))
    for(i in 1:(ncol(which.D1D2)-1)) which.D1D2[diag.si[1,i],diag.si[2,i]] <- which.D1D2[diag.si[2,i],diag.si[1,i]] <- 1
    mat<-mat*which.D1D2
    missing<-100*apply(get(input.seq$data.name, pos=1)$geno[,input.seq$seq.num],2, function(x) sum(x==0))/get(input.seq$data.name, pos=1)$n.ind

    ## Building the data.frame to plot
    mat.LOD <- mat.rf <- mat
    mat.LOD[lower.tri(mat.LOD)] <- t(mat.LOD)[lower.tri(mat.LOD)]
    mat.rf[upper.tri(mat.rf)] <- t(mat.rf)[upper.tri(mat.LOD)]

    if(is(get(input.seq$data.name), "outcross") | is(get(input.seq$data.name), "f2")){
        colnames(LOD$CC) <- rownames(LOD$CC) <- colnames(mat.rf)
        colnames(LOD$CR) <- rownames(LOD$CR) <- colnames(mat.rf)
        colnames(LOD$RC) <- rownames(LOD$RC) <- colnames(mat.rf)
        colnames(LOD$RR) <- rownames(LOD$RR) <- colnames(mat.rf)

        ## Merging all the matrices into one df
        df.graph <- Reduce(function(x, y) merge(x, y, all=TRUE),
                           list(reshape2::melt(round(mat.rf,2), value.name="rf"),
                                reshape2::melt(round(mat.LOD,2), value.name="LOD"),
                                reshape2::melt(round(LOD$CC,2), value.name="CC"),
                                reshape2::melt(round(LOD$CR,2), value.name="CR"),
                                reshape2::melt(round(LOD$RC,2), value.name="RC"),
                                reshape2::melt(round(LOD$RR,2), value.name="RR")))

        colnames(df.graph)[5:8] <- paste0("LOD.",c("CC","CR","RC","RR"))

        
        
    }else{
        df.graph <- merge(reshape2::melt(round(mat.rf,2), value.name="rf"),
                          reshape2::melt(round(mat.LOD,2), value.name="LOD"))
    }
    
    
    colnames(df.graph)[c(1,2)] <- c("x", "y")
    
    if(mrk.axis=="numbers"){
      df.graph$x <- factor(df.graph$x, levels = as.character(input.seq$seq.num))
      df.graph$y <- factor(df.graph$y, levels = as.character(input.seq$seq.num))
    }
    
    missing <- paste0(round(missing,2),"%")

    mrk.type.x <- data.frame(x=colnames(mat.rf),x.type=types)
    mrk.type.y <- data.frame(y=colnames(mat.rf),y.type=types)
    missing.x <- data.frame(x=colnames(mat.rf),x.missing=missing)
    missing.y <- data.frame(y=colnames(mat.rf),y.missing=missing)

    df.graph <- Reduce(function(x, y) merge(x, y, all=TRUE),
                       list(df.graph,
                            mrk.type.x,
                            mrk.type.y,
                            missing.x,
                            missing.y))

    ## Within the df.graph dataframe we plot based on the arguments and data.type
    ## Additional aesthetical (aes) arguments are to be passed to the mouse hover in the interactive plot.
    ## ggplot() just depends on the 'x', 'y', and 'fill' aes arguments
    
    ## If outcross:
    if(is(get(input.seq$data.name), "outcross") | is(get(input.seq$data.name), "f2")){
        if(graph.LOD!=TRUE){
            p <- ggplot(aes(x, y, x.type = x.type, y.type = y.type, x.missing = x.missing, y.missing = y.missing, fill = rf, LOD.CC=LOD.CC, LOD.CR=LOD.CR, LOD.RC=LOD.RC, LOD.RR=LOD.RR), data=df.graph) +
                geom_tile() +
                scale_fill_gradientn(colours = grDevices::rainbow(n.colors), na.value = "white") +
                theme(axis.text.x=element_text(angle=90, hjust=1))
        }else{
            p <- ggplot(aes(x, y, x.type = x.type, y.type = y.type, x.missing = x.missing, y.missing = y.missing, rf=rf, fill = LOD, LOD.CC=LOD.CC, LOD.CR=LOD.CR, LOD.RC=LOD.RC, LOD.RR=LOD.RR), data=df.graph) +
                geom_tile() +
                scale_fill_gradientn(colours = rev(grDevices::rainbow(n.colors)), na.value = "white") +
                theme(axis.text.x=element_text(angle=90, hjust=1))
        }

    ## If inbred:
    }else{
        if(graph.LOD!=TRUE){
            p <- ggplot(aes(x, y, x.missing = x.missing, y.missing = y.missing, fill=rf, LOD=LOD), data=df.graph) +
                geom_tile() +
                scale_fill_gradientn(colours = grDevices::rainbow(n.colors), na.value = "white") +
                theme(axis.text.x=element_text(angle=90, hjust=1))
        }else{
            p <- ggplot(aes(x, y, x.missing = x.missing, y.missing = y.missing, rf=rf, fill=LOD), data=df.graph) +
                geom_tile() +
                scale_fill_gradientn(colours = rev(grDevices::rainbow(n.colors)), na.value = "white") +
                theme(axis.text.x=element_text(angle=90, hjust=1))
                }
    }

    ## Disable lab names:
    if(is.null(lab.xy)){
        p <- p + labs(x = " ", y = " ")
    } else {
      if(length(lab.xy)!=2){
        warning("You should give a character vector with two components to axis labels")
        }else{
          p <- p + labs(x = lab.xy[1], y = lab.xy[2])
        }
    }

    ## Disable markers names:
    if(mrk.axis=="none"){
        p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
     }

    ## Write main
    if(!is.null(main)){
        p <- p + ggtitle(main)
    }

    ## Interactive    
    if(inter){
      if(is.null(html.file)){
        stop("For iteractive mode you must define a name for the outputted html file in 'html.file' argument.")
      }else{
        p <- plotly::ggplotly(p)
        htmlwidgets::saveWidget(p, file = html.file)
        utils::browseURL(html.file)
      }
    }else{
        p #it is a ggplot which can be expanded (+).
    }
}
