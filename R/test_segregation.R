#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: segregations.tests.R                                          #
# Contains: chisq.test.for.segregation.of.markers,                    #
# test.segregation, plot.onemap.segreg.test,                          #
# print.onemap.segreg.test, Bonferroni.alpha, subset.chisq            #
#                                                                     #
# Written by Antonio Augusto Franco Garcia                            #
# copyright (c) 2015 Antonio Augusto Franco Garcia                    #
#                                                                     #
# First version: 2015/04/18                                           #
# Last update: 2015/04/21                                             #
# License: GNU General Public License version 3 or later              #
#                                                                     #
#######################################################################

##' chisq.test.for.segregation.of.markers
##'
##' Applies the chi-square test to check if markers are following the
##' expected segregation pattern, i. e., 1:1:1:1 (A), 1:2:1 (B), 3:1 (C) and 1:1 (D)
##' according to OneMap's notation. It does not use Yate's correction.
##'
##' First, the function select the correct segregation pattern, then it
##' defines the H0 hypothesis, and then tests it, together with percentage of
##' missing data.
##'
##' @param x an object of class bc.onemap, f2.onemap, riself.onemap,
##' risib.onemap or outcross, with data and additional information.
##' @param marker the marker which will be tested for its segregation.
##' 
##' @return a list with the H0 hypothesis being tested, the chi-square statistics,
##' the associated p-values, and the % of individuals genotyped.
##'
##' @examples
##' data(fake.bc.onemap) # Loads a fake backcross dataset installed with onemap
##' chisq.test.for.segregation.of.markers(fake.bc.onemap,1)
##'
##' data(example.out) # Loads a fake outcross dataset installed with onemap
##' chisq.test.for.segregation.of.markers(example.out,1)
chisq.test.for.segregation.of.markers <- function(x, marker) {
    # Segregation pattern for each marker type
    p.a <- rep(1/4, 4)
    p.b <- c(1/4, 1/2, 1/4)
    p.c <- c(3/4, 1/4)
    p.d <- rep(1/2, 2)
    # Counting each category
    count <- table(x$geno[, marker], exclude = 0)
    # Do the chisq test, using the appropriate expected segregation grepl()
    # allows finding the marker type (it has the letter in the argument)
    if (grepl("A", x$segr.type[marker])) {
        qui <- chisq.test(count, p = p.a, correct = FALSE)
        H0 <- "1:1:1:1"
    } else if (grepl("B", x$segr.type[marker])) {
        qui <- chisq.test(count, p = p.b, correct = FALSE)
        H0 <- "1:2:1"
    } else if (grepl("C", x$segr.type[marker])) {
        qui <- chisq.test(count, p = p.c, correct = FALSE)
        H0 <- "3:1"
    } else if (grepl("D", x$segr.type[marker])) {
        qui <- chisq.test(count, p = p.d, correct = FALSE)
        H0 <- "1:1"
    }
    return(list(Hypothesis = H0, qui.quad = qui$statistic, p.val = qui$p.value, 
        perc.genot = 100 * (sum(table(x$geno[, marker], exclude = 0))/x$n.ind)))
}
##'

##' test.segregation
##'
##' Using the OneMap function chisq.test.for.segregation.of.markers,
##' performs the Chi-square test to check if all markers in a dataset are following
##' the expected segregation pattern, i. e., 1:1:1:1 (A), 1:2:1 (B), 3:1 (C) and 1:1 (D)
##' according to OneMap's notation.
##'
##' First, it identifies the correct segregation pattern and corresponding H0 hypothesis,
##' and then tests it.
##'
##' @param x an object of class bc.onemap, f2.onemap, riself.onemap,
##' risib.onemap or outcross, with data and additional information.
##' 
##' @return an object of class onemap.segreg.test, which is a list with marker name,
##' H0 hypothesis being tested, the chi-square statistics,  the associated p-values,
##' and the % of individuals genotyped. To see the object, it is necessary to print
##' it.
##' 
##' @examples
##' data(example.out) # Loads a fake outcross dataset installed with onemap
##' Chi <- test.segregation(example.out) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results
##'
##' @export
test.segregation <- function(x) {
    if (is(x, "bc.onemap") | is(x, "f2.onemap") | is(x, "riself.onemap") | 
        is(x, "risib.onemap") | is(x, "outcross")) {
        y <- list(Marker = dimnames(x$geno)[[2]], Results.of.tests = sapply(1:x$n.mar, 
            function(onemap.object, marker) chisq.test.for.segregation.of.markers(onemap.object, 
                marker), onemap.object = x))
        # sapply iterates from 1 to x$n.mar; x is fixed (onemap object with
        # data)
        class(y) <- c("onemap.segreg.test")
        invisible(y)  #returns y without showing it
    } else stop("This is not a onemap object with raw data")
}
##'

##' print.onemap.segreg.test
##'
##' It shows the results of Chisquare tests performed for all markers in a onemap object
##' of class class bc.onemap, f2.onemap, riself.onemap, risib.onemap or outcross.
##' 
##' @param x an object of class onemap.segreg.test
##' 
##' @return a dataframe with marker name, H0 hypothesis, chi-square statistics,
##' p-values, and % of individuals genotyped.
##' 
##' @examples
##' data(example.out) # Loads a fake outcross dataset installed with onemap
##' Chi <- test.segregation(example.out) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results
##'
##' @export
print.onemap.segreg.test <- function(x) {
    Z <- data.frame(Marker = x$Marker, H0 = unlist(x$Results.of.tests[1, 
        ]), Chi.square = unlist(x$Results.of.tests[2, ]), p.value = unlist(x$Results.of.tests[3, 
        ]), Perc.genot = round(unlist(x$Results.of.tests[4, ]), 2))
    colnames(Z) <- c("Marker", "H0", "Chi-square", "p-value", "% genot.")
    return(Z)
}
##'

##' plot.onemap.segreg.test
##' 
##' Draw a graphic showing the p-values (re-scaled to -log10(p-values)) associated with the
##' chi-square tests for the expected segregation patterns for all markers in a dataset.
##' It includes a vertical line showing the threshold for declaring statistical significance
##' if Bonferroni's correction is considered, as well as the percentage of markers that
##' will be discarded if this criteria is used.
##'
##' @param x an object of class onemap.segreg.test (produced by onemap's function
##' test.segregation()), i. e., after performing segregation tests
##' @param order a variable to define if p-values will be ordered in the plot
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##' 
##' @examples
##' data(fake.bc.onemap) # load OneMap's fake dataset for a backcross population
##' BC.seg <- test.segregation(fake.bc.onemap) # Applies chi-square tests
##' print(BC.seg) # Shows the results
##' plot(BC.seg) # Plot the graph, ordering the p-values
##' plot(BC.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(BC.seg)
##' ggsave('SegregationTests.jpg', g, width=7, height=5, dpi=600)
##'
##' data(example.out) # load OneMap's fake dataset for an outcrossing population
##' Out.seg <- test.segregation(example.out) # Applies chi-square tests
##' print(Out.seg) # Shows the results
##' plot(Out.seg) # Plot the graph, ordering the p-values
##' plot(Out.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(Out.seg)
##' ggsave('SegregationTests.jpg', g, width=7, height=5, dpi=600)
##' 
##' @export
plot.onemap.segreg.test <- function(x, order = TRUE) {
    # Create a data frame
    Z <- data.frame(Marker = x$Marker, X.square = unlist(x$Results.of.tests[2, 
        ]), p.value = unlist(x$Results.of.tests[3, ]))
    Bonf <- -log10(0.05/nrow(Z))  #Bonferroni's threshold'
    Z$signif <- factor(ifelse(-log10(Z$p.value) < Bonf, "non sign.", "sign."))
    Z$order <- 1:nrow(Z)
    # % of distorted
    perc <- 100 * (1 - (table(Z$signif)[1]/nrow(Z)))
    # Keeping markers in their original order (not alphanumeric), or by
    # p-values (default)
    if (order != TRUE) 
        Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)]) else Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$p.value, 
        decreasing = TRUE)])
    # Plotting
    g <- ggplot(data = Z, aes(x = Marker, y = -log10(p.value)))
    g <- g + ylab(expression(-log[10](p - value)))
    g <- g + geom_point(aes(color = signif), stat = "identity", size = 2.5)
    g <- g + scale_colour_manual(name = paste("Bonferroni\n", "(", round(perc, 
        0), "% distorted)", sep = ""), values = c("#46ACC8", "#B40F20"))
    g <- g + geom_hline(yintercept = Bonf, colour = "#E58601", linetype = "longdash")
    g <- g + coord_flip()
    if (nrow(Z) > 30) 
        g <- g + theme(axis.text.y = element_blank())
    g
}
##'


##' Bonferroni.alpha
##'
##' It shows the alpha value to be used in each test to control global type I error
##' for chi-square tests for segregation of all markers if Bonferroni's criteria is applied.
##' 
##' @param x an object of class onemap.segreg.test
##' @param global.alpha the global alpha that is target
##' 
##' @return the alpha value (numeric)
##' 
##' @examples
##' data(fake.bc.onemap) # Loads a fake backcross dataset installed with onemap
##' Chi <- test.segregation(fake.bc.onemap) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results of the Chi-square tests
##' Bonferroni.alpha (Chi) # Shows the global alpha level using Bonferroni's criteria
##'
##' @export
Bonferroni.alpha <- function(x, global.alpha = 0.05) {
    if (!is(x, "onemap.segreg.test")) 
        stop("This is not an object of class onemap.segreg.test")
    alpha.Bonf <- global.alpha/length(x$Marker)
    return(alpha.Bonf)
}
##'

##' subset.chisq
##'
##' A function to shows which marker have segregation distortion if Bonferroni's correction is
##' applied for the Chi-square tests of mendelian segregation.
##'
##' @param x an object of class onemap.segreg.test
##' @param distorted a TRUE/VALUE variable to show distorted or non-distorted markers
##' 
##' @return a vector with marker names, according to the option for 'distorted'
##' 
##' @examples
##' data(fake.bc.onemap) # Loads a fake backcross dataset installed with onemap
##' Chi <- test.segregation(fake.bc.onemap) # Performs the chi-square test for all markers
##' subset.chisq(Chi) # To show non-distorted markers
##' subset.chisq(Chi, distorted=TRUE) # With segregation distortion
##'
##' @export
subset.chisq <- function(x, distorted = FALSE) {
    if (!is(x, "onemap.segreg.test")) 
        stop("This is not an object of class onemap.segreg.test")
    Z <- data.frame(Marker = x$Marker, p.value = unlist(x$Results.of.tests[3, 
        ]))
    if (distorted == FALSE) 
        Z <- subset(Z, p.value >= Bonferroni.alpha(x)) else Z <- subset(Z, p.value < Bonferroni.alpha(x))
    return(as.vector(Z[, 1]))
}
##'
