#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: segregations.tests.R                                          #
# Contains: chisq.test.for.segregation.of.markers,                    #
# test.segregation, plot.onemap.segreg.test,                          #
# print.onemap.segreg.test                                            #
#                                                                     #
# Written by Antonio Augusto Franco Garcia                            #
# copyright (c) 2015 Antonio Augusto Franco Garcia                    #
#                                                                     #
# First version: 2015/04/18                                           #
# Last update: 2015/04/18                                             #
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
    p.a <- rep(1/4, 4); p.b <- c(1/4, 1/2, 1/4); p.c <- c(3/4,1/4); p.d <- rep(1/2, 2)
    # Counting each category
    count <- table(x$geno[,marker], exclude=0)
    # Do the chisq test, using the appropriate expected segregation
    # grepl() allows finding the marker type (it has the letter in the argument)
    if (grepl("A",x$segr.type[marker])) {
        qui <- chisq.test(count, p=p.a, correct = FALSE)
        H0 <- "1:1:1:1" }
    else if (grepl("B",x$segr.type[marker])) {
        qui <- chisq.test(count, p=p.b, correct = FALSE)
        H0 <- "1:2:1" }
    else if (grepl("C",x$segr.type[marker])) {
        qui <- chisq.test(count, p=p.c, correct = FALSE)
        H0 <- "3:1" }
    else if (grepl("D",x$segr.type[marker])) {
        qui <- chisq.test(count, p=p.d, correct = FALSE)
        H0 <- "1:1" }
    return(list(Hypothesis=H0, qui.quad=qui$statistic, p.val=qui$p.value,
                perc.genot=100*(sum(table(x$geno[,marker], exclude=0))/x$n.ind)))
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
    if (is(x,"bc.onemap")|is(x,"f2.onemap")|is(x,"riself.onemap")|
        is(x,"risib.onemap")|is(x,"outcross")) {
        y <- list(Marker=dimnames(x$geno)[[2]],
                     Results.of.tests=sapply(1:x$n.mar, function(onemap.object, marker)
                         chisq.test.for.segregation.of.markers(onemap.object, marker),
                         onemap.object=x))
        # sapply iterates from 1 to x$n.mar; x is fixed (onemap object with data)
        class(y) <- c("onemap.segreg.test")
        invisible(y) #returns y without showing it
    }
    else stop("This is not a onemap object with raw data")
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
    Z <- data.frame(Marker=x$Marker,
                    H0=unlist(x$Results.of.tests[1,]),
                    Chi.square=unlist(x$Results.of.tests[2,]),
                    p.value=unlist(x$Results.of.tests[3,]),
                    Perc.genot=round(unlist(x$Results.of.tests[4,]),2))
    colnames(Z) <- c("Marker","H0","Chi-square","p-value","% genot.")
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
##' ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)
##'
##' data(example.out) # load OneMap's fake dataset for an outcrossing population
##' Out.seg <- test.segregation(example.out) # Applies chi-square tests
##' print(Out.seg) # Shows the results
##' plot(Out.seg) # Plot the graph, ordering the p-values
##' plot(Out.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(Out.seg)
##' ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)
##' 
##' @export
plot.onemap.segreg.test <- function(x, order=TRUE) {
    # Create a data frame
    Z <- data.frame(Marker=x$Marker,
                    X.square=unlist(x$Results.of.tests[2,]),
                    p.value=unlist(x$Results.of.tests[3,]))
    Bonf <- -log10(.05/nrow(Z)) #Bonferroni's threshold'
    Z$signif <- factor(ifelse(-log10(Z$p.value)<Bonf,"non sign.","sign."))
    Z$order <- 1:nrow(Z)
    # % of distorted
    perc <- 100*(1-(table(Z$signif)[1]/nrow(Z)))
    # Keeping markers in their original order (not alphanumeric), or by p-values (default)
    if (order!=TRUE) Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
    else Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$p.value, decreasing=TRUE)])
    # Plotting
    g <- ggplot(data=Z, aes(x=Marker, y=-log10(p.value)))
    g <- g + ylab(expression(-log[10](p-value)))
    g <- g + geom_point(aes(color=signif), stat="identity",size=2.5)
    g <- g + scale_colour_manual(name=paste("Bonferroni\n","(",round(perc,0),"% distorted)",sep=""),
                                 values = c("#46ACC8","#B40F20"))
    g <- g + geom_hline(yintercept = Bonf, colour="#E58601", linetype = "longdash")
    g <- g + coord_flip()
    if (nrow(Z)>30) g <- g + theme(axis.text.y = element_blank())
    g
}
##'

data(fake.bc.onemap) # load OneMap's fake dataset for a backcross population
BC.seg <- test.segregation(fake.bc.onemap) # Applies chi-square tests
print(BC.seg) # Shows the results
plot(BC.seg) # Plot the graph, ordering the p-values
plot(BC.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset

g <- plot(BC.seg)
ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)


data(example.out) # load OneMap's fake dataset for an outcrossing population
Out.seg <- test.segregation(example.out) # Applies chi-square tests
print(Out.seg) # Shows the results
plot(Out.seg) # Plot the graph, ordering the p-values
plot(Out.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
g <- plot(Out.seg)
ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)




chisq.test.for.segregation.of.markers(BC,1)
chisq.test.for.segregation.of.markers(F2,1)
chisq.test.for.segregation.of.markers(example.out,1)
chisq.test.for.segregation.of.markers(fake.bc.onemap,1)
chisq.test.for.segregation.of.markers(fake.f2.onemap,1)
chisq.test.for.segregation.of.markers(RIL1,1)

test.segregation(BC)
test.segregation(F2)
test.segregation(example.out)
test.segregation(fake.bc.onemap)
test.segregation(fake.f2.onemap)
test.segregation(RIL1)

BC. <- test.segregation(BC)
F2. <- test.segregation(F2)
OUT. <- test.segregation(example.out)
FAKE.BC. <- test.segregation(fake.bc.onemap)
FAKE.F2. <- test.segregation(fake.f2.onemap)
RIL. <- test.segregation(RIL1)

plot(BC)

plot(BC.)
plot(F2.)
plot(OUT.)
plot(FAKE.BC.)
plot(FAKE.F2.)
plot(RIL.)

print(BC.)
print(F2.)
print(OUT.)
print(FAKE.BC.)
print(FAKE.F2.)
print(RIL.)


BC.



names(BC)
BC
x <- BC
table(x$geno[,1], exclude=0)
sum(table(x$geno[,1], exclude=0))
BC$n.ind
BC$n.ind/sum(table(x$geno[,1], exclude=0))



pvalue.chisq.markers(BC,1)

?chisq.test
x <- matrix(c(12, 5, 7, 7), ncol = 2)
chisq.test(x)$p.value           # 0.4233
chisq.test(x, simulate.p.value = TRUE, B = 10000)$p.value
chisq.test(x)$statistic



test.segreg(BC)
Zz <- test.segreg(BC)
Zz
names(Zz)
class(Zz)
Zz$Marker
Zz$Values.of.X.and.p[1,]
Zz$Values.of.X.and.p[2,]
Zz$Values.of.X.and.p[[1]][1]


AA <- data.frame(Marker=test.segreg(BC)$Marker,
                X.square=test.segreg(BC)$Values.of.X.and.p[[1]][1],
                 p.value=test.segreg(BC)$Values.of.X.and.p[[2]][1])
AA




foo <- test.segreg(BC)
plot.onemap.segreg.test(foo)
plot(foo)

foo <- test.segreg(F2)
plot(foo)

foo <- test.segreg(example.out)
plot(foo)

foo <- test.segreg(fake.f2.onemap)
plot(foo)

foo <- test.segreg(RIL1)
plot(foo)



plot.chisquare(F2)
plot.chisquare(example.out)
plot.chisquare(fake.f2.onemap)
plot.chisquare(RIL1)





x <- BC
Z <- data.frame(Marker=test.segreg(x)$Marker,
                X.square=unlist(test.segreg(x)$Values.of.X.and.p[1,]),
                p.value=unlist(test.segreg(x)$Values.of.X.and.p[2,]))
Z
Bonf <- -log10(.05/nrow(Z)) #Bonferroni's threshold'
Bonf
Z$signif <- factor(ifelse(-log10(Z$p.value)<Bonf,"non sign.","sign."))
Z$order <- 1:nrow(Z)
Z
perc <- 100*(1-(table(Z$signif)[1]/nrow(Z)))
perc
# Keeping markers in their original order (not alphanumeric), or by p-values (default)
Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
Z
g <- ggplot(data=Z, aes(x=Marker, y=-log10(p.value)))
g <- g + ylab(expression(-log[10](p-value)))
g <- g + geom_point(aes(color=signif), stat="identity",size=2.5)
g <- g + scale_colour_manual(name=paste("Bonferroni\n","(",round(perc,0),"% distorted)",sep=""),
                             values = c("#46ACC8","#B40F20"))
g <- g + geom_hline(yintercept = Bonf, colour="#E58601", linetype = "longdash")
g <- g + coord_flip()
if (nrow(Z)>30) g <- g + theme(axis.text.y = element_blank())
g




x <- test.segreg(BC)
a <- x$Marker
a
b <- unlist(x$Values.of.X.and.p[1,])
b
c <- unlist(x$Values.of.X.and.p[2,])
c
AAA <- data.frame(M=a,X=b,p=c)
AAA

class(b)
b$1





BCseg <- test.segreg(BC)
plot.onemap.segreg.test(BCseg)
plot(AA)

plot.chisquare(F2)
plot.chisquare(example.out)
plot.chisquare(fake.f2.onemap)
plot.chisquare(RIL1)



plot.chisquare(BC)
plot.chisquare(F2)
plot.chisquare(example.out)
plot.chisquare(fake.f2.onemap)
plot.chisquare(RIL1)

plot.chisquare(BC, order=FALSE)
plot.chisquare(F2, order=FALSE)
plot.chisquare(example.out, order=FALSE)
plot.chisquare(fake.f2.onemap, order=FALSE)
plot.chisquare(RIL1, order=FALSE)




library(devtools)
install_github("augusto-garcia/onemap")
library(onemap)
remove.packages("onemap")
install_github("karthik/wesanderson")
library(wesanderson)
wes_palette("FantasticFox")[4]



fake.bc.onemap
names(fake.bc.onemap)
data(fake.bc.onemap)
fake.bc.onemap$geno
ht(fake.bc.onemap$geno)
fake.bc.onemap$geno[,1]
table(fake.bc.onemap$geno[,1])
?table
table(fake.bc.onemap$geno[,1], exclude=0)

x <- table(fake.bc.onemap$geno[,1], exclude=0)
x
length(x)
p=c(.5,.5)
length(p)

prop.test(x, p=c(.5))

BC <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse.raw")
x <- table(BC$geno[,1], exclude=0)
prop.test(x, p=c(.5), correct = FALSE) # Melhor não aplicar a correção de Yates!
chisq.test(x, p=c(.5,.5), correct = FALSE) # Melhor não aplicar a correção de Yates!
names(chisq.test(x, p=c(.5,.5), correct = FALSE))
chisq.test(x, p=c(.5,.5), correct = FALSE)$p.value
names(BC)
BC$segr.type
grepl("C*", BC$segr.type[1])
?grepl
grepl("C", "D1.10")


?chisq.test
BC$segr.type[1]

data(example.out)
names(example.out)



## Versão antiga, sem o sapply
test.segreg1 <- function(x) {
    df <- data.frame(Marker=rep(0,x$n.mar), p.value=rep(0,x$n.mar))
    for (i in 1:x$n.mar) {
        df$Marker[i] <- dimnames(x$geno)[[2]][i]
        df$p.value[i] <- pvalue.chisq.markers(x,i)
    }
    return(df)
}
#




dimnames(BC$geno)[[2]]


BC
m <- matrix(data=cbind(rnorm(30, 0), rnorm(30, 2), rnorm(30, 5)), nrow=30, ncol=3)
sapply(1:3, function(x, y) mean(y[,x]), y=m)

sapply(1:14, function(x,marker) pvalue.chisq.markers(x,marker), x=BC)


test.segreg(BC)
test.segreg(example.out)
test.segreg(fake.f2.onemap)
test.segreg(RIL1)

system.time(test.segreg1(BC))
system.time(test.segreg(BC))
system.time(test.segreg1(example.out))
system.time(test.segreg(example.out))
system.time(test.segreg1(fake.f2.onemap))
system.time(test.segreg(fake.f2.onemap))
system.time(test.segreg1(RIL1))
system.time(test.segreg(RIL1))


# Graphic
data <- Z
data$p.value <- -log10(as.numeric(data$p.value))
data

Z <- data.frame(test.segreg(BC))
Z
is.data.frame(Z)
Z$signif <- factor(ifelse(Z$p.value<.05,1,2))
Z
Z$order <- 1:nrow(Z)
Z    
Z$Marker
Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
Z$Marker  # notice the changed order of factor levels

g <- ggplot(data=Z, aes(x=Marker, y=p.value))
g <- g + geom_point(aes(color=signif), stat="identity",size=2.5)
g <- g + scale_colour_manual(name="With Bonferroni",
                             values = c("red", "darkblue"),
                             labels = c("sign.","non sign."))
g <- g + geom_hline(yintercept = 0.05, colour="darkred", linetype = "longdash")
#g <- g + coord_flip()
g









library(doParallel)
detectCores()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
getDoParWorkers()

library(foreach)
install.packages("doSNOW")

library(foreach) 
library(doSNOW) 
getDoParWorkers() 
getDoParName() 
registerDoSNOW(makeCluster(8, type = "SOCK")) 
getDoParWorkers() 
getDoParName() 

##
test.segreg2 <- function(x) {
    p.value <- foreach(i = 1:x$n.mar,.combine=c,.export="pvalue.chisq.markers") %dopar% {
        pvalue.chisq.markers(x,i)
    }
    return(data.frame(cbind(Marker=dimnames(x$geno)[[2]],p.value)))
}
#

x <- BC
x$Marker <- dimnames(x$geno)[[2]]
x
dimnames(x$geno)[[2]]

print(test.segreg2(BC))
test.segreg2(example.out)
test.segreg2(fake.f2.onemap)
test.segreg2(RIL1)
test.segreg(RIL1)



system.time(test.segreg(RIL1))
system.time(test.segreg2(RIL1))
system.time(test.segreg(fake.f2.onemap))
system.time(test.segreg2(fake.f2.onemap))
system.time( replicate(500, test.segreg(example.out)))
system.time( replicate(500, test.segreg2(example.out)))





x <- example.out
p <- foreach(i = 1:x$n.mar, .combine = c) %dopar% pvalue.chisq.markers(x,i)
p


library(onemap)



test.segreg(BC)
test.segreg(example.out)
test.segreg(fake.f2.onemap)
is.matrix(BC$geno)
class(BC$geno)
dimnames(BC$geno)
dimnames(BC$geno)[[2]][1]

RIL1 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="feijao.raw")
plot(RIL1)
test.segreg(RIL1)
z <- test.segreg(RIL1)
z

data(fake.f2.onemap)
plot(BC)

BC$geno[,1]
BC$geno.mmk
BC[[1]]
BC$1
names(BC)


m0 <- matrix(NA, 4, 0)
rownames(m0)



chisq.markers(BC, "M1")
chisq.markers(BC, 1)
chisq.markers(BC, 5)

F2 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="maize.raw")

chisq.markers(F2, 1)



data(example.out)
print(example.out)

chisq.markers(example.out, 1)
chisq.markers(example.out, 2)
chisq.markers(example.out, 3)
chisq.markers(example.out, 4)
chisq.markers(example.out, 5)
chisq.markers(example.out, 13)


table(example.out$segr.type[13])
example.out$segr.type[13]
table(example.out$geno[,13], exclude=0)


