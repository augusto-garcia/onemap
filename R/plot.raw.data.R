#####
plot.bc.onemap <- function(df) {
    # Creating the data frame
    df.BC <- data.frame(df$geno)
    df.BC <- melt(df.BC) # function from package reshape
    df.BC <- cbind(ind=1:df$n.ind, df.BC)
    df.BC$value <- factor(df.BC$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.BC$value)==c("0","1","2")))) labels.bc <- c("-","AA","AB")
    else if (all(levels(df.BC$value)==c("1","2"))) labels.bc <- c("AA","AB")
    # Plotting
    g <- ggplot(data=df.BC, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype",labels=labels.bc, values=c("#F21A00","#3B9AB2","#EBCC2A")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###

#####
plot.riself.onemap <- function(df) {
    # Creating the data frame
    df.RIL <- data.frame(df$geno)
    df.RIL <- melt(df.RIL) # function from package reshape
    df.RIL <- cbind(ind=1:df$n.ind, df.RIL)
    df.RIL$value <- factor(df.RIL$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.RIL$value)==c("0","1","2")))) labels.ril <- c("-","AA","BB")
    else if (all(levels(df.RIL$value)==c("1","2"))) labels.ril <- c("AA","BB")
    # Plotting
    g <- ggplot(data=df.RIL, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",labels=labels.ril, values=c("#F21A00","#3B9AB2","#EBCC2A")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###

#####
plot.risib.onemap <- function(df) {
    # Creating the data frame
    df.RIL <- data.frame(df$geno)
    df.RIL <- melt(df.RIL) # function from package reshape
    df.RIL <- cbind(ind=1:df$n.ind, df.RIL)
    df.RIL$value <- factor(df.RIL$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.RIL$value)==c("0","1","2")))) labels.ril <- c("-","AA","BB")
    else if (all(levels(df.RIL$value)==c("1","2"))) labels.ril <- c("AA","BB")
    # Plotting
    g <- ggplot(data=df.RIL, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",labels=labels.ril, values=c("#F21A00","#3B9AB2","#EBCC2A")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###

#####
plot.f2.onemap <- function(df) {
    # Creating the data frame
    df.F2 <- data.frame(df$geno.mmk[[1]])
    df.F2 <- melt(df.F2) # function from package reshape
    df.F2[is.na(df.F2)] <- 0 # To avoid problems with NAs
    df.F2 <- cbind(ind=1:df$n.ind, df.F2)
    df.F2$value <- factor(df.F2$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4","5")))) labels.f2 <-
        c("-","AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4","5")))) labels.f2 <- c("AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3")))) labels.f2 <- c("-","AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3")))) labels.f2 <- c("AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","5")))) labels.f2 <- c("-","AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","5")))) labels.f2 <- c("AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4")))) labels.f2 <- c("-","AA","AB","BB","not BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4")))) labels.f2 <- c("AA","AB","BB","not BB")
    # Plotting
    g <- ggplot(data=df.F2, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype", labels=labels.f2, values=c("#000000", "#00A08A", "#5BBCD6", "#F2AD00", "#F98400", "#FF0000")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###


###
# FUNCIONANDO!
# input: x (an object of class outcross)
# output: a dataframe suitable for ggplot2
create.dataframe.outcross.plot <- function(x) {
    for (i in 1:4) {
    if (length(which(x$segr.type==paste("A.",i,sep=""))!=0)) {
        F1.A <- data.frame(x$geno[,which(x$segr.type==paste("A.",i,sep=""))])
        colnames(F1.A) <- colnames(x$geno)[which(x$segr.type==paste("A.",i,sep=""))]
        F1.A <- melt(F1.A)
        F1.A <- cbind(ind=1:x$n.ind, F1.A, Mk.type=paste("A.",i,sep=""))
        F1.A$value <- factor(F1.A$value)
        if (i==1 ) F1.A. <- data.frame(F1.A)
        else F1.A. <- rbind(F1.A.,F1.A)
        }
    }
    class(F1.A.) <- "data.frame"
    return(F1.A.)
}
###

ZZ <- create.dataframe.outcross.plot(example.out)
class(ZZ)

##
plot.outcross <- function(df, all=TRUE) {
    df <- create.dataframe.outcross.plot(df)
    g <- ggplot(data=df, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype",
                          values=c("#00A08A", "#5BBCD6",  "#F2AD00", "#F98400", "#FF0000")) 
    if (length(levels(df$variable))>20) g <- g +
                          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (all==TRUE) g
    else g + facet_grid( . ~ Mk.type)
    #else g + facet_wrap( ~ Mk.type, ncol=3)

}
##

plot(example.out)
plot(example.out,all=FALSE)


###
# NOVA TENTATIVA, INCLUINDO B, C, D1 e D2
# SE FUNCIONAR, É A FUNÇÃO DESEJADA!!!!
# input: x (an object of class outcross)
# output: a dataframe suitable for ggplot2
create.dataframe.outcross.plot <- function(x) {
    # Markers of A type
    for (i in 1:4) {
    if (length(which(x$segr.type==paste("A.",i,sep=""))!=0)) {
        F1.A <- data.frame(x$geno[,which(x$segr.type==paste("A.",i,sep=""))])
        colnames(F1.A) <- colnames(x$geno)[which(x$segr.type==paste("A.",i,sep=""))]
        F1.A <- melt(F1.A)
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
        F1.B <- melt(F1.B)
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
        F1.C <- melt(F1.C)
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
        F1.D1 <- melt(F1.D1)
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
        F1.D2 <- melt(F1.D2)
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
###

create.dataframe.outcross.plot(example.out)


plot(example.out)
plot(example.out,all=FALSE)


example.out$segr.type
paste("A.",1,2,sep="")
i <- 1
paste("B",i,".",i+4,sep="")
rm(i)

###
## atenção: usada apenas para testes; apagar depois
# input: x (an object of class outcross)
# output: a dataframe suitable for ggplot2
create.dataframe.outcross.plot <- function(x) {
    A.1 <- 0; A.2 <- 0; A.3 <- 0; A.4 <- 0
    if (length(which(x$segr.type=="A.1"))!=0) {
        F1.A1 <- data.frame(x$geno[,which(x$segr.type=="A.1")])
        colnames(F1.A1) <- colnames(x$geno)[which(x$segr.type=="A.1")]
        F1.A1 <- melt(F1.A1)
        F1.A1 <- cbind(ind=1:x$n.ind, F1.A1, Mk.type="A.1")
        F1.A1$value <- factor(F1.A1$value)
        A.1 <- 1
    }
    if (length(which(x$segr.type=="A.2"))!=0) {
        F1.A2 <- data.frame(x$geno[,which(x$segr.type=="A.2")])
        colnames(F1.A2) <- colnames(x$geno)[which(x$segr.type=="A.2")]
        F1.A2 <- melt(F1.A2)
        F1.A2 <- cbind(ind=1:x$n.ind, F1.A2, Mk.type="A.2")
        F1.A2$value <- factor(F1.A2$value)
        A.2 <- 1
    }
    if (length(which(x$segr.type=="A.3"))!=0) {
        F1.A3 <- data.frame(x$geno[,which(x$segr.type=="A.3")])
        colnames(F1.A3) <- colnames(x$geno)[which(x$segr.type=="A.3")]
        F1.A3 <- melt(F1.A3)
        F1.A3 <- cbind(ind=1:x$n.ind, F1.A3, Mk.type="A.3")
        F1.A3$value <- factor(F1.A3$value)
        A.3 <- 1
    }
    if (length(which(x$segr.type=="A.4"))!=0) {
        F1.A4 <- data.frame(x$geno[,which(x$segr.type=="A.4")])
        colnames(F1.A4) <- colnames(x$geno)[which(x$segr.type=="A.4")]
        F1.A4 <- melt(F1.A4)
        F1.A4 <- cbind(ind=1:x$n.ind, F1.A4, Mk.type="A.4")
        F1.A4$value <- factor(F1.A4$value)
        A.4 <- 1
    }
    for (i in 1:4){
    if (paste("A.",i,sep="")==1) F1.A <- data.frame(rbind(paste("F1.A",i,sep="")))
    }
F1.A
}
###



ZZ <- create.dataframe.outcross.plot(example.out)
ZZ

g <- ggplot(data=ZZ, aes(x=ind, y=variable, fill=factor(value)))
g <- g + geom_tile()
g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",
                      values=c("#00A08A", "#5BBCD6",  "#F2AD00", "#F98400", "#FF0000")) 
if (length(levels(ZZ$variable))>20) g <- g +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

g # All together

g <- g + facet_grid(.~Mk.type)
g # Panels


for (i in 1:4){
    if (paste("A.",i,sep="")==1) F1.A <- data.frame(rbind(paste("F1.A",i,sep="")))
}
F1.A


paste("A.",1,sep="")

?paste
F1.A <- data.frame(rbind(F1.A1,F1.A2,F1.A3,F1.A4))



data(example.out)
print(example.out)
names(example.out)
# Creating the data frame
example.out$geno
colnames(example.out$geno)[12]
example.out$segr.type.num
example.out$segr.type

# Incluindo dados perdidos
example.out$geno[1,4] <- 0
example.out$geno[2,4] <- 0
example.out$geno[3,4] <- 0






# Separating A.1 type
which(example.out$segr.type=="A.1")
example.out$geno[,which(example.out$segr.type=="A.1")]
F1.A1 <- data.frame(example.out$geno[,which(example.out$segr.type=="A.1")])
colnames(F1.A1) <- colnames(example.out$geno)[which(example.out$segr.type=="A.1")]

F1.A1
F1.A1 <- melt(F1.A1)
ht(F1.A1)
F1.A1 <- cbind(ind=1:example.out$n.ind, F1.A1, Mk.type="A.1")
F1.A1$value <- factor(F1.A1$value)
ht(F1.A1)

# Separating A.2 type
which(example.out$segr.type=="A.2")
example.out$geno[,which(example.out$segr.type=="A.2")]
# Usar o código abaixo quando houver um único marcador selecionado
F1.A2 <- data.frame(example.out$geno[,which(example.out$segr.type=="A.2")])
colnames(F1.A2) <- colnames(example.out$geno)[which(example.out$segr.type=="A.2")]

F1.A2
F1.A2 <- melt(F1.A2)
ht(F1.A2)
F1.A2 <- cbind(ind=1:example.out$n.ind, F1.A2, Mk.type="A.2")
F1.A2$value <- factor(F1.A2$value)
ht(F1.A2)

# Separating A.3 type
length(which(example.out$segr.type=="A.3"))

example.out$geno[,which(example.out$segr.type=="A.3")]
F1.A3 <- data.frame(example.out$geno[,which(example.out$segr.type=="A.3")])
F1.A3
F1.A3 <- melt(F1.A3)
ht(F1.A3)
F1.A3 <- cbind(ind=1:example.out$n.ind, F1.A3, Mk.type="A.3")
F1.A3$value <- factor(F1.A3$value)
ht(F1.A3)

# Separating A.4 type
which(example.out$segr.type=="A.4")
example.out$geno[,which(example.out$segr.type=="A.4")]
F1.A4 <- data.frame(example.out$geno[,which(example.out$segr.type=="A.4")])
F1.A4
F1.A4 <- melt(F1.A4)
ht(F1.A4)
F1.A4 <- cbind(ind=1:example.out$n.ind, F1.A4, Mk.type="A.4")
F1.A4$value <- factor(F1.A4$value)
ht(F1.A4)


## All

F1.A <- data.frame(rbind(F1.A1,F1.A2,F1.A4))
F1.A
head(F1.A)

F1.A$n.mar
levels(F1.A$variable)
length(levels(F1.A$variable))

g <- ggplot(data=F1.A, aes(x=ind, y=variable, fill=factor(value)))
g <- g + geom_tile()
g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",
                      values=c("#00A08A", "#5BBCD6",  "#F2AD00", "#F98400", "#FF0000")) 
if (length(levels(F1.A$variable))>20) g <- g +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

g # All together

g <- g + facet_grid(.~Mk.type)
g # Panels



# Plotting
g <- ggplot(data=F1.A1, aes(x=ind, y=variable, fill=factor(value)))
g <- g + geom_tile()
g







library(devtools)
install_github("augusto-garcia/onemap")
library(onemap)
remove.packages("reshape")
install.packages("reshape2")
library(reshape2)


###

# Read data using OneMap
BC <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse.raw")

print(BC)

# Read data using OneMap
BC1 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse.raw")
BC2 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse2.raw")
BC3 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse3.raw")
BC4 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse4.raw")

#####
print.onemap.raw <- function(df) {
    # Creating the data frame
    df.BC <- data.frame(df$geno)
    df.BC <- melt(df.BC) # function from package reshape
    df.BC <- cbind(ind=1:df$n.ind, df.BC)
    df.BC$value <- factor(df.BC$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.BC$value)==c("0","1","2")))) labels.bc <- c("-","AA","AB")
    else if (all(levels(df.BC$value)==c("1","2"))) labels.bc <- c("AA","AB")
    # Plotting
    g <- ggplot(data=df.BC, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype",labels=labels.bc, values=c("#F21A00","#3B9AB2","#EBCC2A")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###

RIL1 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="feijao.raw")
RIL1

#####
print.onemap.raw2 <- function(df) {
    # Creating the data frame
    df.RIL <- data.frame(df$geno)
    df.RIL <- melt(df.RIL) # function from package reshape
    df.RIL <- cbind(ind=1:df$n.ind, df.RIL)
    df.RIL$value <- factor(df.RIL$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.RIL$value)==c("0","1","2")))) labels.ril <- c("-","AA","BB")
    else if (all(levels(df.RIL$value)==c("1","2"))) labels.ril <- c("AA","BB")
    # Plotting
    g <- ggplot(data=df.RIL, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",labels=labels.ril, values=c("#F21A00","#3B9AB2","#EBCC2A")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###

#####
print.onemap.raw3 <- function(df) {
    # Creating the data frame
    df.F2 <- data.frame(df$geno.mmk[[1]])
    df.F2 <- melt(df.F2) # function from package reshape
    df.F2[is.na(df.F2)] <- 0 # To avoid problems with NAs
    df.F2 <- cbind(ind=1:df$n.ind, df.F2)
    df.F2$value <- factor(df.F2$value)
    # Defining the label for genotypes
    if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4","5")))) labels.f2 <-
        c("-","AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4","5")))) labels.f2 <- c("AA","AB","BB","not BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3")))) labels.f2 <- c("-","AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3")))) labels.f2 <- c("AA","AB","BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","5")))) labels.f2 <- c("-","AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","5")))) labels.f2 <- c("AA","AB","BB","not AA")
    else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4")))) labels.f2 <- c("-","AA","AB","BB","not BB")
    else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4")))) labels.f2 <- c("AA","AB","BB","not BB")
    # Plotting
    g <- ggplot(data=df.F2, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
        scale_fill_manual(name="Genotype", labels=labels.f2, values=c("#000000", "#00A08A", "#5BBCD6", "#F2AD00", "#F98400", "#FF0000")) 
    if (df$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g
}
###


print.onemap.raw3(fake.f2.onemap)

F2 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="maize.raw")
print(F2)
print.onemap.raw3(F2)

F2. <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="maize2.raw")
print(F2.)
print.onemap.raw3(F2.)



RIL1 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="feijao.raw")
print.onemap.raw2(RIL1)

RIL2 <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="feijao2.raw")
print.onemap.raw2(RIL2)

data(fake.f2.onemap)
print(fake.f2.onemap)

names(fake.f2.onemap)
fake.f2.onemap$segr.type
fake.f2.onemap$geno.mmk


    # Creating the data frame
    df.F2 <- data.frame(fake.f2.onemap$geno.mmk[[1]])
    df.F2 <- melt(df.F2) # function from package reshape
    df.F2[is.na(df.F2)] <- 0 # To avoid problems with NAs
    df.F2 <- cbind(ind=1:fake.f2.onemap$n.ind, df.F2)
    df.F2$value <- factor(df.F2$value)

ht(df.F2)
df.F2[1:50,]


    # Defining the label for genotypes
if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4","5")))) labels.f2 <- c("-","AA","AB","BB","not BB","not AA")

else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4","5")))) labels.f2 <- c("AA","AB","BB","not BB","not AA")
else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3")))) labels.f2 <- c("-","AA","AB","BB")
else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3")))) labels.f2 <- c("AA","AB","BB")
else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","5")))) labels.f2 <- c("-","AA","AB","BB","not AA")
else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","5")))) labels.f2 <- c("AA","AB","BB","not AA")
else if (suppressWarnings(all(levels(df.F2$value)==c("0","1","2","3","4")))) labels.f2 <- c("-","AA","AB","BB","not BB")
else if (suppressWarnings(all(levels(df.F2$value)==c("1","2","3","4")))) labels.f2 <- c("AA","AB","BB","not BB")

labels.f2
levels(df.F2$value)


    # Plotting
    g <- ggplot(data=df.F2, aes(x=ind, y=variable, fill=factor(value)))
    g <- g + geom_tile()
    g <- g + xlab("Individual") + ylab("Marker") +
    scale_fill_manual(name="Genotype",labels=labels.f2, values=c("#000000","#FF0000","#00A08A","#F2AD00","#F98400","#5BBCD6")) 
    if (fake.f2.onemap$n.mar>20) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    g


wes_palette("Darjeeling")[5]
wes_palette("Darjeeling2")[5]


?all
?melt
# Colors
install.packages("wesanderson")
library(wesanderson)
names(wes_palette)
wes_palette("Zissou")

?wes_palette

df.BC3 <- data.frame(BC2$geno)
df.BC3
df.BC3 <- melt(df.BC) # function from package reshape
df.BC3 <- cbind(ind=1:df$n.ind, df.BC)
df.BC3$value <- factor(df.BC$value)

head(BC1$geno)
head(BC2$geno)
head(BC3$geno)
head(BC4$geno)


print.onemap.raw(BC1)
print.onemap.raw(BC2)
print.onemap.raw(BC3)
print.onemap.raw(BC4)


# See atributes of object BC
class(BC)
str(BC)
names(BC)
BC$geno
BC$n.ind
BC$n.mar
BC$segr.type

# Convert genotype data to a data frame with long format
data.BC <- data.frame(BC$geno)
data.BC <- melt(data.BC)
data.BC <- cbind(ind=1:BC$n.ind, data.BC)
data.BC$value <- factor(data.BC$value)
head(data.BC)

levels(data.BC$value)
labels.bc <- 0

    if (levels(data.BC$value)==c("1","2")) labels.bc <- c("Bla","Ble")

if (all(levels(data.BC$value)==c("1","2"))) labels.bc <- c("Bla","Ble")

labels.bc
labels=c("A","H","."))

colors <- c("blue", "green", "yellow")
# Build the graphic
g <- ggplot(data=data.BC, aes(x=ind, y=variable, fill=factor(value)))
g <- g + geom_tile()
g <- g + scale_fill_manual(name="Genotype",values=colors,labels=c("Ble","Bliu"))
g <- g + xlab("Individual") + ylab("Marker") + scale_fill_discrete(name="Genotype",
                                                                   labels=c("Ble","Bliu"))
g <- g + theme(axis.text.y = element_blank())
g

head(data.BC)


##########

#### Trying to use class S3

print(BC)
BC







plot(BC)
