library(devtools)
install_github("augusto-garcia/onemap")
library(onemap)
install.packages("reshape")
library(reshape)

###

# Read data using OneMap
BC <- read.mapmaker(dir="~/Dropbox/Disciplinas/LGN\ 5830\ -\ Biometria\ de\ Marcadores\ Genéticos/Dados\ experimentais", file="mouse.raw")

# See atributes of object BC
class(BC)
str(BC)
names(BC)
BC$geno
BC$n.ind

# Convert genotype data to a data frame with long format
data.BC <- data.frame(BC$geno)
data.BC <- melt(data.BC)
data.BC <- cbind(ind=1:BC$n.ind, data.BC)

# Build the graphic
g <- ggplot(data=data.BC, aes(x=ind, y=variable, fill=factor(value)))
g <- g + geom_tile()
g <- g + xlab("Individual") + ylab("Marker") + scale_fill_discrete(name="Genotype",
                                                                   labels=c("A","H","."))
g <- g + theme(axis.text.y = element_blank())
g

