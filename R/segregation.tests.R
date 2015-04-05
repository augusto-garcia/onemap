fake.bc.onemap
names(fake.bc.onemap)
library(onemap)
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

# x: onemap object with data, marker: marker to test
pvalue.chisq.markers <- function(x, marker) {
    # Segregation pattern for each marker type
    p.a <- rep(1/4, 4); p.b <- c(1/4, 1/2, 1/4); p.c <- c(3/4,1/4); p.d <- rep(1/2, 2)
    # Counting each category
    count <- table(x$geno[,marker], exclude=0)
    # Do the chisq test, using the appropriate expected segregation
    # grepl() allows finding the marker type (it has the letter in the argument)
    if (grepl("A",x$segr.type[marker]))
        p.val <- chisq.test(count, p=p.a, correct = FALSE)$p.value
    else if (grepl("B",x$segr.type[marker]))
        p.val <- chisq.test(count, p=p.b, correct = FALSE)$p.value
    else if (grepl("C",x$segr.type[marker]))
        p.val <- chisq.test(count, p=p.c, correct = FALSE)$p.value
    else if (grepl("D",x$segr.type[marker]))
        p.val <- chisq.test(count, p=p.d, correct = FALSE)$p.value
    return(p.val)
}

##
test.segreg <- function(x) {
    df <- data.frame(Marker=rep(0,x$n.mar), p.value=rep(0,x$n.mar))
    for (i in 1:x$n.mar) {
        df$Marker[i] <- dimnames(x$geno)[[2]][i]
        df$p.value[i] <- pvalue.chisq.markers(x,i)
    }
    return(df)
}
#


test.segreg(BC)
test.segreg(example.out)
test.segreg(fake.f2.onemap)
test.segreg(RIL1)


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


