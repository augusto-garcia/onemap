fake.bc.onemap
names(fake.bc.onemap)
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



