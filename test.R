###########
# Packages
###########

library(vcfR)
library(polyRAD)
library(updog)

####################################
# Testing for outcrossing species
####################################

##########
# Using GQ
##########

vcf.out <- read.vcfR("vcf_example_out.vcf")
out <- onemap_read_vcfR(vcfR.object = vcf.out, cross = "outcross", parent1 = "P1", parent2 = "P2")

segr <- test_segregation(out)
plot(segr)

twopts <- rf_2pts(out)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

genotypes_errors <- extract_depth(vcfR.object = vcf.out, onemap.object = out, vcf.par = "GQ", parent1 = "P1", 
                                 parent2 = "P2", mean_phred = 20, recovering = FALSE)

new.errors <- create_probs(out, 
                           genotypes_errors = genotypes_errors, 
                           global_error = NULL, 
                           genotypes_probs = NULL)

twopts <- rf_2pts(new.errors)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,2)
map.lg1.1 <- map(lg1)

########
# updog
########

old.updog <- updog_error(vcfR.object = vcf.out, onemap.object = out, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

segr <- test_segregation(old.updog)
plot(segr)

twopts <- rf_2pts(old.updog)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1 <- make_seq(lgs,2)
map.lg1 <- map(lg1)

##########
# PolyRAD
##########

poly.test <- VCF2RADdata("vcf_example_out.vcf", phaseSNPs = FALSE, 
                         min.ind.with.reads = 0,
                         min.ind.with.minor.allele = 0)


poly.test <- SetDonorParent(poly.test, "P1")
poly.test <- SetRecurrentParent(poly.test, "P2")

mydata2 <- PipelineMapping2Parents(poly.test, 
                                   freqAllowedDeviation = 0.06,
                                   useLinkage = FALSE,
                                   minLikelihoodRatio = 2)

Export_MAPpoly(mydata2, "test")

genotypes <- read.table("test", skip=12)

pos <- sapply(strsplit(as.character(genotypes$V1), split = "_"),"[",1)

pos.onemap <- colnames(out$geno)
genotypes <- genotypes[which(pos%in%pos.onemap),]
keep.mks <- which(pos.onemap%in%pos)

# Atualizar geno
out$geno <- out$geno[,keep.mks]

new.geno <- vector()
for(i in 1:dim(genotypes)[1]){
  if(which.max(genotypes[i,3:5]) == 3){
    new.geno[i] <- 3
  }else if(which.max(genotypes[i,3:5]) == 2){
    new.geno[i] <- 2
  } else if(which.max(genotypes[i,3:5]) == 1){
    new.geno[i] <- 1
  }
}

new.geno <- matrix(new.geno,nrow = out$n.ind, ncol = length(keep.mks) )
colnames(new.geno) <- colnames(out$geno)
rownames(new.geno) <- rownames(out$geno)

# Mudando a ordem
genotypes <- genotypes[order(genotypes$V2),]

# Quanto mudou
sum(new.geno == out$geno)/length(new.geno)

out$geno <- new.geno
# Remover marcadores
out$n.mar <- length(keep.mks)
out$segr.type <- out$segr.type[keep.mks]
out$segr.type.num <- out$segr.type.num[keep.mks]
out$CHROM <- out$CHROM[keep.mks]
out$POS <- out$POS[keep.mks]

polyrad.one <- create_probs(onemap.obj = out, genotypes_probs =  genotypes[,3:5])
head(polyrad.one$error)

twopts <- rf_2pts(polyrad.one)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)


####################################
# Testing for f2 populations
####################################
# The example is without the F1
fix.f2 <- read.table("vcf_example_f2.vcf")
head(fix.f2)
new.f2 <- cbind(fix.f2[,1:11], rep("0/1:7,4:11:99:111,0,219",length(fix.f2[,1])), fix.f2[,12:dim(fix.f2)[2]])
write.table(new.f2, file = "new.f2", quote = F, sep = "\t", col.names = F, row.names = F)

f2.vcf <- read.vcfR("vcf_example_f2.new.vcf")


###########
# Using GQ
###########

# fazer teste com exemplos

data("vcf_example_f2")
segr <- test_segregation(vcf_example_f2)
plot(segr)
test.df <- create_probs(vcf_example_f2)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_riself")
segr <- test_segregation(vcf_example_riself)
plot(segr)

test.df <- create_probs(vcf_example_riself)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_out")
test.df <- create_probs(vcf_example_out)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_bc")
test.df <- create_probs(vcf_example_bc)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)


