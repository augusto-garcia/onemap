library(vcfR)

vcf.gatk <- read.vcfR("family1_gatk.vcf")
vcf.stacks <- read.vcfR("populations.snps.vcf")
vcf.freebayes <- read.vcfR("family1_freebayes.vcf")

gatk <- onemap_read_vcfR(vcfR.object = vcf.gatk, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")
str(gatk)

table(gatk$geno)

genotypes_errors <- extract_depth(vcfR.object = vcf.gatk, onemap.object = gatk, vcf.par = "GQ", parent1 = "P1", 
                                 parent2 = "P2", f1 = "F1", mean_phred = 20, recovering = FALSE)


onemap.obj <- gatk
old.errors <- create_probs(onemap.obj)

new.errors <- create_probs(onemap.obj, 
                           genotypes_errors = genotypes_errors, 
                           global_error = NULL, 
                           genotypes_probs = NULL)

twopts <- rf_2pts(old.errors)
twopts <- rf_2pts(new.errors)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg1  <- make_seq(lgs,1)
map.lg1.1 <- map(lg1)


old.updog <- updog_error(vcfR.object = vcf.gatk, onemap.object = gatk, vcf.par = "AD", parent1 = "P1", 
                         parent2 = "P2", f1 = "F1", recovering = TRUE, mean_phred = 20, cores = 3, 
                         depths = NULL)

# fazer teste com exemplos

data("vcf_example_f2")
test.df <- create_probs(vcf_example_f2)
twopts <- rf_2pts(test.df)
seq1 <- make_seq(twopts, "all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)
map.lg2 <- map(lg2)

data("vcf_example_riself")
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


