context("Graphics utilities")

library(testthat)
library(vcfR)

test_that("create depth profile", {
  test_depth_profile <- function(df, cross, parent1, parent2, n.genos, n.genos.alt.ref){
    eval(bquote(onemap.obj <- onemap_read_vcfR(.(df), cross = .(cross), parent1 = .(parent1), parent2 = .(parent2))))
    temp.file <- tempfile()
    eval(bquote(p <- create_depths_profile(onemap.obj, vcf = .(df), 
                                           parent1 = .(parent1), parent2 = .(parent2), 
                                           rds.file = temp.file)))
    df <- readRDS(temp.file)
    file.remove(temp.file)
    eval(bquote(expect_equal(as.numeric(table(df$gt.onemap)), .(n.genos))))
    eval(bquote(expect_equal(as.numeric(table(df$gt.vcf)), .(n.genos))))
    eval(bquote(expect_equal(as.numeric(table(df$gt.onemap.alt.ref)), .(n.genos.alt.ref))))
    eval(bquote(expect_equal(as.numeric(table(df$gt.vcf.alt.ref)), .(n.genos.alt.ref))))
  }
  
  test_depth_profile(df = system.file("extdata/vcf_example_f2.vcf.gz", package = "onemap"),
                     cross = "f2 intercross", parent1 = "P1", parent2 = "P2",
                     n.genos = c(676, 987, 1937,1250), n.genos.alt.ref = c(1937, 1250, 987,676))
  
  test_depth_profile(df = system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"),
                     cross = "outcross", parent1 = "P1", parent2 = "P2", 
                     n.genos = c(4,774,1077,401), n.genos.alt.ref = c(1077,401,774,4))
})