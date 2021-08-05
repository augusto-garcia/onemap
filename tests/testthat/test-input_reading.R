context("Reading input files")
library(vcfR)

test_that("reading files",{
  expect_values_equal <- function(segr.type1.4, 
                                  segr.type.num1.4, n.phe, pheno1.3, dim.geno, table.geno, error1.4){
    
    eval(bquote(expect_equal(data$segr.type[1:4], .(segr.type1.4))))
    eval(bquote(expect_equal(data$segr.type.num[1:4], .(segr.type.num1.4))))
    eval(bquote(expect_equal(data$n.phe, .(n.phe))))
    eval(bquote(expect_equal(data$pheno[1:3,1], .(pheno1.3))))
    eval(bquote(expect_equal(dim(data$geno), .(dim.geno))))
    eval(bquote(expect_equal(as.vector(table(data$geno)), .(table.geno))))
    eval(bquote(expect_equal(as.vector(data$error[1:4,1]), .(error1.4))))
  }
  
  data <- read_mapmaker(system.file("extdata/mapmaker_example_f2.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = c("A.H.B", "C.A", "D.B", "C.A"),
                      segr.type.num1.4 = c(4,7,6,7),
                      n.phe = 1,
                      pheno1.3 = c(37.5892, 36.3664, 37.2230),
                      dim.geno = c(200,66), 
                      table.geno = c(1980, 5192,4338,1690),
                      error1.4 = rep(0.99999,4))
  
  data <- read_mapmaker(system.file("extdata/mapmaker_example_bc.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H",4),
                      segr.type.num1.4 = rep(8,4),
                      n.phe = 1,
                      pheno1.3 = c(40.7594, 39.5339, 37.9111),
                      dim.geno = c(150,67), 
                      table.geno = c(1507, 4411,4132),
                      error1.4 = c(rep(0.99999,3), 0.00001))
  
  data <- read_onemap(system.file("extdata/onemap_example_f2.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = c("A.H.B", "C.A", "D.B", "C.A"),
                      segr.type.num1.4 = c(4,7,6,7),
                      n.phe = 1,
                      pheno1.3 = c(37.5892, 36.3664, 37.2230),
                      dim.geno = c(200,66), 
                      table.geno = c(1980, 5192,4338,1690),
                      error1.4 = c(rep(0.99999,4)))
  
  data <- read_onemap(system.file("extdata/onemap_example_bc.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H",4),
                      segr.type.num1.4 = rep(8,4),
                      n.phe = 1,
                      pheno1.3 = c(40.7594, 39.5339, 37.9111),
                      dim.geno = c(150,67), 
                      table.geno = c(1507, 4411,4132),
                      error1.4 = c(rep(0.99999,3), 0.00001))
  
  data <- read_onemap(system.file("extdata/onemap_example_out.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = c("B3.7", "D2.18", "D1.13", "A.4"),
                      segr.type.num1.4 = c(4,7,6,1),
                      n.phe = 3,
                      pheno1.3 =  c(43, 12, 20),
                      dim.geno = c(100,30), 
                      table.geno = c(1221, 1182,390, 207),
                      error1.4 = c(rep(0.00001,3), 0.99999))
  
  
  data <- read_onemap(system.file("extdata/onemap_example_riself.raw", package = "onemap"))
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.B",4),
                      segr.type.num1.4 = rep(9,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(100,68), 
                      table.geno = c(597, 3229,2974),
                      error1.4 = c(0.00001, rep(0.99999,3)))
  
  vcfR.obj <- read.vcfR(system.file("extdata/vcf_example_bc.vcf", package = "onemap"))
  data <- onemap_read_vcfR(vcfR.obj, cross = "f2 backcross", parent1 = "P1", parent2 = "P2")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H",4),
                      segr.type.num1.4 = rep(8,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(142,25), 
                      table.geno = c(462, 1569,1519),
                      error1.4 = rep(10^(-5),4))
  
  vcfR.obj <- read.vcfR(system.file("extdata/vcf_example_f2.vcf", package = "onemap"))
  data <- onemap_read_vcfR(vcfR.obj, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H.B",4),
                      segr.type.num1.4 = rep(4,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(192,25), 
                      table.geno = c(676, 962,1937, 1225),
                      error1.4 = c(rep(0.00001,2),1, 0.99999))
  
  data <- onemap_read_vcfR(vcfR.obj, cross = "f2 intercross", parent1 = "P2", parent2 = "P1", f1 = "F1")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H.B",4),
                      segr.type.num1.4 = rep(4,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(192,25), 
                      table.geno = c(676, 1225,1937, 962),
                      error1.4 = c(rep(0.00001,2),1, 3.333333e-06))
  
  vcfR.obj <- read.vcfR(system.file("extdata/vcf_example_out.vcf", package = "onemap"))
  data <- onemap_read_vcfR(vcfR.object = vcfR.obj, cross = "outcross", parent1 = "P1", parent2 = "P2")
  expect_equal(check_data(data), 0)
  # expect_values_equal(segr.type1.4 = rep("B3.7",4),
  #                     segr.type.num1.4 = rep(4,4),
  #                     n.phe = 0,
  #                     pheno1.3 = NULL,
  #                     dim.geno =  c(92,24), 
  #                     table.geno = c(16, 761,1030, 401),
  #                     error1.4 = c(rep(0.99999,4)))
  
  data <- onemap_read_vcfR(vcfR.obj, cross = "outcross", parent1 = "P2", parent2 = "P1")
  expect_equal(check_data(data), 0)
  # expect_values_equal(segr.type1.4 = rep("B3.7",4),
  #                     segr.type.num1.4 = rep(4,4),
  #                     n.phe = 0,
  #                     pheno1.3 = NULL,
  #                     dim.geno =  c(92,24),
  #                     table.geno = c(16, 761,1030, 401),
  #                     error1.4 = c(rep(0.99999,4)))
  
  vcfR.obj <- read.vcfR(system.file("extdata/vcf_example_riself.vcf", package = "onemap"))
  data <- onemap_read_vcfR(vcfR.obj, cross = "ri self", parent1 = "P1", parent2 = "P2")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.B",4),
                      segr.type.num1.4 = rep(9,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(92,25), 
                      table.geno = c(87, 1092,1121),
                      error1.4 = c(rep(10^(-5),4)))
  
})

test_that("writting files", {
  expect_values_equal <- function(segr.type1.4, 
                                  segr.type.num1.4, n.phe, pheno1.3, dim.geno, table.geno, error1.4){
    
    eval(bquote(expect_equal(data$segr.type[1:4], .(segr.type1.4))))
    eval(bquote(expect_equal(data$segr.type.num[1:4], .(segr.type.num1.4))))
    eval(bquote(expect_equal(data$n.phe, .(n.phe))))
    eval(bquote(expect_equal(data$pheno[1:3,1], .(pheno1.3))))
    eval(bquote(expect_equal(dim(data$geno), .(dim.geno))))
    eval(bquote(expect_equal(as.vector(table(data$geno)), .(table.geno))))
    eval(bquote(expect_equal(as.vector(data$error[1:4,1]), .(error1.4))))
  }
  
  data(onemap_example_out)
  write_onemap_raw(onemap.obj = onemap_example_out, file.name = "test_out.raw")
  data <- read_onemap("test_out.raw")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = c("B3.7", "D2.18", "D1.13", "A.4"),
                      segr.type.num1.4 = c(4,7,6,1),
                      n.phe = 3,
                      pheno1.3 =  c(43, 12, 20),
                      dim.geno = c(100,30), 
                      table.geno = c(1221, 1182,390, 207),
                      error1.4 = c(rep(0.00001,3), 0.99999))
  
  file.remove("test_out.raw")
  
  data(onemap_example_f2)
  write_onemap_raw(onemap.obj = onemap_example_f2, file.name = "test_f2.raw")
  data <- read_onemap("test_f2.raw")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = c("A.H.B", "C.A", "D.B", "C.A"),
                      segr.type.num1.4 = c(4,7,6,7),
                      n.phe = 1,
                      pheno1.3 = c(37.5892, 36.3664, 37.2230),
                      dim.geno = c(200,66), 
                      table.geno = c(1980, 5192,4338,1690),
                      error1.4 = c(rep(0.99999,4)))
  
  file.remove("test_f2.raw")
  
  data("onemap_example_bc")
  write_onemap_raw(onemap.obj = onemap_example_bc, file.name = "test_bc.raw")
  data <- read_onemap("test_bc.raw")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.H",4),
                      segr.type.num1.4 = rep(8,4),
                      n.phe = 1,
                      pheno1.3 = c(40.7594, 39.5339, 37.9111),
                      dim.geno = c(150,67), 
                      table.geno = c(1507, 4411,4132),
                      error1.4 = c(rep(0.99999,3), 0.00001))
  
  file.remove("test_bc.raw")
  
  data("onemap_example_riself")
  write_onemap_raw(onemap.obj = onemap_example_riself, file.name = "test_rils.raw")
  data <- read_onemap("test_rils.raw")
  expect_equal(check_data(data), 0)
  expect_values_equal(segr.type1.4 = rep("A.B",4),
                      segr.type.num1.4 = rep(9,4),
                      n.phe = 0,
                      pheno1.3 = NULL,
                      dim.geno =  c(100,68), 
                      table.geno = c(597, 3229,2974),
                      error1.4 = c(0.00001, rep(0.99999,3)))
  file.remove("test_rils.raw")
})


