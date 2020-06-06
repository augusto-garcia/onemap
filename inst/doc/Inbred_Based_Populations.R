## ----knitr_init, echo=FALSE, cache=FALSE--------------------------------------
library(knitr)
library(rmarkdown)

knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      fig.width = 6,
                      fig.height = 6,
                      fig.align = "center",
                      dev = "png",
                      dpi = 36,
                      cache = TRUE)


## -----------------------------------------------------------------------------
library(onemap)

## ---- eval=FALSE--------------------------------------------------------------
#  save.image("C:/.../yourfile.RData")

## ---- eval=FALSE--------------------------------------------------------------
#  mapmaker_example_f2 <- read_mapmaker(dir="C:/workingdirectory",
#                                  file="your_data_file.raw")

## ---- eval=FALSE--------------------------------------------------------------
#  mapmaker_example_f2 <- read_mapmaker(file= system.file("extdata/mapmaker_example_f2.raw",
#                                                         package = "onemap"))

## ---- load_data---------------------------------------------------------------
data("mapmaker_example_f2")

## -----------------------------------------------------------------------------
mapmaker_example_f2

## ---- eval=FALSE--------------------------------------------------------------
#  onemap_example_f2 <- read_onemap(dir="C:/workingdirectory",
#                                  inputfile = "your_data_file.raw")

## ---- eval=FALSE--------------------------------------------------------------
#  onemap_example_f2 <- read_onemap(inputfile= system.file("extdata/onemap_example_f2.raw",
#                                                         package = "onemap"))

## -----------------------------------------------------------------------------
data("onemap_example_f2")

## -----------------------------------------------------------------------------
onemap_example_f2

## ---- eval=FALSE--------------------------------------------------------------
#  library(vcfR)
#  vcfR.object <- read.vcfR(system.file("extdata/vcf_example_f2.vcf", package = "onemap"))

## ---- eval=FALSE--------------------------------------------------------------
#  vcf_example_f2 <- onemap_read_vcfR(vcfR.object = vcfR.object,
#                                     parent1 = "P1",
#                                     parent2 = "P2",
#                                     cross = "f2 intercross")

## ---- eval=FALSE--------------------------------------------------------------
#  save(vcfR.object, file = "vcfR.object.RData")
#  rm(vcfR.object)

## ---- echo=FALSE--------------------------------------------------------------
data(vcf_example_f2)

## ---- class_of_object---------------------------------------------------------
class(onemap_example_f2)
class(vcf_example_f2)

## ---- plot_raw_data-----------------------------------------------------------
plot(onemap_example_f2)
plot(vcf_example_f2)

## ---- eval=FALSE--------------------------------------------------------------
#  ?plot.onemap

## ---- plot_by_type------------------------------------------------------------
plot_by_segreg_type(onemap_example_f2)
plot_by_segreg_type(vcf_example_f2)

## -----------------------------------------------------------------------------
comb_example <- combine_onemap(onemap_example_f2, vcf_example_f2)
comb_example

## -----------------------------------------------------------------------------
plot(comb_example)

## ---- eval=FALSE--------------------------------------------------------------
#  write_onemap_raw(comb_example, file.name = "new_dataset.raw", cross="f2 intercross")

## ---- chi_square--------------------------------------------------------------
f2_test <- test_segregation(comb_example)

## ---- class_of_chisquare_test-------------------------------------------------
class(f2_test)

## ---- print_chi1--------------------------------------------------------------
f2_test

## ---- print_chi2--------------------------------------------------------------
print(f2_test)

## ---- Bonferroni--------------------------------------------------------------
Bonferroni_alpha(f2_test)

## ---- plot_chisq--------------------------------------------------------------
plot(f2_test)

## ---- select_nondistorted-----------------------------------------------------
select_segreg(f2_test)

## ---- select_distorted--------------------------------------------------------
select_segreg(f2_test, distorted = TRUE)

## -----------------------------------------------------------------------------
no_dist <- select_segreg(f2_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

dist <- select_segreg(f2_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist


## ---- two_point_tests, results="hide"-----------------------------------------
twopts_f2 <- rf_2pts(comb_example)

## ---- Suggest_a_LOD-----------------------------------------------------------
(LOD_sug <- suggest_lod(comb_example))

## ---- print_2Mks_two_points---------------------------------------------------
print(twopts_f2, c("M12", "M42"))

## ---- class_of_twopoint-------------------------------------------------------
class(twopts_f2)

## ---- print_all_two_points----------------------------------------------------
print(twopts_f2)

## ---- subset_all--------------------------------------------------------------
mark_all_f2 <- make_seq(twopts_f2, "all")

## ---- class_subset_all--------------------------------------------------------
class(mark_all_f2)

## ---- subset_3mks-------------------------------------------------------------
mrk_subset <- make_seq(twopts_f2, c(1, 3, 7))

## ---- without segregation distortion------------------------------------------
mark_no_dist_f2 <- make_seq(twopts_f2, no_dist)

## ---- group1------------------------------------------------------------------
LGs_f2 <- group(mark_all_f2)
LGs_f2

## ---- group2------------------------------------------------------------------
(LGs_f2 <- group(mark_all_f2, LOD = LOD_sug, max.rf = 0.5))

## ---- class_group-------------------------------------------------------------
class(LGs_f2)

## ---- haldane, eval=FALSE-----------------------------------------------------
#  set_map_fun(type = "haldane")

## ---- kosambi, eval=FALSE-----------------------------------------------------
#  set_map_fun(type = "kosambi")

## -----------------------------------------------------------------------------
LG2_f2 <- make_seq(LGs_f2, 2)

## -----------------------------------------------------------------------------
LG2_f2

## ---- class_lg----------------------------------------------------------------
class(LG2_f2)

## ---- results="hide"----------------------------------------------------------
LG2_ser_f2 <- seriation(LG2_f2)
LG2_rcd_f2 <- rcd(LG2_f2)
LG2_rec_f2 <- record(LG2_f2)
LG2_ug_f2 <- ug(LG2_f2)

## ---- order_seq---------------------------------------------------------------
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3)
                        

## ---- show_order_seq----------------------------------------------------------
LG2_f2_ord

## ---- safe, results="hide"----------------------------------------------------
LG2_f2_safe <- make_seq(LG2_f2_ord, "safe")

## ---- force-------------------------------------------------------------------
(LG2_f2_all <- make_seq(LG2_f2_ord, "force"))

## ---- touchdown---------------------------------------------------------------
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ---- lg2_final---------------------------------------------------------------
(LG2_f2_final <- make_seq(LG2_f2_ord, "force"))

## ---- ripple_lg2_final, results="hide"----------------------------------------
ripple_seq(LG2_f2_final, ws = 5, LOD = 3)

## -----------------------------------------------------------------------------
LG2_f2_final

## -----------------------------------------------------------------------------
LG1_f2 <- make_seq(LGs_f2, 1)

## ---- results='hide'----------------------------------------------------------
LG1_f2_ord <- order_seq(input.seq = LG1_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

