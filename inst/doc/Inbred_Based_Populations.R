## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
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


## ------------------------------------------------------------------------
library(onemap)

## ---- eval=FALSE---------------------------------------------------------
#  save.image("C:/.../yourfile.RData")

## ---- eval=FALSE---------------------------------------------------------
#  mapmaker_example_f2 <- read_mapmaker(dir="C:/workingdirectory",
#                                  file="your_data_file.raw")

## ---- eval=FALSE---------------------------------------------------------
#  mapmaker_example_f2 <- read_mapmaker(file= system.file("extdata/mapmaker_example_f2.raw",
#                                                         package = "onemap"))

## ---- load_data----------------------------------------------------------
data("mapmaker_example_f2")

## ------------------------------------------------------------------------
mapmaker_example_f2

## ---- eval=FALSE---------------------------------------------------------
#  onemap_example_f2 <- read_onemap(dir="C:/workingdirectory",
#                                  inputfile = "your_data_file.raw")

## ---- eval=FALSE---------------------------------------------------------
#  onemap_example_f2 <- read_onemap(inputfile= system.file("extdata/onemap_example_f2.raw",
#                                                         package = "onemap"))

## ------------------------------------------------------------------------
data("onemap_example_f2")

## ------------------------------------------------------------------------
onemap_example_f2

## ---- eval=FALSE---------------------------------------------------------
#  library(vcfR)
#  vcfR.object <- read.vcfR(system.file("extdata/vcf_example_f2.vcf", package = "onemap"))

## ---- eval=FALSE---------------------------------------------------------
#  vcf_example_f2 <- onemap_read_vcfR(vcfR.object = vcfR.object,
#                                     parent1 = "P1",
#                                     parent2 = "P2",
#                                     cross = "f2 intercross")

## ---- eval=FALSE---------------------------------------------------------
#  save(vcfR.object, file = "vcfR.object.RData")
#  rm(vcfR.object)

## ---- echo=FALSE---------------------------------------------------------
data(vcf_example_f2)

## ---- class_of_object----------------------------------------------------
class(onemap_example_f2)
class(vcf_example_f2)

## ---- plot_raw_data------------------------------------------------------
plot(onemap_example_f2)
plot(vcf_example_f2)

## ---- eval=FALSE---------------------------------------------------------
#  ?plot.onemap

## ---- plot_by_type-------------------------------------------------------
plot_by_segreg_type(onemap_example_f2)
plot_by_segreg_type(vcf_example_f2)

## ------------------------------------------------------------------------
comb_example <- combine_onemap(onemap_example_f2, vcf_example_f2)
comb_example

## ------------------------------------------------------------------------
plot(comb_example)

## ---- eval=FALSE---------------------------------------------------------
#  write_onemap_raw(comb_example, file.name = "new_dataset.raw", cross="f2 intercross")

## ---- chi_square---------------------------------------------------------
f2_test <- test_segregation(comb_example)

## ---- class_of_chisquare_test--------------------------------------------
class(f2_test)

## ---- print_chi1---------------------------------------------------------
f2_test

## ---- print_chi2---------------------------------------------------------
print(f2_test)

## ---- Bonferroni---------------------------------------------------------
Bonferroni_alpha(f2_test)

## ---- plot_chisq---------------------------------------------------------
plot(f2_test)

## ---- select_nondistorted------------------------------------------------
select_segreg(f2_test)

## ---- select_distorted---------------------------------------------------
select_segreg(f2_test, distorted = TRUE)

## ------------------------------------------------------------------------
no_dist <- select_segreg(f2_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

dist <- select_segreg(f2_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist


## ---- two_point_tests, results="hide"------------------------------------
twopts_f2 <- rf_2pts(comb_example)

## ---- Suggest_a_LOD------------------------------------------------------
(LOD_sug <- suggest_lod(comb_example))

## ---- print_2Mks_two_points----------------------------------------------
print(twopts_f2, c("M12", "M42"))

## ---- class_of_twopoint--------------------------------------------------
class(twopts_f2)

## ---- print_all_two_points-----------------------------------------------
print(twopts_f2)

## ---- subset_all---------------------------------------------------------
mark_all_f2 <- make_seq(twopts_f2, "all")

## ---- class_subset_all---------------------------------------------------
class(mark_all_f2)

## ---- subset_3mks--------------------------------------------------------
mrk_subset <- make_seq(twopts_f2, c(1, 3, 7))

## ---- without segregation distortion-------------------------------------
mark_no_dist_f2 <- make_seq(twopts_f2, no_dist)

## ---- group1-------------------------------------------------------------
LGs_f2 <- group(mark_all_f2)
LGs_f2

## ---- group2-------------------------------------------------------------
(LGs_f2 <- group(mark_all_f2, LOD = LOD_sug, max.rf = 0.5))

## ---- class_group--------------------------------------------------------
class(LGs_f2)

## ---- haldane, eval=FALSE------------------------------------------------
#  set_map_fun(type = "haldane")

## ---- kosambi, eval=FALSE------------------------------------------------
#  set_map_fun(type = "kosambi")

## ------------------------------------------------------------------------
LG2_f2 <- make_seq(LGs_f2, 2)

## ------------------------------------------------------------------------
LG2_f2

## ---- class_lg-----------------------------------------------------------
class(LG2_f2)

## ---- results="hide"-----------------------------------------------------
LG2_ser_f2 <- seriation(LG2_f2)
LG2_rcd_f2 <- rcd(LG2_f2)
LG2_rec_f2 <- record(LG2_f2)
LG2_ug_f2 <- ug(LG2_f2)

## ---- order_seq----------------------------------------------------------
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3)
                        

## ---- show_order_seq-----------------------------------------------------
LG2_f2_ord

## ---- safe, results="hide"-----------------------------------------------
LG2_f2_safe <- make_seq(LG2_f2_ord, "safe")

## ---- force--------------------------------------------------------------
(LG2_f2_all <- make_seq(LG2_f2_ord, "force"))

## ---- touchdown----------------------------------------------------------
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ---- lg2_final----------------------------------------------------------
(LG2_f2_final <- make_seq(LG2_f2_ord, "force"))

## ---- ripple_lg2_final, results="hide"-----------------------------------
ripple_seq(LG2_f2_final, ws = 5, LOD = 3)

## ------------------------------------------------------------------------
LG2_f2_final

## ------------------------------------------------------------------------
LG1_f2 <- make_seq(LGs_f2, 1)

## ---- results='hide'-----------------------------------------------------
LG1_f2_ord <- order_seq(input.seq = LG1_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ------------------------------------------------------------------------
(LG1_f2_final <- make_seq(LG1_f2_ord, "force"))

## ---- results="hide"-----------------------------------------------------
ripple_seq(LG1_f2_final, ws = 5)

## ------------------------------------------------------------------------
LG1_f2_final

## ------------------------------------------------------------------------
LG3_f2 <- make_seq(LGs_f2, 3)

## ---- order_LG3, results='hide'------------------------------------------
LG3_f2_ord <- order_seq(input.seq = LG3_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ---- LG3_force----------------------------------------------------------
(LG3_f2_final <- make_seq(LG3_f2_ord, "force"))

## ---- ripple_LG3, results="hide"-----------------------------------------
ripple_seq(LG3_f2_final, ws = 5)

## ---- LG3_final----------------------------------------------------------
LG3_f2_final

## ---- remove_M38, results="hide"-----------------------------------------
temp_seq <- drop_marker(LG3_f2_final, 38)

## ---- add_M38_end_LG-----------------------------------------------------
(temp_seq <- add_marker(temp_seq, 38))
(LG3_f2_wrong <- map(temp_seq))

## ---- rec_matrix---------------------------------------------------------
rf_graph_table(LG3_f2_wrong)

## ---- rec_matrix_inter, eval=FALSE---------------------------------------
#  rf_graph_table(LG3_f2_wrong, inter = TRUE, html.file = "LG3_f2_wrong.html")

## ------------------------------------------------------------------------
rf_graph_table(LG3_f2_wrong, n.colors = 7, main="LG3", lab.xy = c("markers", "markers"), mrk.axis = "numbers")

## ---- width=9, height=9--------------------------------------------------
temp_seq <- drop_marker(LG3_f2_wrong, 38)
temp_map <- map(temp_seq)
(temp_try <- try_seq(temp_map, 38))

## ------------------------------------------------------------------------
(LG3_f2_final <- make_seq(temp_try, 4))

## ------------------------------------------------------------------------
rf_graph_table(LG1_f2_final)

rf_graph_table(LG2_f2_final)

## ------------------------------------------------------------------------
CHR1 <- make_seq(twopts_f2, "1")
CHR1
CHR2 <- make_seq(twopts_f2, "2")
CHR3 <- make_seq(twopts_f2, "3")


## ------------------------------------------------------------------------
CHR_mks <- group_seq(input.2pts = twopts_f2, seqs = "CHROM", unlink.mks = mark_all_f2,
                      repeated = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  CHR_mks <- group_seq(input.2pts = twopts_f2, seqs = list(CHR1=CHR1, CHR2=CHR2, CHR3=CHR3),
#                        unlink.mks = mark_all_f2, repeated = FALSE)

## ------------------------------------------------------------------------
CHR_mks

## ------------------------------------------------------------------------
CHR_mks$repeated

## ------------------------------------------------------------------------
CHR_mks$sequences$CHR1
# or
CHR_mks$sequences[[1]]

## ---- results='hide'-----------------------------------------------------
CHR1_ord <- order_seq(CHR_mks$sequences$CHR1)
CHR1_frame <- make_seq(CHR1_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR1_frame) # graphic not showed

## ------------------------------------------------------------------------
CHR1_test_seq <- drop_marker(CHR1_frame, 58)
CHR1_test_map <- map(CHR1_test_seq)
CHR1_add58_seq <- try_seq(CHR1_test_map, 58)
CHR1_add58 <- make_seq(CHR1_add58_seq, 20) # marker 58 was placed at the same position as before

## ------------------------------------------------------------------------
CHR1_test_map

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR1_test_map) # graphic not showed

## ------------------------------------------------------------------------
CHR1_add58

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR1_add58)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR1_add58, inter = FALSE)

## ------------------------------------------------------------------------
CHR1_final <- CHR1_add58

## ------------------------------------------------------------------------
ripple_seq(CHR1_final)

## ---- results='hide'-----------------------------------------------------
CHR2_ord <- order_seq(CHR_mks$sequences$CHR2)
CHR2_frame <- make_seq(CHR2_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR2_frame) # graphic not showed

## ------------------------------------------------------------------------
CHR2_test_seq <- drop_marker(CHR2_frame, 20)
CHR2_test_map <- map(CHR2_test_seq)
CHR2_add20_seq <- try_seq(CHR2_test_map, 20)
CHR2_add20 <- make_seq(CHR2_add20_seq, 20) # marker 20 was placed at the same position as before

## ------------------------------------------------------------------------
CHR2_test_map

## ------------------------------------------------------------------------
rf_graph_table(CHR1_test_map)

## ------------------------------------------------------------------------
CHR2_add20

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR2_add20)  # graphic not showed

## ------------------------------------------------------------------------
CHR2_final <- CHR2_test_map

## ------------------------------------------------------------------------
ripple_seq(CHR2_final)

## ---- results='hide'-----------------------------------------------------
CHR3_ord <- order_seq(CHR_mks$sequences$CHR3)
CHR3_frame <- make_seq(CHR3_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR3_frame) # graphic not showed

## ---- results='hide'-----------------------------------------------------
CHR3_test_seq <- drop_marker(CHR3_frame, 32)
CHR3_test_ord <- order_seq(CHR3_test_seq)
CHR3_test_map <- make_seq(CHR3_test_ord, "force")

## ---- results='hide'-----------------------------------------------------
CHR3_add32_seq <- try_seq(CHR3_test_map, 32)
CHR3_add32 <- make_seq(CHR3_add32_seq, 13) # Marker 32 keeped in the map

## ------------------------------------------------------------------------
CHR3_final <- CHR3_add32
rf_graph_table(CHR3_final, inter = FALSE)

## ------------------------------------------------------------------------
ripple_seq(CHR3_final)

## ------------------------------------------------------------------------
LG3seq_f2 <- make_seq(twopts_f2, c(47, 38, 59, 16, 62, 21, 20, 48, 22))
(LG3seq_f2_map <- map(LG3seq_f2))

## ---- markers_names_and_numbers------------------------------------------
marker_type(LG3seq_f2_map)

## ---- add_marker---------------------------------------------------------
(LG3seq_f2_map <- add_marker(LG3seq_f2_map, c(18, 56, 50)))

## ---- drop_marker--------------------------------------------------------
(LG3seq_f2_map <- drop_marker(LG3seq_f2_map, c(59, 21)))

## ------------------------------------------------------------------------
maps_list <- list(LG1_f2_final, LG2_f2_final, LG3_f2_final)

## ------------------------------------------------------------------------
draw_map(maps_list, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ------------------------------------------------------------------------
draw_map(LG1_f2_final, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ------------------------------------------------------------------------
map_list_all <- list(LG1_f2_final, CHR1_final, LG2_f2_final, CHR3_final, CHR2_final, LG3_f2_final)
draw_map(map_list_all, names = TRUE, grid = TRUE, cex.mrk = 0.7)


## ---- eval=FALSE, results='hide', eval=FALSE-----------------------------
#  draw_map2(LG1_f2_final, LG2_f2_final, LG3_f2_final, main = "Only linkage information",
#            group.names = c("LG1", "LG2", "LG3"))

## ---- echo=FALSE, results='hide', echo=FALSE-----------------------------
draw_map2(LG1_f2_final, LG2_f2_final, LG3_f2_final, main = "Only linkage information", 
          group.names = c("LG1", "LG2", "LG3"), output = "map.png")

## ---- results='hide', eval=FALSE-----------------------------------------
#  draw_map2(LG1_f2_final, col.group = "#58A4B0", col.mark = "#335C81", output = "map_LG1.pdf")

## ---- results='hide', echo=FALSE-----------------------------------------
draw_map2(LG1_f2_final, col.group = "#58A4B0", col.mark = "#335C81", output = "map_LG1.png")

## ---- results='hide', eval=FALSE-----------------------------------------
#  draw_map2(LG1_f2_final, CHR1_final, LG2_f2_final, CHR3_final, CHR2_final, LG3_f2_final,
#            tag = c("M18", "M59", "M38", "M47", "M1"),
#            output = "map_all.pdf")
#  

## ---- results='hide', echo=FALSE-----------------------------------------
draw_map2(LG1_f2_final, CHR1_final, LG2_f2_final, CHR3_final, CHR2_final, LG3_f2_final,
          tag = c("M18", "M59", "M38", "M47", "M1"),
          output = "map_all.png")


## ------------------------------------------------------------------------
sessionInfo()

