## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)
library(rmdformats)

## Global options
options(max.print="75")
opts_knit$set(width=75)

## ---- globalsetup, echo=FALSE, results='hide', cache=FALSE---------------
#opts_chunk$set(cache=TRUE, autodep=TRUE)

## ---- echo=FALSE, results='hide'-----------------------------------------
library(onemap)

## ---- eval=FALSE---------------------------------------------------------
#  system.file(package = "onemap")

## ---- eval=FALSE, message='hide'-----------------------------------------
#  vcf2raw(input = system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"),
#          output = "vcf_example_out.raw", parent1 = "P1", parent2 = "P2", cross = "outcross")

## ---- eval=FALSE---------------------------------------------------------
#  example_out <- read_onemap(dir = "C:/workingdirectory", inputfile = "example_out.raw")

## ---- eval=FALSE---------------------------------------------------------
#  example_out <- read_onemap(inputfile = "example_out.raw")

## ---- echo=FALSE---------------------------------------------------------
example_out <- read_onemap(inputfile = system.file("extdata/example_out.raw", package = "onemap"))

## ------------------------------------------------------------------------
data(example_out)

## ---- eval=FALSE---------------------------------------------------------
#  example_out

## ---- echo=FALSE---------------------------------------------------------
example_out

## ------------------------------------------------------------------------
plot(example_out)

## ------------------------------------------------------------------------
plot(example_out, all = FALSE)

## ------------------------------------------------------------------------
plot_by_segreg_type(example_out)

## ---- eval=FALSE---------------------------------------------------------
#  vcf_example_out <- read_onemap("C:/workingdirectory", "vcf_example_out.raw")

## ------------------------------------------------------------------------
data(vcf_example_out)

## ------------------------------------------------------------------------
vcf_example_out

## ------------------------------------------------------------------------
plot(vcf_example_out, all = FALSE)

## ------------------------------------------------------------------------
comb_example <- combine_onemap(example_out, vcf_example_out)
comb_example

## ------------------------------------------------------------------------
plot(comb_example)

## ------------------------------------------------------------------------
bins <- find_bins(comb_example, exact = FALSE)
bins

## ------------------------------------------------------------------------
bins_example <- create_data_bins(comb_example, bins)
bins_example

## ------------------------------------------------------------------------
test_segregation_of_a_marker(bins_example, 4)

## ------------------------------------------------------------------------
segreg_test <- test_segregation(bins_example)
print(segreg_test)

## ------------------------------------------------------------------------
select_segreg(segreg_test, distorted = TRUE) #to show the markers names with segregation distortion

select_segreg(segreg_test, distorted = FALSE) #to show the markers names without segregation distortion

## ------------------------------------------------------------------------
dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist

no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

## ------------------------------------------------------------------------
plot(segreg_test)

## ------------------------------------------------------------------------
LOD_sug <- suggest_lod(bins_example)
LOD_sug

## ---- twopoint, results='hide'-------------------------------------------
twopts <- rf_2pts(bins_example, LOD = LOD_sug, max.rf = 0.5)

## ---- eval=FALSE---------------------------------------------------------
#  twopts

## ---- echo=FALSE---------------------------------------------------------
twopts

## ------------------------------------------------------------------------
print(twopts, c("M1", "M3"))

## ------------------------------------------------------------------------
mark_no_dist <- make_seq(twopts, c(no_dist))

## ------------------------------------------------------------------------
marker_type(mark_no_dist)

## ------------------------------------------------------------------------
LGs <- group(mark_no_dist)

## ---- eval=FALSE---------------------------------------------------------
#  LGs

## ---- echo=FALSE---------------------------------------------------------
LGs

## ------------------------------------------------------------------------
print(LGs, detailed = FALSE)

## ------------------------------------------------------------------------
LGs <- group(mark_no_dist, LOD = 6)
LGs

## ------------------------------------------------------------------------
LGs <- group(mark_no_dist, LOD = LOD_sug, max.rf = 0.4)
LGs

## ---- eval=FALSE---------------------------------------------------------
#  set_map_fun(type = "haldane")

## ---- eval=FALSE---------------------------------------------------------
#  set_map_fun(type = "kosambi")

## ------------------------------------------------------------------------
LG3 <- make_seq(LGs, 3)

## ---- eval=FALSE---------------------------------------------------------
#  LG3

## ---- echo=FALSE---------------------------------------------------------
LG3

## ----  results='hide'----------------------------------------------------
LG3_ser <- seriation(LG3)
LG3_rcd <- rcd(LG3)
LG3_rec <- record(LG3)
LG3_ug <- ug(LG3)

## ---- eval=TRUE,  results='hide'-----------------------------------------
LG3_comp <- compare(LG3)

## ------------------------------------------------------------------------
LG3_comp

## ------------------------------------------------------------------------
LG3_final <- make_seq(LG3_comp, 1, 1)

## ------------------------------------------------------------------------
LG3_final <- make_seq(LG3_comp)

## ------------------------------------------------------------------------
LG3_final

## ------------------------------------------------------------------------
rf_graph_table(LG3_final, inter = FALSE)

## ------------------------------------------------------------------------
LG2 <- make_seq(LGs, 2)
LG2

## ------------------------------------------------------------------------
LG2_rcd <- rcd(LG2)
LG2_rcd

## ---- , echo=TRUE--------------------------------------------------------
marker_type(LG2)

## ---- results='hide'-----------------------------------------------------
LG2_init <- make_seq(twopts, c(4, 19, 23, 48,49,50))

## ---- results='hide'-----------------------------------------------------
LG2_comp <- compare(LG2_init)

## ------------------------------------------------------------------------
LG2_frame <- make_seq(LG2_comp)

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG2_frame)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG2_frame, inter = FALSE)

## ---- results='hide'-----------------------------------------------------
LG2_extend <- try_seq(LG2_frame, 51)

## ------------------------------------------------------------------------
LG2_extend

## ------------------------------------------------------------------------
print(LG2_extend, 7)

## ---- echo=TRUE----------------------------------------------------------
LG2_test <- make_seq(LG2_extend, 7, 1)

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG2_test)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG2_test, inter=FALSE)

## ------------------------------------------------------------------------
LG2_frame <- LG2_test

## ---- results='hide'-----------------------------------------------------
LG2_extend <- try_seq(LG2_frame, 8)
LG2_frame <- make_seq(LG2_extend, 3)
LG2_extend <- try_seq(LG2_frame, 15)
LG2_frame <- make_seq(LG2_extend, 1)
LG2_extend <- try_seq(LG2_frame, 20)
LG2_frame <- make_seq(LG2_extend, 4)
LG2_extend <- try_seq(LG2_frame, 22)
LG2_frame <- make_seq(LG2_extend, 5)
LG2_extend <- try_seq(LG2_frame, 26)
LG2_frame <- make_seq(LG2_extend, 1)
LG2_extend <- try_seq(LG2_frame, 28)
LG2_frame <- make_seq(LG2_extend, 12)
LG2_extend <- try_seq(LG2_frame, 44)
LG2_frame <- make_seq(LG2_extend, 7)
LG2_extend <- try_seq(LG2_frame, 45)
LG2_frame <- make_seq(LG2_extend, 6)
LG2_extend <- try_seq(LG2_frame, 46)
LG2_frame <- make_seq(LG2_extend, 6)
LG2_extend <- try_seq(LG2_frame, 47)
LG2_final <- make_seq(LG2_extend, 6)

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG2_final)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG2_final, inter=FALSE)

## ------------------------------------------------------------------------
LG2_ord <- order_seq(LG2, n.init = 5, THRES = 3)

## ------------------------------------------------------------------------
LG2_ord

## ------------------------------------------------------------------------
LG2_safe <- make_seq(LG2_ord, "safe")

## ------------------------------------------------------------------------
LG2_all <- make_seq(LG2_ord, "force")
LG2_all

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG2_all)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG2_all, inter = FALSE)

## ---- results='hide'-----------------------------------------------------
LG2_ord <- order_seq(LG2, n.init = 5, THRES = 3, touchdown = TRUE)

## ------------------------------------------------------------------------
LG2_ord

## ------------------------------------------------------------------------
ripple_seq(LG2_all, ws = 4, LOD = LOD_sug)

## ------------------------------------------------------------------------
LG2_test_seq <- drop_marker(LG2_all, c(22,50))

## ---- eval=FALSE---------------------------------------------------------
#  (LG2_test_map <- map(LG2_test_seq))

## ---- echo=FALSE---------------------------------------------------------
(LG2_test_map <- onemap::map(LG2_test_seq))

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG2_test_map)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG2_test_map, inter=FALSE)

## ------------------------------------------------------------------------
LG2_test_seq <- try_seq(LG2_test_map, 50)
LG2_test_50 <- make_seq(LG2_test_seq, 15)
LG2_test_seq <- try_seq(LG2_test_map, 22)
LG2_test_50_22 <- make_seq(LG2_test_seq, 6)

## ------------------------------------------------------------------------
LG2_final <- LG2_test_50
LG2_final

## ------------------------------------------------------------------------
LG1 <- make_seq(LGs, 1)

## ---- results='hide'-----------------------------------------------------
LG1_ord <- order_seq(LG1, n.init = 6, touchdown = TRUE)

## ------------------------------------------------------------------------
LG1_ord

## ------------------------------------------------------------------------
(LG1_frame <- make_seq(LG1_ord, "force"))

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG1_frame)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG1_frame, inter = FALSE)

## ------------------------------------------------------------------------
ripple_seq(LG1_frame)

## ---- eval=FALSE---------------------------------------------------------
#  LG1_test_seq <- drop_marker(LG1_frame, c(9,10,27,41))
#  LG1_test_map <- map(LG1_test_seq)

## ---- echo=FALSE---------------------------------------------------------
LG1_test_seq <- drop_marker(LG1_frame, c(9,10,27,41))
LG1_test_map <- onemap::map(LG1_test_seq)

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG1_test_map)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG1_test_map, inter=FALSE)

## ---- results='hide'-----------------------------------------------------
LG1_extend <- try_seq(LG1_test_map,9)
LG1_test <- make_seq(LG1_extend,15) # We choose to remove this marker
LG1_extend <- try_seq(LG1_test_map,10)
LG1_test <- make_seq(LG1_extend,22) # We choose to remove this marker
LG1_extend <- try_seq(LG1_test_map,27)
LG1_test <- make_seq(LG1_extend,15) 
LG1_extend <- try_seq(LG1_test,41)
LG1_final <- make_seq(LG1_extend,17) 

## ------------------------------------------------------------------------
LG1_final

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG1_final)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(LG1_final, inter=FALSE)

## ------------------------------------------------------------------------
LG1_ser <- seriation(LG1)
LG1_rcd <- rcd(LG1)
LG1_rec <- record(LG1)
LG1_ug  <- ug(LG1)

## ------------------------------------------------------------------------
CHR1 <- make_seq(twopts, "1")
CHR1
CHR2 <- make_seq(twopts, "2")
CHR3 <- make_seq(twopts, "3")


## ------------------------------------------------------------------------
CHR_mks <- group_seq(input.2pts = twopts, seqs = "CHROM", unlink.mks = mark_no_dist,
                      rm.repeated = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  CHR_mks <- group_seq(input.2pts = twopts, seqs = list(CHR1=CHR1, CHR2=CHR2, CHR3=CHR3),
#                        unlink.mks = mark_no_dist, rm.repeated = TRUE)

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
#  rf_graph_table(CHR1_frame)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR1_frame, inter=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  CHR1_test_seq <- drop_marker(CHR1_frame, 10)
#  CHR1_test_map <- map(CHR1_test_seq)
#  CHR1_add10_seq <- try_seq(CHR1_test_map, 10)
#  CHR1_add10 <- make_seq(CHR1_add10_seq, 25) # marker 10 was placed at the same position as before

## ---- echo=FALSE, results='hide'-----------------------------------------
CHR1_test_seq <- drop_marker(CHR1_frame, 10)
CHR1_test_map <- onemap::map(CHR1_test_seq)
CHR1_add10_seq <- try_seq(CHR1_test_map, 10)
CHR1_add10 <- make_seq(CHR1_add10_seq, 25) # marker 10 was placed at the same position as before

## ------------------------------------------------------------------------
CHR1_test_map

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR1_test_map)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR1_test_map, inter = FALSE)

## ------------------------------------------------------------------------
CHR1_add10

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR1_add10)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR1_add10, inter = FALSE)

## ------------------------------------------------------------------------
CHR1_final <- CHR1_test_map

## ------------------------------------------------------------------------
ripple_seq(CHR1_final)

## ---- results='hide'-----------------------------------------------------
CHR2_ord <- order_seq(CHR_mks$sequences$CHR2)
CHR2_frame <- make_seq(CHR2_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR2_frame)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR2_frame, inter=FALSE)

## ------------------------------------------------------------------------
CHR2_final <- CHR2_frame

## ---- results='hide'-----------------------------------------------------
CHR3_ord <- order_seq(CHR_mks$sequences$CHR3)
CHR3_frame <- make_seq(CHR3_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR3_frame)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR3_frame, inter=FALSE)

## ---- results='hide'-----------------------------------------------------
CHR3_test_seq <- drop_marker(CHR3_frame, c(22,26, 28))
CHR3_test_ord <- order_seq(CHR3_test_seq)
CHR3_test_map <- make_seq(CHR3_test_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR3_test_map)

## ---- echo=FALSE---------------------------------------------------------
rf_graph_table(CHR3_test_map, inter = FALSE)

## ---- results='hide'-----------------------------------------------------
CHR3_add22_seq <- try_seq(CHR3_test_map, 22)
CHR3_add22 <- make_seq(CHR3_add22_seq, 6) # Marker 22 removed of the map
CHR3_add26_seq <- try_seq(CHR3_test_map, 26)
CHR3_add26 <- make_seq(CHR3_add26_seq, 1) # Marker 26 was better positioned
CHR3_add28_seq <- try_seq(CHR3_add26, 28)
CHR3_add28 <- make_seq(CHR3_add28_seq, 12) # Marker 28 increase the map size disproportionately, it was removed from the map

## ------------------------------------------------------------------------
CHR3_final <- CHR3_add26
rf_graph_table(CHR3_final, inter = FALSE)

## ------------------------------------------------------------------------
ripple_seq(CHR3_final)

## ---- echo=TRUE, fig=TRUE------------------------------------------------
map1 <- list(LG1_final, LG2_final, LG3_final)
draw_map(map1, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- echo=TRUE, fig=TRUE------------------------------------------------
map2 <- list(CHR1_final, CHR2_final, CHR3_final)
draw_map(map2, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- width = 9.5, height = 9.5------------------------------------------
CHR1_comp <- list(LG1_final, CHR1_final)
draw_map(CHR1_comp, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- width = 9.5, height = 9.5------------------------------------------
CHR2_comp <- list(LG3_final, CHR2_final)
draw_map(CHR2_comp, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- width = 9.5, height = 9.5------------------------------------------
CHR3_comp <- list(LG2_final, CHR3_final)
draw_map(CHR3_comp, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ------------------------------------------------------------------------
draw_map(CHR1_final, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- eval=FALSE---------------------------------------------------------
#  any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2))
#  (any_seq_map <- map(any_seq))

## ---- echo=FALSE---------------------------------------------------------
any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2))
(any_seq_map <- onemap::map(any_seq))

## ---- eval=FALSE---------------------------------------------------------
#  any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2), phase = c(4, 1, 4, 3))
#  (any_seq_map <- map(any_seq))

## ---- echo=FALSE---------------------------------------------------------
any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2), phase = c(4, 1, 4, 3))
(any_seq_map <- onemap::map(any_seq))

## ------------------------------------------------------------------------
(any_seq <- add_marker(any_seq, 4:8))

## ------------------------------------------------------------------------
(any_seq <- drop_marker(any_seq, c(3, 4, 5, 12, 30)))

## ------------------------------------------------------------------------
sessionInfo()

