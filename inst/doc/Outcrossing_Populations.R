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


## ---- echo=FALSE, results='hide'-----------------------------------------
library(onemap)

## ---- eval=FALSE---------------------------------------------------------
#  system.file(package = "onemap")

## ---- eval=FALSE---------------------------------------------------------
#  onemap_example_out <- read_onemap(dir = "C:/workingdirectory", inputfile = "onemap_example_out.raw")

## ---- eval=FALSE---------------------------------------------------------
#  onemap_example_out <- read_onemap(inputfile = "onemap_example_out.raw")

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  onemap_example_out <- read_onemap(inputfile = system.file("extdata/onemap_example_out.raw", package = "onemap"))

## ------------------------------------------------------------------------
data("onemap_example_out")

## ---- eval=FALSE---------------------------------------------------------
#  onemap_example_out

## ---- echo=FALSE---------------------------------------------------------
onemap_example_out

## ------------------------------------------------------------------------
plot(onemap_example_out)

## ------------------------------------------------------------------------
plot(onemap_example_out, all = FALSE)

## ------------------------------------------------------------------------
plot_by_segreg_type(onemap_example_out)

## ---- eval=FALSE---------------------------------------------------------
#  library(vcfR)
#  vcfR.object <- read.vcfR(system.file("extdata/vcf_example_out.vcf", package = "onemap"))

## ---- eval=FALSE---------------------------------------------------------
#  vcf_example_out <- onemap_read_vcfR(vcfR.object = vcfR.object,
#                                     parent1 = "P1",
#                                     parent2 = "P2",
#                                     cross = "outcross")

## ---- echo=FALSE---------------------------------------------------------
data(vcf_example_out)

## ---- eval=FALSE---------------------------------------------------------
#  save(vcfR.object, file = "vcfR.object.RData")
#  rm(vcfR.object)

## ------------------------------------------------------------------------
comb_example <- combine_onemap(onemap_example_out, vcf_example_out)
comb_example

## ------------------------------------------------------------------------
plot(comb_example)

## ------------------------------------------------------------------------
bins <- find_bins(comb_example, exact = FALSE)
bins

## ------------------------------------------------------------------------
bins_example <- create_data_bins(comb_example, bins)
bins_example

## ---- eval=FALSE---------------------------------------------------------
#  write_onemap_raw(comb_example, file.name = "new_dataset.raw", cross="outcross")

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
rf_graph_table(LG3_final)

## ------------------------------------------------------------------------
rf_graph_table(LG3_final, graph.LOD = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(LG3_final, inter = TRUE, html.file = "LG3.html")

## ------------------------------------------------------------------------
rf_graph_table(LG3_final, inter = FALSE, mrk.axis = "numbers")

## ------------------------------------------------------------------------
rf_graph_table(LG3_final, main = "LG3", inter = FALSE, n.colors = 7, lab.xy = c("markers", "markers"))

## ------------------------------------------------------------------------
LG2 <- make_seq(LGs, 2)
LG2

## ------------------------------------------------------------------------
LG2_rcd <- rcd(LG2)
LG2_rcd

## ---- , echo=TRUE--------------------------------------------------------
marker_type(LG2)

## ---- results='hide'-----------------------------------------------------
LG2_init <- make_seq(twopts, c(4, 20, 24, 49,50,51))

## ---- results='hide'-----------------------------------------------------
LG2_comp <- compare(LG2_init)

## ------------------------------------------------------------------------
LG2_frame <- make_seq(LG2_comp)

## ------------------------------------------------------------------------
rf_graph_table(LG2_frame, mrk.axis = "numbers")

## ---- results='hide'-----------------------------------------------------
LG2_extend <- try_seq(LG2_frame, 52)

## ------------------------------------------------------------------------
LG2_extend

## ------------------------------------------------------------------------
print(LG2_extend, 7)

## ---- echo=TRUE----------------------------------------------------------
LG2_test <- make_seq(LG2_extend, 7, 1)

## ------------------------------------------------------------------------
rf_graph_table(LG2_test, mrk.axis = "numbers")

## ------------------------------------------------------------------------
LG2_frame <- LG2_test

## ---- results='hide'-----------------------------------------------------
LG2_extend <- try_seq(LG2_frame, 9)
LG2_frame <- make_seq(LG2_extend, 3)
LG2_extend <- try_seq(LG2_frame, 16)
LG2_frame <- make_seq(LG2_extend, 1)
LG2_extend <- try_seq(LG2_frame, 21)
LG2_frame <- make_seq(LG2_extend, 4)
LG2_extend <- try_seq(LG2_frame, 23)
LG2_frame <- make_seq(LG2_extend, 5)
LG2_extend <- try_seq(LG2_frame, 27)
LG2_frame <- make_seq(LG2_extend, 1)
LG2_extend <- try_seq(LG2_frame, 29)
LG2_frame <- make_seq(LG2_extend, 12)
LG2_extend <- try_seq(LG2_frame, 45)
LG2_frame <- make_seq(LG2_extend, 8)
LG2_extend <- try_seq(LG2_frame, 46)
LG2_frame <- make_seq(LG2_extend, 7)
LG2_extend <- try_seq(LG2_frame, 47)
LG2_frame <- make_seq(LG2_extend, 7)
LG2_extend <- try_seq(LG2_frame, 48)
LG2_final <- make_seq(LG2_extend, 10)

## ------------------------------------------------------------------------
rf_graph_table(LG2_final)

## ------------------------------------------------------------------------
LG2_ord <- order_seq(LG2, n.init = 5, THRES = 3)

## ------------------------------------------------------------------------
LG2_ord

## ------------------------------------------------------------------------
LG2_safe <- make_seq(LG2_ord, "safe")

## ------------------------------------------------------------------------
LG2_all <- make_seq(LG2_ord, "force")
LG2_all

## ---- results='hide'-----------------------------------------------------
LG2_ord <- order_seq(LG2, n.init = 5, THRES = 3, touchdown = TRUE)

## ------------------------------------------------------------------------
LG2_ord

## ------------------------------------------------------------------------
rf_graph_table(LG2_all, mrk.axis = "numbers")

## ------------------------------------------------------------------------
ripple_seq(LG2_all, ws = 4, LOD = LOD_sug)

## ------------------------------------------------------------------------
LG2_test_seq <- drop_marker(LG2_all, c(23,29))

## ------------------------------------------------------------------------
(LG2_test_map <- map(LG2_test_seq))

## ------------------------------------------------------------------------
rf_graph_table(LG2_test_map, mrk.axis = "numbers")

## ------------------------------------------------------------------------
LG2_test_seq <- try_seq(LG2_test_map, 23)
LG2_test_23 <- make_seq(LG2_test_seq, 6)
LG2_test_seq <- try_seq(LG2_test_23, 29)
LG2_test_23_29 <- make_seq(LG2_test_seq, 12)

rf_graph_table(LG2_test_23_29, mrk.axis = "numbers")

## ------------------------------------------------------------------------
LG2_final <- LG2_test_23
LG2_final

rf_graph_table(LG2_final, mrk.axis = "numbers")

## ------------------------------------------------------------------------
LG1 <- make_seq(LGs, 1)

## ---- results='hide'-----------------------------------------------------
LG1_ord <- order_seq(LG1, n.init = 6, touchdown = TRUE)

## ------------------------------------------------------------------------
LG1_ord

## ------------------------------------------------------------------------
(LG1_frame <- make_seq(LG1_ord, "force"))

## ------------------------------------------------------------------------
rf_graph_table(LG1_frame, mrk.axis = "numbers")

## ------------------------------------------------------------------------
ripple_seq(LG1_frame)

## ------------------------------------------------------------------------
LG1_test_seq <- drop_marker(LG1_frame, c(10,11,28,42))
LG1_test_map <- map(LG1_test_seq)

## ------------------------------------------------------------------------
rf_graph_table(LG1_test_map, mrk.axis = "numbers")

## ---- results='hide'-----------------------------------------------------
LG1_extend <- try_seq(LG1_test_map,10)
LG1_test_map <- make_seq(LG1_extend,15) 
LG1_extend <- try_seq(LG1_test_map,11)
LG1_test <- make_seq(LG1_extend,23) # We choose to remove this marker
LG1_extend <- try_seq(LG1_test_map,28)
LG1_test_map <- make_seq(LG1_extend,16) 
LG1_extend <- try_seq(LG1_test_map,42)
LG1_final <- make_seq(LG1_extend,17) 

## ------------------------------------------------------------------------
LG1_final

## ------------------------------------------------------------------------
rf_graph_table(LG1_final)

## ---- eval=FALSE---------------------------------------------------------
#  LG1_ser <- seriation(LG1)
#  LG1_rcd <- rcd(LG1)
#  LG1_rec <- record(LG1)
#  LG1_ug  <- ug(LG1)

## ------------------------------------------------------------------------
CHR1 <- make_seq(twopts, "1")
CHR1
CHR2 <- make_seq(twopts, "2")
CHR3 <- make_seq(twopts, "3")

## ------------------------------------------------------------------------
CHR_mks <- group_seq(input.2pts = twopts, seqs = "CHROM", unlink.mks = mark_no_dist,
                      repeated = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  CHR_mks <- group_seq(input.2pts = twopts, seqs = list(CHR1=CHR1, CHR2=CHR2, CHR3=CHR3),
#                        unlink.mks = mark_no_dist, repeated = FALSE)

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
CHR1_test_seq <- drop_marker(CHR1_frame, 11)
CHR1_test_map <- map(CHR1_test_seq)
CHR1_add11_seq <- try_seq(CHR1_test_map, 11)
CHR1_add11 <- make_seq(CHR1_add11_seq, 25) # marker 11 was placed at the same position as before

## ------------------------------------------------------------------------
CHR1_test_map

## ------------------------------------------------------------------------
rf_graph_table(CHR1_test_map) 

## ------------------------------------------------------------------------
CHR1_final <- CHR1_test_map

## ------------------------------------------------------------------------
ripple_seq(CHR1_final)

## ---- results='hide'-----------------------------------------------------
CHR2_ord <- order_seq(CHR_mks$sequences$CHR2)
CHR2_frame <- make_seq(CHR2_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR2_frame) # graphic not showed

## ------------------------------------------------------------------------
CHR2_final <- CHR2_frame

## ---- results='hide'-----------------------------------------------------
CHR3_ord <- order_seq(CHR_mks$sequences$CHR3)
CHR3_frame <- make_seq(CHR3_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR3_frame, mrk.axis = "numbers") # graphic not showed

## ---- results='hide'-----------------------------------------------------
CHR3_test_seq <- drop_marker(CHR3_frame, c(29))
CHR3_test_ord <- order_seq(CHR3_test_seq)
CHR3_test_map <- make_seq(CHR3_test_ord, "force")

## ---- eval=FALSE---------------------------------------------------------
#  rf_graph_table(CHR3_test_map, mrk.axis = "numbers") #graphic not showed

## ---- results='hide'-----------------------------------------------------
CHR3_add29_seq <- try_seq(CHR3_test_map, 29)
CHR3_add29 <- make_seq(CHR3_add29_seq, 12) # Marker 29 increase the map size disproportionately, it was removed from the map

## ------------------------------------------------------------------------
CHR3_final <- CHR3_test_map
rf_graph_table(CHR3_final, inter = FALSE)

## ------------------------------------------------------------------------
ripple_seq(CHR3_final)

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  data("example.RData")

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  segr <- test_segregation(example)
#  plot(segr) # all markers follow the expected Mendelian segregation pattern
#  
#  twopts <- rf_2pts(example)
#  seq <- make_seq(twopts, "all")
#  
#  lgs <- group(seq)
#  lg1 <- make_seq(lgs,1)
#  

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  #  MDS method
#  start <- Sys.time()
#  map.mds3 <- mds_onemap(lg1, out.file = "mds.input.file.txt", displaytext = F, mds.graph.file = "graph.mds.pdf", hmm = T)
#  end <- Sys.time()
#  end - start # It will take about 25 minutes
#  rf_graph_table(map.mds)
#  
#  # order_seq method
#  start <- Sys.time()
#  map.order <- order_seq(lg1)
#  end <- Sys.time()
#  end - start # It will take about 25 mins
#  map.seq <- make_seq(map.order, "force")
#  rf_graph_table(map.seq)
#  
#  # real order
#  real.ord <- make_seq(twopts, lg1$seq.num)
#  start <- Sys.time()
#  real.map <- map(real.ord)
#  end <- Sys.time()
#  end - start # It will take about 25 mins
#  rf_graph_table(real.map)

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

## ---- echo=TRUE, results='hide', eval=FALSE------------------------------
#  draw_map2(LG1_final, LG2_final, LG3_final, main = "Only with linkage information",
#            group.names = c("LG1", "LG2", "LG3"))

## ---- echo=FALSE, results='hide'-----------------------------------------
draw_map2(LG1_final, LG2_final, LG3_final, main = "Only with linkage information", 
          group.names = c("LG1", "LG2", "LG3"), output = "map1.png")

## ---- echo=TRUE, results='hide', eval=FALSE------------------------------
#  draw_map2(CHR1_final, CHR2_final, CHR3_final, output= "map_ref.pdf",
#            col.group = "#58A4B0",
#           col.mark= "#335C81")

## ---- echo=TRUE, results='hide', echo=FALSE------------------------------
draw_map2(CHR1_final, CHR2_final, CHR3_final, output= "map_ref.png", 
          col.group = "#58A4B0",
         col.mark= "#335C81")

## ---- results='hide', eval=FALSE-----------------------------------------
#  draw_map2(LG1_final, CHR1_final, output = "map_comp.pdf", tag = c("M1","SNP2"))

## ---- results='hide', echo=FALSE-----------------------------------------
draw_map2(LG1_final, CHR1_final, output = "map_comp.png", tag = c("M1","SNP2"))

## ---- results='hide', eval=FALSE-----------------------------------------
#  draw_map2(LG2_final, CHR3_final, tag= c("SNP17", "SNP18", "M29"), main = "Chromosome 3",
#            group.names = c("Only linkage", "With genome"), centered = TRUE, output = "map_comp2.pdf")

## ---- results='hide', echo=FALSE-----------------------------------------
draw_map2(LG2_final, CHR3_final, tag= c("SNP17", "SNP18", "M29"), main = "Chromosome 3", 
          group.names = c("Only linkage", "With genome"), centered = TRUE, output = "map_comp2.png")

## ------------------------------------------------------------------------
any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2))
(any_seq_map <- map(any_seq))

## ---- eval=FALSE---------------------------------------------------------
#  any_seq <- make_seq(twopts, c(30, 12, 3, 14, 2), phase = c(4, 1, 4, 3))
#  (any_seq_map <- map(any_seq))

## ------------------------------------------------------------------------
(any_seq <- add_marker(any_seq, 4:8))

## ------------------------------------------------------------------------
(any_seq <- drop_marker(any_seq, c(3, 4, 5, 12, 30)))

## ------------------------------------------------------------------------
sessionInfo()

