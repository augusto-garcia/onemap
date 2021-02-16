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


## ---- warning=FALSE, message=FALSE--------------------------------------------
library(onemap)

## ---- eval=FALSE--------------------------------------------------------------
#  save.image("C:/.../yourfile.RData")

## ---- eval= FALSE, message=F,warning=F----------------------------------------
#  # For f2 intercross population
#  run_pedsim(chromosome = c("Chr1", "Chr10"), n.marker = c(24,15),
#             tot.size.cm = c(100,150), centromere = c(50, 75),
#             n.ind = 200, mk.types = c("A.H.B", "C.A", "D.B"),
#             n.types = c(10,14,15), pop = "F2", path.pedsim = "~/onemap/",
#             name.mapfile = "mapfile.txt", name.founderfile="founderfile.gen",
#             name.chromfile="sim.chrom", name.parfile="sim.par",
#             name.out="sim_F2")

## ---- eval=FALSE, message=F,warning=F-----------------------------------------
#  # For f2 intercross population
#  pedsim2raw(cross="f2 intercross", genofile = "sim_F2_genotypes.dat",
#             parent1 = "P1", parent2 = "P2", f1 = "F1",
#             out.file = "sim_F2.example.raw", miss.perc = 10)

## ---- eval=FALSE, message=F,warning=F-----------------------------------------
#  # Dominant markers are not supported, then, we simulate other dataset with only codominant markers
#  run_pedsim(chromosome = c("Chr1", "Chr10"), n.marker = c(15,15),
#             tot.size.cm = c(100,150), centromere = c(50, 75),
#             n.ind = 200, mk.types = c("A.H.B"),
#             n.types = c(30), pop = "F2", path.pedsim = "~/onemap/",
#             name.mapfile = "mapfile_f2.txt", name.founderfile="founderfile.gen",
#             name.chromfile="sim_f2.chrom", name.parfile="sim.par",
#             name.out="sim_cod_F2")

## ---- eval=FALSE, message=F,warning=F-----------------------------------------
#  # For outcrossing population
#  pedsim2vcf(inputfile = "sim_cod_F2_genotypes.dat",
#             map.file = "mapfile_f2.txt",
#             chrom.file = "sim_f2.chrom",
#             out.file = "simu_cod_F2.vcf",
#             miss.perc = 0,
#             counts = TRUE,
#             mean.depth = 100,
#             p.mean.depth = 100,
#             chr.mb = 10,
#             method = "updog",
#             mean.phred = 20,
#             bias=1,
#             od=0.00001,
#             pos=NULL,
#             chr=NULL,
#             phase = FALSE,
#             disper.par=2)

## ---- message=F,warning=F, echo=FALSE-----------------------------------------
# For outcrossing population
pedsim2vcf(inputfile = system.file("extdata/sim_cod_F2_genotypes.dat", package = "onemap"),
           map.file = system.file("extdata/mapfile_f2.txt", package = "onemap"), 
           chrom.file = system.file("extdata/sim_f2.chrom", package = "onemap"), 
           out.file = "simu_cod_F2.vcf",
           miss.perc = 0, 
           counts = TRUE, 
           mean.depth = 100, 
           p.mean.depth = 100, 
           chr.mb = 10, 
           method = "updog", 
           mean.phred = 20,
           bias=1, 
           od=0.00001,
           pos=NULL,
           chr=NULL,
           phase = FALSE,
           disper.par=2)

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

## -----------------------------------------------------------------------------
library(vcfR)
vcfR.object <- read.vcfR(system.file("extdata/vcf_example_f2.vcf", package = "onemap"))

## ---- eval=FALSE--------------------------------------------------------------
#  vcf_example_f2 <- onemap_read_vcfR(vcfR.object = vcfR.object,
#                                     parent1 = "P1",
#                                     parent2 = "P2",
#                                     cross = "f2 intercross")

## ---- eval=FALSE--------------------------------------------------------------
#  save(vcfR.object, file = "vcfR.object.RData")
#  #rm(vcfR.object)

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

## ---- message=F,warning=F-----------------------------------------------------
vcf_simu_f2 <- read.vcfR(system.file("extdata/simu_cod_F2.vcf", package="onemap"))
simu_f2_obj <- onemap_read_vcfR(vcfR.object = vcf_simu_f2, cross = "f2 intercross", parent1 = "P1", parent2 = "P2", f1 = "F1")

create_depths_profile(onemap.obj = simu_f2_obj,
                      vcfR.object = vcf_simu_f2, 
                      parent1 = "P1", 
                      parent2 = "P2", 
                      f1 = "F1",
                      vcf.par = "AD", 
                      recovering = FALSE, 
                      mks = NULL, 
                      inds = NULL, 
                      GTfrom = "vcf", 
                      alpha = 0.1,
                      rds.file = "depths_f2.rds")

## -----------------------------------------------------------------------------
comb_example <- combine_onemap(onemap_example_f2, vcf_example_f2)
comb_example

## -----------------------------------------------------------------------------
plot(comb_example)

## -----------------------------------------------------------------------------
bins <- find_bins(comb_example, exact = FALSE)
bins

## -----------------------------------------------------------------------------
bins_example <- create_data_bins(comb_example, bins)
bins_example

## ---- eval=FALSE--------------------------------------------------------------
#  write_onemap_raw(bins_example, file.name = "new_dataset.raw", cross="f2 intercross")

## ---- chi_square--------------------------------------------------------------
f2_test <- test_segregation(bins_example)

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
twopts_f2 <- rf_2pts(input.obj = bins_example)

## ---- Suggest_a_LOD-----------------------------------------------------------
(LOD_sug <- suggest_lod(bins_example))

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
LG1_f2 <- make_seq(LGs_f2, 1)

## -----------------------------------------------------------------------------
LG1_f2

## ---- class_lg----------------------------------------------------------------
class(LG1_f2)

## ---- results="hide"----------------------------------------------------------
LG1_rcd_f2 <- rcd(LG1_f2)
LG1_rec_f2 <- record(LG1_f2)
LG1_ug_f2 <- ug(LG1_f2)

## ---- eval=FALSE--------------------------------------------------------------
#  LG1_ser_f2 <- seriation(LG1_f2) # Will return an error (can not be used in this case)

## -----------------------------------------------------------------------------
LG1_mds_f2 <- mds_onemap(input.seq = LG1_f2, rm_unlinked = TRUE)

## -----------------------------------------------------------------------------
rf_graph_table(LG1_rcd_f2)
rf_graph_table(LG1_rec_f2)
rf_graph_table(LG1_ug_f2)
rf_graph_table(LG1_mds_f2)

## ---- eval=FALSE--------------------------------------------------------------
#  rf_graph_table(LG1_ug_f2, inter = TRUE, html.file = "test.html")

## -----------------------------------------------------------------------------
# New LG1 will be this separated group
pos11 <- which(LG1_ug_f2$seq.num == 11) # Find position of marker 11
mksLG1 <- LG1_ug_f2$seq.num[1:pos11] # From marker 89 to 11

# LG2
pos23 <- which(LG1_ug_f2$seq.num == 23)
pos55 <- which(LG1_ug_f2$seq.num == 55)
mksLG2 <- LG1_ug_f2$seq.num[pos23:pos55]

# LG3
pos39 <- which(LG1_ug_f2$seq.num == 39)
mksLG3 <- LG1_ug_f2$seq.num[pos39:length(LG1_ug_f2$seq.num)] # use the position to find the 39 marker and take all the markers from there to the end of sequence

# Ordering again LG1
LG1 <- make_seq(twopts_f2, mksLG1)
LG1_ug2_f2 <- ug(LG1)

rf_graph_table(LG1_ug2_f2) # Now it is better

# Ordering LG2
LG2 <- make_seq(twopts_f2, mksLG2)
LG2_ug_f2 <- ug(LG2)

rf_graph_table(LG2_ug_f2)

# Ordering LG3
LG3 <- make_seq(twopts_f2, mksLG3)
LG3_ug_f2 <- ug(LG3)

rf_graph_table(LG3_ug_f2)

## ---- order_seq---------------------------------------------------------------
LG1_f2_ord <- order_seq(input.seq = LG1_ug2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3)
                        

## ---- show_order_seq----------------------------------------------------------
LG1_f2_ord

## ---- safe, results="hide"----------------------------------------------------
LG1_f2_safe <- make_seq(LG1_f2_ord, "safe")

## ---- force-------------------------------------------------------------------
(LG1_f2_all <- make_seq(LG1_f2_ord, "force"))

## ---- touchdown---------------------------------------------------------------
LG1_f2_ord <- order_seq(input.seq = LG1_ug2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ---- lg2_final---------------------------------------------------------------
(LG1_f2_final <- make_seq(LG1_f2_ord, "force"))

rf_graph_table(LG1_f2_final)

## ---- ripple_lg2_final, results="hide"----------------------------------------
ripple_seq(LG1_f2_final, ws = 5, LOD = 3)

## ---- results='hide'----------------------------------------------------------
LG2_f2_ord <- order_seq(input.seq = LG2_ug_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE, rm_unlinked = TRUE)

## -----------------------------------------------------------------------------
(LG2_f2_final <- make_seq(LG2_f2_ord, "force"))

rf_graph_table(LG2_f2_final)

## -----------------------------------------------------------------------------
LG2_edit <- drop_marker(LG2_f2_final, 23) # removing marker 23

## ---- eval=FALSE--------------------------------------------------------------
#  LG2_edit_map <- map(LG2_edit)

## -----------------------------------------------------------------------------
library(stringr)
LG2_edit_map <- onemap::map(LG2_edit) 

## -----------------------------------------------------------------------------
rf_graph_table(LG2_edit_map)

## -----------------------------------------------------------------------------
(LG2_temp <- try_seq(input.seq = LG2_edit_map, mrk = 23))

## -----------------------------------------------------------------------------
LG2_f2_final <- make_seq(LG2_temp, 1)

rf_graph_table(LG2_f2_final)

## ---- results="hide"----------------------------------------------------------
ripple_seq(LG2_f2_final, ws = 5)

## -----------------------------------------------------------------------------
LG2_f2_final

## ---- order_LG3, results='hide'-----------------------------------------------
LG3_f2_ord <- order_seq(input.seq = LG3_ug_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)

## ---- LG3_force---------------------------------------------------------------
(LG3_f2_final <- make_seq(LG3_f2_ord, "force"))

## -----------------------------------------------------------------------------
rf_graph_table(LG3_f2_final)

## -----------------------------------------------------------------------------
LG3_edit <- drop_marker(LG3_f2_final, c(34,39,50,56, 20,24, 64))
LG3_edit_map <- order_seq(LG3_edit) # We remove several markers maybe it's better to order again
LG3_edit_map <- make_seq(LG3_edit_map, "force")

rf_graph_table(LG3_edit_map)

## -----------------------------------------------------------------------------
LG3_edit <- try_seq(LG3_edit_map, 34)
LG3_edit_temp <- make_seq(LG3_edit, 1) # Not included
LG3_edit <- try_seq(LG3_edit_map, 39)
LG3_edit_temp <- make_seq(LG3_edit, 3) 
LG3_edit_map <- LG3_edit_temp # include
LG3_edit <- try_seq(LG3_edit_map, 50)
LG3_edit_temp <- make_seq(LG3_edit, 22) # Not included
LG3_edit <- try_seq(LG3_edit_map, 56)
LG3_edit_temp <- make_seq(LG3_edit, 22) # Not included 
LG3_edit <- try_seq(LG3_edit_map, 20)
LG3_edit_temp <- make_seq(LG3_edit, 22) 
LG3_edit_map <- LG3_edit_temp # include
LG3_edit <- try_seq(LG3_edit_map, 24)
LG3_edit_temp <- make_seq(LG3_edit, 22) 
LG3_edit_map <- LG3_edit_temp # include
LG3_edit <- try_seq(LG3_edit_map, 64)
LG3_edit_temp <- make_seq(LG3_edit, 23) # Not included

LG3_f2_final <- LG3_edit_map

## ---- ripple_LG3, results="hide"----------------------------------------------
ripple_seq(LG3_f2_final, ws = 5)

## -----------------------------------------------------------------------------
idx <- which(LG3_f2_final$seq.num == 59) 
new_seq <- LG3_f2_final$seq.num
new_seq[idx:(idx+4)] <- c(59, 49, 76, 80, 28)
LG3_edit_seq <- make_seq(twopts_f2, new_seq)

## -----------------------------------------------------------------------------
LG3_edit_map <- onemap::map(LG3_edit_seq)

## ---- LG3_final---------------------------------------------------------------
LG3_f2_final <- LG3_edit_map

rf_graph_table(LG3_f2_final)

## -----------------------------------------------------------------------------
CHR1 <- make_seq(twopts_f2, "1")
CHR1
CHR2 <- make_seq(twopts_f2, "2")
CHR3 <- make_seq(twopts_f2, "3")

## -----------------------------------------------------------------------------
CHR_mks <- group_seq(input.2pts = twopts_f2, seqs = "CHROM", unlink.mks = mark_all_f2,
                      repeated = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  CHR_mks <- group_seq(input.2pts = twopts_f2, seqs = list(CHR1=CHR1, CHR2=CHR2, CHR3=CHR3),
#                        unlink.mks = mark_all_f2, repeated = FALSE)

## -----------------------------------------------------------------------------
CHR_mks

## -----------------------------------------------------------------------------
CHR_mks$repeated

## -----------------------------------------------------------------------------
CHR_mks$sequences$CHR1
# or
CHR_mks$sequences[[1]]

## -----------------------------------------------------------------------------
LG3seq_f2 <- make_seq(twopts_f2, c(47, 38, 59, 16, 62, 21, 20, 48, 22))


## ---- eval=FALSE--------------------------------------------------------------
#  LG3seq_f2_map <- map(LG3seq_f2)

## -----------------------------------------------------------------------------
library(stringr)
(LG3seq_f2_map <- onemap::map(LG3seq_f2))

## ---- markers_names_and_numbers-----------------------------------------------
marker_type(LG3seq_f2_map)

## ---- add_marker--------------------------------------------------------------
(LG3seq_f2_map <- add_marker(LG3seq_f2_map, c(18, 56, 50)))

## ---- drop_marker-------------------------------------------------------------
(LG3seq_f2_map <- drop_marker(LG3seq_f2_map, c(59, 21)))

## -----------------------------------------------------------------------------
maps_list <- list(LG1_f2_final, LG2_f2_final, LG3_f2_final)

## -----------------------------------------------------------------------------
draw_map(maps_list, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## -----------------------------------------------------------------------------
draw_map(LG1_f2_final, names = TRUE, grid = TRUE, cex.mrk = 0.7)

## ---- eval=FALSE, results='hide', eval=FALSE----------------------------------
#  draw_map2(LG1_f2_final, LG2_f2_final, LG3_f2_final, main = "Only linkage information",
#            group.names = c("LG1", "LG2", "LG3"))

## ---- echo=FALSE, results='hide', echo=FALSE----------------------------------
draw_map2(LG1_f2_final, LG2_f2_final, LG3_f2_final, main = "Only linkage information", 
          group.names = c("LG1", "LG2", "LG3"), output = "map.png")

## ---- results='hide', eval=FALSE----------------------------------------------
#  draw_map2(LG1_f2_final, col.group = "#58A4B0", col.mark = "#335C81", output = "map_LG1.pdf")

## ---- results='hide', echo=FALSE----------------------------------------------
draw_map2(LG1_f2_final, col.group = "#58A4B0", col.mark = "#335C81", output = "map_LG1.png")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  data(parallel_results_out)
#  time_spent <-time_spent/(60*60)
#  colnames(time_spent) <- c("Without parallelization (h)", "With parallelization (h)" )
#  knitr::kable(time_spent)

## ---- eval=FALSE--------------------------------------------------------------
#  run_pedsim(chromosome = "Chr1", n.marker = 300, tot.size.cm = 100, centromere = 50,
#             n.ind = 200, mk.types = c("A.H.B", "C.A", "D.B"),
#             n.types = rep(100,3), pop = "F2", path.pedsim = "./",
#             name.mapfile = "mapfile.txt", name.founderfile="founderfile.gen",
#             name.chromfile="sim.chrom", name.parfile="sim.par",
#             name.out="simParall_f2")
#  
#  # Do the conversion
#  
#  pedsim2raw(cross="f2 intercross", genofile = "simParall_f2_genotypes.dat",
#             parent1 = "P1", parent2 = "P2", out.file = "simParall_f2.raw",
#             miss.perc = 25)
#  
#  # Import to R environment as onemap object
#  
#  simParallel <- read_onemap("simParall_f2.raw")
#  plot(simParallel)

## ---- eval=FALSE--------------------------------------------------------------
#  # Calculates two-points recombination fractions
#  twopts <- rf_2pts(simParallel)
#  
#  seq_all <- make_seq(twopts, "all")
#  
#  # There are no redundant markers
#  find_bins(simParallel)
#  
#  # There are no distorted markers
#  print(test_segregation(simParallel)) # Not shown

## ---- eval=FALSE--------------------------------------------------------------
#  batch_size <- pick_batch_sizes(input.seq = seq_all,
#                                 size = 80,
#                                 overlap = 30,
#                                 around = 10)
#  
#  batch_size

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  time_spent <- data.frame("without-parallelization"= rep(0,5), "with-parallelization" =rep(0,5))
#  rownames(time_spent) <- c("rcd", "record_map", "ug_map", "mds_onemap", "map")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # Without parallelization
#  time <- system.time(rcd_map <- rcd(input.seq = seq_all))
#  time_spent$without.parallelization[1] <- time[3]
#  
#  # With parallelization
#  time <- system.time(rcd_map_par <- rcd(input.seq = seq_all,
#                                         phase_cores = 4,
#                                         size = batch_size,
#                                         overlap = 30))
#  
#  time_spent$with.parallelization[1] <- time[3]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Without parallelization
#  rcd_map <- rcd(input.seq = seq_all)
#  
#  # With parallelization
#  rcd_map_par <- rcd(input.seq = seq_all,
#                     phase_cores = 4,
#                     size = batch_size,
#                     overlap = 30)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a <- rf_graph_table(rcd_map, mrk.axis = "none")
#  b <- rf_graph_table(rcd_map_par, mrk.axis = "none")
#  
#  p <- ggarrange(a,b , common.legend = TRUE,
#            labels = c("rcd", "rcd + parallel"),
#            vjust = 0.2,
#            hjust= -1.4,
#            font.label = list(size = 10),
#            ncol=2, nrow=1)
#  
#  ggsave(p, filename = "rcd.jpg")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # Without parallelization
#  time <- system.time(record_map <- record(input.seq = seq_all))
#  time_spent$without.parallelization[2] <- time[3]
#  
#  # With parallelization
#  time <- system.time(record_map_par <- record(input.seq = seq_all,
#                                               phase_cores = 4,
#                                               size = batch_size,
#                                               overlap = 30))
#  time_spent$with.parallelization[2] <- time[3]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Without parallelization
#  record_map <- record(input.seq = seq_all)
#  
#  # With parallelization
#  record_map_par <- record(input.seq = seq_all,
#                           phase_cores = 4,
#                           size = batch_size,
#                           overlap = 30)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a <- rf_graph_table(record_map, mrk.axis = "none")
#  b <- rf_graph_table(record_map_par, mrk.axis = "none")
#  
#  p <- ggarrange(a,b , common.legend = TRUE,
#            labels = c("record", "record + parallel"),
#            vjust = 0.2,
#            hjust= -0.8,
#            font.label = list(size = 10),
#            ncol=2, nrow=1)
#  
#  ggsave(p, filename = "record.jpg")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # Without parallelization
#  time <- system.time(ug_map <- ug(input.seq = seq_all))
#  time_spent$without.parallelization[3] <- time[3]
#  
#  # With parallelization
#  time <- system.time(ug_map_par <- ug(input.seq = seq_all,
#                                       phase_cores = 4,
#                                       size = batch_size,
#                                       overlap = 30))
#  
#  time_spent$with.parallelization[3] <- time[3]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Without parallelization
#  ug_map <- ug(input.seq = seq_all)
#  
#  # With parallelization
#  ug_map_par <- ug(input.seq = seq_all,
#                   phase_cores = 4,
#                   size = batch_size,
#                   overlap = 30)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a <- rf_graph_table(ug_map, mrk.axis = "none")
#  b <- rf_graph_table(ug_map_par, mrk.axis = "none")
#  
#  p <- ggarrange(a,b , common.legend = TRUE,
#            labels = c("ug", "ug + parallel"),
#            vjust = 0.2,
#            hjust= -1.6,
#            font.label = list(size = 10),
#            ncol=2, nrow=1)
#  
#  ggsave(p, filename = "ug.jpg")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # Without parallelization ok
#  time <- system.time(map_mds <- mds_onemap(input.seq = seq_all))
#  time_spent$without.parallelization[4] <- time[3]
#  
#  # With parallelization
#  time <- system.time(map_mds_par <- mds_onemap(input.seq = seq_all,
#                                                phase_cores = 4,
#                                                size = batch_size,
#                                                overlap = 30))
#  
#  time_spent$with.parallelization[4] <- time[3]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Without parallelization ok
#  map_mds <- mds_onemap(input.seq = seq_all)
#  
#  # With parallelization
#  map_mds_par <- mds_onemap(input.seq = seq_all,
#                            phase_cores = 4,
#                            size = batch_size,
#                            overlap = 30)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a <- rf_graph_table(map_mds, mrk.axis = "none")
#  b <- rf_graph_table(map_mds_par, mrk.axis = "none")
#  
#  p <- ggarrange(a,b , common.legend = TRUE,
#            labels = c("mds", "mds + parallel"),
#            vjust = 0.2,
#            hjust= -1,
#            font.label = list(size = 10),
#            ncol=2, nrow=1)
#  
#  ggsave(p, filename = "mds.jpg")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  batch_map <- map_overlapping_batches(input.seq = seq_all,
#                                       size = batch_size,
#                                       phase_cores = 4,
#                                       overlap = 30,
#                                       rm_unlinked = TRUE)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # Without parallelization
#  time <- system.time(batch_map <- map_avoid_unlinked(input.seq = seq_all))
#  time_spent$without.parallelization[5] <- time[3]
#  
#  # With parallelization
#  time <- system.time(batch_map_par <- map_avoid_unlinked(input.seq = seq_all,
#                                                          size = batch_size,
#                                                          phase_cores = 4,
#                                                          overlap = 30))
#  
#  time_spent$with.parallelization[5] <- time[3]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # Without parallelization
#  batch_map <- map_avoid_unlinked(input.seq = seq_all)
#  
#  # With parallelization
#  batch_map_par <- map_avoid_unlinked(input.seq = seq_all,
#                                      size = batch_size,
#                                      phase_cores = 4,
#                                      overlap = 30)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  a <- rf_graph_table(batch_map, mrk.axis = "none")
#  b <- rf_graph_table(batch_map_par, mrk.axis = "none")
#  
#  p <- ggarrange(a,b , common.legend = TRUE,
#            labels = c("map", "map + parallel"),
#            vjust = 0.2,
#            hjust= -1,
#            font.label = list(size = 10),
#            ncol=2, nrow=1)
#  
#  ggsave(p, filename = "map.jpg")

## ---- eval=FALSE--------------------------------------------------------------
#  time <- system.time(mds_ripple <- ripple_ord(input.seq = map_mds_par,
#                                               ws=5,
#                                               tol=10E-4,
#                                               start = 1,
#                                               method = "one",
#                                               n = NULL,
#                                               pref = NULL,
#                                               no_reverse = TRUE,
#                                               phase_cores = 2,
#                                               ripple_cores = 6,
#                                               verbose = c("time")))
#  

## -----------------------------------------------------------------------------
(progeny_haplot <- progeny_haplotypes(LG2_f2_final, most_likely = TRUE, ind = 2, group_names = "LG2_final"))

## -----------------------------------------------------------------------------
plot(progeny_haplot, position = "stack")

plot(progeny_haplot, position = "split")

## -----------------------------------------------------------------------------
sessionInfo()

