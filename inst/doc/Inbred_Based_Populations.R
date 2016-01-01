## ---- globalsetup, echo=FALSE, results='hide', cache=FALSE---------------
#opts_chunk$set(cache=TRUE, autodep=TRUE)

## ------------------------------------------------------------------------
library(onemap)

## ---- eval=FALSE---------------------------------------------------------
#  save.image("C:/.../yourfile.RData")

## ---- eval=FALSE---------------------------------------------------------
#  fake.f2.onemap <- read.mapmaker(dir="C:/workingdirectory",
#                                  file="your_data_file.raw")

## ---- load_data----------------------------------------------------------
data(fake.f2.onemap)

## ------------------------------------------------------------------------
fake.f2.onemap

## ---- class_of_object----------------------------------------------------
class(fake.f2.onemap)

## ---- plot_raw_data------------------------------------------------------
plot(fake.f2.onemap)

## ---- eval=FALSE---------------------------------------------------------
#  ?plot.f2.onemap

## ---- plot_by_type-------------------------------------------------------
plot_by_segreg_type(fake.f2.onemap)

## ---- chi_square---------------------------------------------------------
f2.test <- test_segregation(fake.f2.onemap)

## ---- class_of_chisquare_test--------------------------------------------
class(f2.test)

## ---- print_chi1---------------------------------------------------------
f2.test

## ---- print_chi2---------------------------------------------------------
print(f2.test)

## ---- Bonferroni---------------------------------------------------------
Bonferroni_alpha(f2.test)

## ---- plot_chisq---------------------------------------------------------
plot(f2.test)

## ---- select_nondistorted------------------------------------------------
(non.dist <- select_segreg(f2.test))

## ---- select_distorted---------------------------------------------------
(distorted <- select_segreg(f2.test, distorted=TRUE))

## ---- two_point_tests, results="hide"------------------------------------
twopts.f2 <- rf.2pts(fake.f2.onemap)

## ---- Suggest_a_LOD------------------------------------------------------
suggest_lod(fake.f2.onemap)

## ---- print_2Mks_two_points----------------------------------------------
print(twopts.f2, "M12", "M42")

## ---- class_of_twopoint--------------------------------------------------
class(twopts.f2)

## ---- print_all_two_points-----------------------------------------------
print(twopts.f2)

## ---- subset_all---------------------------------------------------------
mark.all.f2 <- make.seq(twopts.f2, "all")

## ---- class_subset_all---------------------------------------------------
class(mark.all.f2)

## ---- subset_3mks--------------------------------------------------------
mrk.subset<-make.seq(twopts.f2, c(1,3,7))

## ---- group1-------------------------------------------------------------
LGs.f2 <- group(mark.all.f2)
LGs.f2

## ---- group2-------------------------------------------------------------
(LGs.f2 <- group(mark.all.f2, LOD=4, max.rf=0.5))

## ---- class_group--------------------------------------------------------
class(LGs.f2)

## ---- haldane, eval=FALSE------------------------------------------------
#  set.map.fun(type="haldane")

## ---- kosambi, eval=FALSE------------------------------------------------
#  set.map.fun(type="kosambi")

## ------------------------------------------------------------------------
LG2.f2 <- make.seq(LGs.f2, 2)

## ------------------------------------------------------------------------
LG2.f2

## ---- class_lg-----------------------------------------------------------
class(LG2.f2)

## ---- results="hide"-----------------------------------------------------
LG2.ser.f2 <- seriation(LG2.f2)
LG2.rcd.f2 <- rcd(LG2.f2)
LG2.rec.f2 <- record(LG2.f2)
LG2.ug.f2 <- ug(LG2.f2)

## ---- order_seq, eval=FALSE----------------------------------------------
#  LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5,
#                          subset.search = "twopt",
#                          twopt.alg = "rcd", THRES = 3,
#                          draw.try = TRUE, wait = 1)
#  

## ---- order_seq2, echo=FALSE---------------------------------------------
LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        draw.try = FALSE, wait = 1)
                        

## ---- show_order_seq-----------------------------------------------------
LG2.f2.ord

## ---- safe, results="hide"-----------------------------------------------
LG2.f2.safe <- make.seq(LG2.f2.ord,"safe")

## ---- force--------------------------------------------------------------
(LG2.f2.all <- make.seq(LG2.f2.ord,"force"))

## ---- touchdown, eval=FALSE----------------------------------------------
#  LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5,
#                          subset.search = "twopt",
#                          twopt.alg = "rcd", THRES = 3,
#                          draw.try = TRUE, wait = 1,
#                          touchdown=TRUE)

## ---- touchdown2, echo=FALSE---------------------------------------------
LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        draw.try = FALSE, wait = 1,
                        touchdown=TRUE)

## ---- lg2_final----------------------------------------------------------
(LG2.f2.final<-make.seq(LG2.f2.ord, "force"))

## ---- ripple_lg2_final, results="hide"-----------------------------------
ripple_seq(LG2.f2.final, ws=5, LOD=3)

## ------------------------------------------------------------------------
LG2.f2.final

## ------------------------------------------------------------------------
LG1.f2 <- make.seq(LGs.f2, 1)

## ------------------------------------------------------------------------
LG1.f2.ord <- order.seq(input.seq=LG1.f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        draw.try = FALSE, wait = 1,
                        touchdown=TRUE)

## ------------------------------------------------------------------------
(LG1.f2.final <- make.seq(LG1.f2.ord,"force"))

## ---- results="hide"-----------------------------------------------------
ripple_seq(ws=5, LG1.f2.final)

## ------------------------------------------------------------------------
LG1.f2.final

## ------------------------------------------------------------------------
LG3.f2 <- make.seq(LGs.f2, 3)

## ---- order_LG3----------------------------------------------------------
LG3.f2.ord <- order.seq(input.seq=LG3.f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        draw.try = FALSE, wait = 1,
                        touchdown=TRUE)

## ---- LG3_force----------------------------------------------------------
(LG3.f2.final <- make.seq(LG3.f2.ord,"force"))

## ---- ripple_LG3, results="hide"-----------------------------------------
ripple_seq(ws=5, LG3.f2.final)

## ---- LG3_final----------------------------------------------------------
LG3.f2.final

## ---- subset_twopoint----------------------------------------------------
LG3seq.f2 <- make.seq(twopts.f2,c(47,38,59,16,62,21,20,48,22))
(LG3seq.f2.map <- map(LG3seq.f2))

## ---- markers_names_and_numbers------------------------------------------
marker.type(LG3seq.f2.map)

## ---- add_marker---------------------------------------------------------
(LG3seq.f2.map <- add.marker(LG3seq.f2.map, c(18,56,50)))

## ---- drop_marker--------------------------------------------------------
(LG3seq.f2.map <- drop.marker(LG3seq.f2.map, c(59,21)))

## ---- remove_M38, results="hide"-----------------------------------------
temp.seq<-drop.marker(LG3.f2.final, 38)

## ---- add_M38_end_LG-----------------------------------------------------
(temp.seq<-add.marker(temp.seq, 38))
(LG3.f2.wrong<-map(temp.seq))

## ---- rec_matrix---------------------------------------------------------
rf.graph.table(LG3.f2.wrong, inter=FALSE)

## ---- rec_matrix_inter, eval=FALSE---------------------------------------
#  rf.graph.table(LG3.f2.wrong)

## ---- try_M38, width=9, height=9-----------------------------------------
temp.seq <- drop.marker(LG3.f2.wrong,38)
temp.map <- map(temp.seq)
temp.try <- try_seq(temp.map, 38, draw.try=TRUE)

## ------------------------------------------------------------------------
(LG3.f2.final<-make.seq(temp.try, 4))

## ------------------------------------------------------------------------
maps.list<-list(LG1.f2.final, LG2.f2.final, LG3.f2.final)

## ------------------------------------------------------------------------
draw.map(maps.list, names= TRUE, grid=TRUE, cex.mrk=0.7)

## ------------------------------------------------------------------------
draw.map(LG1.f2.final, names= TRUE, grid=TRUE, cex.mrk=0.7)

## ---- write_map----------------------------------------------------------
write.map(maps.list, "fake.f2.onemap.map")

## ---- install_qtl, eval=FALSE--------------------------------------------
#  install.packages("qtl")

## ---- load_qtl-----------------------------------------------------------
library(qtl)

## ---- where_is_raw-------------------------------------------------------
raw.file<-paste(system.file("example",package="onemap"),
                "fake.f2.onemap.raw", sep="/")

## ---- read_data_in_qtl---------------------------------------------------
fake.f2.qtl <- read.cross("mm", file=raw.file, mapfile="fake.f2.onemap.map")

## ---- map_in_qtl---------------------------------------------------------
newmap <- est.map(fake.f2.qtl, tol=1e-6, map.function="kosambi")

## ---- compare_qtl_onemap-------------------------------------------------
plot.map(fake.f2.qtl, newmap)

## ---- IM_qtl-------------------------------------------------------------
fake.f2.qtl <- calc.genoprob(fake.f2.qtl, step=2)
out.em <- scanone(fake.f2.qtl, method="em")
out.hk <- scanone(fake.f2.qtl, method="hk")
plot(out.em, out.hk, col=c("blue","red"))

## ------------------------------------------------------------------------
write.cross(fake.f2.qtl, format="qtlcart", filestem="fake.f2.onemap")

