context("test ordering and HMM algorithms")

test_that("ordering and HMM test", {
  ordering_func <- function(example_data, which.group, 
                            ord.ser, size.ser, 
                            ord.rcd, size.rcd, 
                            ord.rec, size.rec,
                            ord.ug, size.ug,
                            ord.mds, size.mds, 
                            ord.order, size.order){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    all_mark <- make_seq(twopt,"all")
    lod_sug <- suggest_lod(all_mark)
    groups <- group(all_mark, LOD = lod_sug)
    LG <- make_seq(groups,which.group)
    LG <- make_seq(LG$twopt, LG$seq.num[1:5])
    set.seed(2020)
    LG.ser <- seriation(LG)
    eval(bquote(expect_equal(LG.ser$seq.num[1:5], .(ord.ser))))
    size <- cumsum(kosambi(LG.ser$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.ser), tolerance = 0.001)))
    LG.rcd <- rcd(LG)
    eval(bquote(expect_equal(LG.rcd$seq.num[1:5], .(ord.rcd))))
    size <- cumsum(kosambi(LG.rcd$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.rcd), tolerance = 0.001)))
    LG.rec <- record(LG)
    eval(bquote(expect_equal(LG.rec$seq.num[1:5], .(ord.rec))))
    size <- cumsum(kosambi(LG.rec$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.rec), tolerance = 0.001)))
    LG.ug <- ug(LG)
    eval(bquote(expect_equal(LG.ug$seq.num[1:5], .(ord.ug))))
    size <- cumsum(kosambi(LG.ug$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.ug), tolerance = 0.001)))
    LG.mds <- mds_onemap(LG)
    eval(bquote(expect_equal(LG.mds$seq.num[1:5], .(ord.mds))))
    size <- cumsum(kosambi(LG.mds$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.mds),tolerance = 0.001)))
    set.seed(2021)
    LG.order_seq <- order_seq(LG, n.init = 3, twopt.alg = "rcd")
    LG.order <- make_seq(LG.order_seq, "force")
    eval(bquote(expect_equal(LG.order$seq.num[1:5], .(ord.order))))
    size <- cumsum(kosambi(LG.order$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.order), tolerance = 0.001)))
  }
  
  ordering_func("onemap_example_out", 3,
                c(7,22,13,8,18), 46.32121,
                c(7,22,13,8,18), 46.32121394,
                c(7,13,18,8,22), 94.87187, 
                c(7,22,13,8,18), 46.32121394,
                c(7,22,13,8,18), 46.32121394,
                c(22,18,8,13,7), 93.03766)
  
  ordering_func("onemap_example_f2", 1,
                c(2,1,3,5,4), 50.3,
                c(2,1,3,5,4), 50.3,
                c(2,1,3,5,4), 50.3,
                c(2,1,3,5,4), 50.3,
                c(2,1,3,5,4), 50.3,
                c(4,5,3,1,2), 50.3)
  
  ordering_func("onemap_example_bc", 1,
                c(15,1,2,18,22), 55.8,
                c(15,1,2,18,22), 55.8,
                c(15,1,2,18,22), 55.8,
                c(15,1,2,18,22), 55.8,
                c(15,1,2,18,22), 55.8,
                c(15,1,2,18,22), 55.8)
  
  ordering_func("onemap_example_riself", 1,
                c(7,1,10,15,16), 40.61501,
                c(7,1,15,10,16), 40.61501,
                c(7,1,10,15,16), 40.61501,
                c(7,1,16,15,10), 42.48526,
                c(7,1,16,10,15), 42.48526,
                c(7,1,16,10,15), 42.5)
})

test_that("ordering and HMM parallel", {
  ordering_func <- function(example_data, 
                            ord.rcd, size.rcd, 
                            ord.rec, size.rec,
                            ord.ug, size.ug,
                            ord.mds, size.mds){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    LG <- make_seq(twopt, 1:25)
    batch_size <- pick_batch_sizes(LG,
                                   size = 5,
                                   overlap = 3,
                                   around = 1)
    set.seed(2020)
    eval(bquote(expect_error(seriation(input.seq = LG, size = batch_size, phase_cores = 4, overlap = 3), "There are too many ties in the ordination process - please, consider using another ordering algorithm.")))
    LG.rcd <- rcd(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.rcd$seq.num[1:5], .(ord.rcd))))
    size <- cumsum(kosambi(LG.rcd$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.rcd), tolerance = 0.001)))
    LG.rec <- record(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.rec$seq.num[1:5], .(ord.rec))))
    size <- cumsum(kosambi(LG.rec$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.rec), tolerance = 0.001)))
    LG.ug <- ug(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.ug$seq.num[1:5], .(ord.ug))))
    size <- cumsum(kosambi(LG.ug$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.ug), tolerance = 0.001)))
    LG.mds <- mds_onemap(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.mds$seq.num[1:5], .(ord.mds))))
    size <- cumsum(kosambi(LG.mds$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(size.mds), tolerance = 0.001)))
  }
  
  ordering_func("onemap_example_out", 
                c(6,12,17,10,1), 622.95326,
                c(22,8,18,13,7), 516.82254,
                c(25,15,5,11,6), 487.86333,
                c(10,2,23,14,17), 927.20040)
  
  ordering_func("onemap_example_f2",
                c(17,13,5,3,2), 383.32966,
                c(13,17,21,16,19), 330.57134,
                c(7,9,24,10,6), 632.506725,
                c(18,2,22,1,14), 1331.10242)
  
  ordering_func("onemap_example_bc",
                c(9,7,14,3,23), 512.690549,
                c(4,6,19,17,5), 804.31139,
                c(9,7,14,3,12), 774.573013,
                c(4,6,8,20,14), 2701.98197)
  
  ordering_func("onemap_example_riself",
                c(7,1,19,10,15), 313.01,
                c(7,1,19,15,10), 286.21189,
                c(21,6,14,3,13), 979.791324,
                c(7,17,11,16,25), 1106.8291)
})