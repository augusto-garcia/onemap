context("test ordering algorithms")

test_that("ordering test", {
  ordering_func <- function(example_data, which.group, ord.ser, ord.rcd, ord.rec, ord.ug, ord.mds, ord.order){
    eval(bquote(data(.(example_data))))
    twopt <- eval(bquote(rf_2pts(get(.(example_data)))))
    all_mark <- make_seq(twopt,"all")
    lod_sug <- suggest_lod(all_mark)
    groups <- group(all_mark, LOD = lod_sug)
    LG <- make_seq(groups,which.group)
    LG <- make_seq(LG$twopt, LG$seq.num[1:5])
    set.seed(2020)
    LG.ser <- seriation(LG)
    eval(bquote(expect_equal(LG.ser$seq.num[1:5], .(ord.ser))))
    LG.rcd <- rcd(LG)
    eval(bquote(expect_equal(LG.rcd$seq.num[1:5], .(ord.rcd))))
    LG.rec <- record(LG)
    eval(bquote(expect_equal(LG.rec$seq.num[1:5], .(ord.rec))))
    LG.ug <- ug(LG)
    eval(bquote(expect_equal(LG.ug$seq.num[1:5], .(ord.ug))))
    LG.mds <- mds_onemap(LG)
    eval(bquote(expect_equal(LG.mds$seq.num[1:5], .(ord.mds))))
    set.seed(2021)
    LG.order_seq <- order_seq(LG, n.init = 3, twopt.alg = "rcd")
    LG.order <- make_seq(LG.order_seq, "force")
    eval(bquote(expect_equal(LG.order$seq.num[1:5], .(ord.order))))
  }
  
  ordering_func("onemap_example_out", 3,
                c(7,22,13,8,18), 
                c(7,22,13,8,18),
                c(7,13,18,8,22),
                c(7,22,13,8,18), 
                c(7,22,13,8,18), 
                c(22,18,8,13,7))
  
  ordering_func("onemap_example_f2", 1,
                c(4,1,3,5,2), 
                c(2,5,3,1,4),
                c(4,1,3,5,2),
                c(4,1,3,2,5), 
                c(1,3,4,2,5), 
                c(3,1,4,5,2))
  
  ordering_func("onemap_example_bc", 1,
                c(20,2,1,14,28), 
                c(20,2,1,14,28),
                c(20,1,2,14,28),
                c(20,1,2,14,28), 
                c(28,14,2,1,20), 
                c(20,2,1,14,28))
  
  ordering_func("onemap_example_riself", 1,
                c(8,1,11,16,17), 
                c(8,1,16,11,17),
                c(8,1,11,16,17),
                c(8,1,17,16,11), 
                c(8,1,17,11,16), 
                c(8,1,17,11,16))
})

test_that("ordering test parallel", {
  ordering_func <- function(example_data, ord.rcd, ord.rec, ord.ug, ord.mds){
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
    LG.rec <- record(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.rec$seq.num[1:5], .(ord.rec))))
    LG.ug <- ug(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.ug$seq.num[1:5], .(ord.ug))))
    LG.mds <- mds_onemap(LG, size = batch_size, phase_cores = 4, overlap = 3)
    eval(bquote(expect_equal(LG.mds$seq.num[1:5], .(ord.mds))))
  }
  
  ordering_func("onemap_example_out", 
                c(6,12,17,10,1), 
                c(22,8,18,13,7),
                c(25,15,5,11,6),
                c(10,2,23,14,17))
  
  ordering_func("onemap_example_f2",
                c(17,13,5,3,2), 
                c(13,17,21,16,19),
                c(7,9,24,10,6),
                c(18,2,22,1,14))
  
  ordering_func("onemap_example_bc",
                c(9,7,14,3,23),
                c(4,6,19,17,5), 
                c(9,7,14,3,12), 
                c(4,6,8,20,14))
  
  ordering_func("onemap_example_riself",
                c(7,1,19,10,15), 
                c(7,1,19,15,10),
                c(21,6,14,3,13),
                c(7,17,11,16,25))
})