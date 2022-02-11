context("test ordering and HMM algorithms")

test_that("ordering and HMM test", {
  ordering_func <- function(example_data, right.order, right.size, tol.order, tol.size){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    all_mark <- make_seq(twopt,"all")
    LG <- make_seq(twopt, 1:27)
    if(!is(onemap_mis, "outcross")){
      LG.ser <- seriation(LG)
      eval(bquote(expect_equal(LG.ser$seq.num, .(right.order), tolerance = tol.order)))
      size <- cumsum(kosambi(LG.ser$seq.rf))
      eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    }
    LG.rcd <- rcd(LG)
    eval(bquote(expect_equal(LG.rcd$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.rcd$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.rec <- record(LG)
    eval(bquote(expect_equal(LG.rec$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.rec$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.ug <- ug(LG)
    eval(bquote(expect_equal(LG.ug$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.ug$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.mds <- mds_onemap(LG)
    eval(bquote(expect_equal(LG.mds$seq.num, .(right.order), tolerance = tol.order))) # mds makes local rearrangements
    size <- cumsum(kosambi(LG.mds$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size),tolerance = tol.size)))
    LG.order <- order_seq(LG)
    LG.order <- make_seq(LG.order, "force")
    eval(bquote(expect_equal(LG.order$seq.num, .(right.order), tolerance = tol.order))) # mds makes local rearrangements
    size <- cumsum(kosambi(LG.mds$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size),tolerance = tol.size)))
  }
  
  ordering_func("simu_example_out", 1:27, 100, 3, 10)
  ordering_func("simu_example_f2", right.order = 1:27, right.size = 100, tol.order = 3, 10)
  ordering_func("simu_example_bc",1:27, 100, 3,10)
  
})

test_that("ordering and HMM parallel", {
  ordering_func <- function(example_data, right.order, right.size, tol.order, tol.size, n.mar){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    LG <- make_seq(twopt, 1:27)
    batch_size <- pick_batch_sizes(LG,
                                   size = 5,
                                   overlap = 3,
                                   around = 1)

    if(!is(onemap_mis, "outcross")){
      LG.ser <- seriation(LG, size = batch_size, phase_cores = 1, overlap = 3)
      eval(bquote(expect_equal(LG.ser$seq.num, .(right.order), tolerance = tol.order)))
      size <- cumsum(kosambi(LG.ser$seq.rf))
      eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    }
    LG.rcd <- rcd(LG, size = batch_size, phase_cores = 1, overlap = 3)
    eval(bquote(expect_equal(LG.rcd$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.rcd$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.rec <- record(LG, size = batch_size, phase_cores = 1, overlap = 3)
    eval(bquote(expect_equal(LG.rec$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.rec$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.ug <- ug(LG, size = batch_size, phase_cores = 1, overlap = 3)
    eval(bquote(expect_equal(LG.ug$seq.num, .(right.order), tolerance = tol.order)))
    size <- cumsum(kosambi(LG.ug$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))
    LG.mds <- mds_onemap(LG, size = batch_size, phase_cores = 1, overlap = 3)
    eval(bquote(expect_equal(LG.mds$seq.num, .(right.order), tolerance = tol.order))) # mds makes local rearrangements
    size <- cumsum(kosambi(LG.mds$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size),tolerance = tol.size)))

    # Testing with only one core
    LG.map <- map(make_seq(LG.mds$twopt, LG.mds$seq.num))
    size <- cumsum(kosambi(LG.map$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))

    LG.map.avoid <- map_avoid_unlinked(make_seq(LG.mds$twopt, LG.mds$seq.num))
    size <- cumsum(kosambi(LG.map.avoid$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))

    LG.map.save <- onemap:::map_save_ram(input.seq = make_seq(LG.mds$twopt, LG.mds$seq.num))
    size <- cumsum(kosambi(LG.map.save$seq.rf))
    eval(bquote(expect_equal(size[length(size)], .(right.size), tolerance = tol.size)))

    # Testing try_seq_by_seq
    new_seq <- try_seq_by_seq(LG.map.save, make_seq(LG.map.save$twopt, 28:31)$seq.num)
    eval(bquote(expect_equal(length(new_seq$seq.num), .(n.mar))))
  }

  ordering_func("simu_example_out", 1:27, 100, 3, 10,27)
  ordering_func("simu_example_f2", 1:27, 100, 3, 10,27)
  ordering_func("simu_example_bc",1:27, 100, 3,10,27)
})
