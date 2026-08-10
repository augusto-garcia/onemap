context("testing two-points estimations")

test_that("two-points tests",{
  expect_twopts_phases <- function(example_data, phase1 , phase2, phase3, phase4, sum_all, nmks1, nmks2){
    eval(bquote(data(.(example_data))))
    twopts <- eval(bquote(rf_2pts(get(.(example_data)), verbose = TRUE)))
    expect_equal(check_twopts(twopts),0)
    eval(bquote(expect_equal(as.vector(twopts$analysis[[1]][1:4,1]), .(phase1), tolerance = 0.01)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[2]][1:4,1]), .(phase2), tolerance = 0.01)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[3]][1:4,1]), .(phase3), tolerance = 0.01)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[4]][1:4,1]), .(phase4), tolerance = 0.01)))
    eval(bquote(expect_equal(sum(twopts$analysis[[1]]), .(sum_all), tolerance = 0.01)))
    
    ld_sug <- eval(bquote(suggest_lod(get(.(example_data)))))
    lgs <- group(make_seq(twopts, "all"), LOD = ld_sug)
    lg1 <- make_seq(lgs, 1)
    if(lgs$n.groups > 1){
      lg2 <- make_seq(lgs, 2)
      seq1 <- make_seq(twopts, c(lg1$seq.num, lg2$seq.num[1:5]))
    } else {
      seq1 <- make_seq(twopts, c(lg1$seq.num))
    }
    seq1 <- mds_onemap(seq1, hmm = F)
    # Filters
    filt_seq1 <- rf_snp_filter_onemap(input.seq = seq1, probs = c(0.25,1))
    eval(bquote(expect_equal(length(filt_seq1$seq.num), nmks1, tolerance = 2)))
    filt_seq1 <- filter_2pts_gaps(seq1, max.gap = 20)
    eval(bquote(expect_equal(length(filt_seq1$seq.num), nmks2, tolerance = 2)))
  }
  
  expect_twopts_phases(example_data = "onemap_example_out", 
                       c(0, 0.1818151, 0.2954514, 0.573),
                       c(0, 0.8181849, 0.2954514, 0.5089259),
                       c(0, 0.1818151, 0.7045486, 0.4910741),
                       c(0, 0.8181849, 0.7045486, 0.4253372),
                       891, 15,10)
  
  expect_twopts_phases(example_data = "onemap_example_f2", 
                       c(0, 0.447, 0.159, 0.772),
                       c(0, 0.5, 0.5, 0.5),
                       c(0, 0.5,0.5,0.5),
                       c(0, 0.553, 0.841, 0.228),
                       7853, 26, 34)
  
  expect_twopts <- function(example_data, values, sum_all, nmks1, nmks2){
    eval(bquote(data(.(example_data))))
    onemap.mis <- eval(bquote(filter_missing(get(.(example_data)), 0.25)))
    twopts <- eval(bquote(rf_2pts(onemap.mis)))
    expect_equal(check_twopts(twopts),0)
    eval(bquote(expect_equal(as.vector(twopts$analysis[1:4,1]), .(values), tolerance = 0.01)))
    eval(bquote(expect_equal(sum(twopts$analysis), .(sum_all), tolerance = 0.01)))
    
    lgs <- group(make_seq(twopts, "all"))
    lg1 <- make_seq(lgs, 1)
    lg2 <- make_seq(lgs, 2)
    seq1 <- mds_onemap(make_seq(twopts, c(lg1$seq.num, lg2$seq.num[1:5])), hmm = F)
    # Filters
    filt_seq1 <- rf_snp_filter_onemap(seq1, probs = c(0.25,1))
    eval(bquote(expect_equal(length(filt_seq1$seq.num), nmks1)))
    filt_seq1 <- filter_2pts_gaps(seq1, max.gap = 20)
    eval(bquote(expect_equal(length(filt_seq1$seq.num), nmks2)))
  }
  
  expect_twopts(example_data = "onemap_example_bc",
                c(0.00000000, 0.0364, 0.50000000, 0.495000), 6715, 16, 16)
  
  expect_twopts(example_data = "onemap_example_riself",
                c(0.00000000, 0.495, 0.5000000, 0.5000000), 5546, 14, 14)
})
