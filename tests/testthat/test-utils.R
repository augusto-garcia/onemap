context("Utils functions")


test_that("Combine and split datasets", {
  check_combine <- function(data1, data2, n.ind, n.mks){
    
    eval(bquote(data(.(data1))))
    eval(bquote(data(.(data2))))
    
    split.mks <- eval(bquote(colnames(get(.(data1))$geno)))
    split.mks2 <- eval(bquote(colnames(get(.(data2))$geno)))
    
    out <- eval(bquote(combine_onemap(get(.(data1)), get(.(data2)))))
  
    out_filt <- filter_missing(out, threshold = 0.5)
    twopts <- rf_2pts(out_filt)
    
    expect_equal(check_data(out), 0)
    eval(bquote(expect_equal(out$n.ind, .(n.ind))))
    eval(bquote(expect_equal(out$n.mar, .(n.mks))))
    
    out_sep <- split_onemap(out, split.mks)
    expect_equal(check_data(out_sep), 0)
    
    eval(bquote(expect_equal(out_sep$n.ind, get(.(data1))$n.ind)))
    eval(bquote(expect_equal(out_sep$n.mar, get(.(data1))$n.mar)))
    
    out_2pts <- split_2pts(twopts, split.mks2)
    expect_equal(check_twopts(out_2pts), 0)
    
    eval(bquote(expect_equal(out_2pts$n.mar, get(.(data2))$n.mar)))
    
  }
  
  check_combine("onemap_example_out", "vcf_example_out", 100, 54)
  check_combine("onemap_example_f2", "vcf_example_f2", 200, 91)
  check_combine("onemap_example_bc", "vcf_example_bc", 150, 92)
  check_combine("onemap_example_riself", "vcf_example_riself", 100, 93)
})


test_that("Edit onemap object", {
  data("onemap_example_out")
  new_data <- remove_inds(onemap_example_out, rm.ind = c("IND1","IND5"))
  expect_equal(check_data(new_data),0)
  expect_equal(onemap_example_out$n.ind - new_data$n.ind, 2)

  data("vcf_example_out")
  new_data <- sort_by_pos(vcf_example_out)
  expect_equal(check_data(new_data),0)
})

