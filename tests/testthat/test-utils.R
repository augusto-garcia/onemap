context("Utils functions")


test_that("Combine and split datasets", {
  check_combine <- function(data1, data2, n.ind, n.mks, n_mk.end, rm_ind){
    
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
    
    seq1 <- make_seq(out_2pts, c(1:14))
    seq2 <- make_seq(out_2pts, c(14:23))
    list.sequences <- list(seq1, seq2)

    new.seqs <- keep_only_selected_mks(list.sequences)
    n_mk <- dim(new.seqs[[1]]$twopt$data.name$geno)[2]
    
    eval(bquote(expect_equal(n_mk, .(n_mk.end))))
    
    obj_up <- eval(bquote(remove_inds(out_filt, rm.ind = .(rm_ind))))
    
    obj_up <- eval(bquote(remove_inds(rm.ind = .(rm_ind), list.seqs = list.sequences)))
  }
  
  check_combine(data1 = "onemap_example_out", 
                data2 = "vcf_example_out", 
                n.ind = 100, 
                n.mks = 54,
                n_mk.end = 23,
                rm_ind = "IND2")
  
  check_combine(data1 = "onemap_example_f2", 
                data2 = "vcf_example_f2",
                n.ind =  200, 
                n.mks = 91, 
                n_mk.end = 23,
                rm_ind = "IND2")
  
  check_combine(data1 = "onemap_example_bc", 
                data2 = "vcf_example_bc", 
                n.ind = 150, 
                n.mks = 92, 
                n_mk.end = 23,
                rm_ind = "ID2")
  
  check_combine(data1 = "onemap_example_riself", 
                data2 = "vcf_example_riself", 
                n.ind = 100, 
                n.mks = 93, 
                n_mk.end = 23,
                rm_ind = "ID2")
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

