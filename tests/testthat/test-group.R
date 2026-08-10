context("testing group functions")

test_that("group function",{
  check_group <- function(example_data, table.groups){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.25)))
    twopts <- rf_2pts(onemap_mis)
    lod_sug <- suggest_lod(onemap_mis)
    seq1 <- make_seq(twopts, "all")
    lgs <- group(seq1, lod_sug)
    expect_equal(as.vector(table(lgs$groups)), table.groups)
  }
  
  check_group("onemap_example_out", c(15,8,5,2))
  check_group(example_data = "onemap_example_f2", c(27,16,23))
  check_group("onemap_example_bc", c(17,28,22))
  check_group("onemap_example_riself", c(2,14,33,13))
})

test_that("group_seq function", {
  check_group_seq <- function(example_data1,example_data2, n.unlinked, n.repeated){
    eval(bquote(data(.(example_data1))))
    eval(bquote(data(.(example_data2))))
    comb <-  eval(bquote(combine_onemap(get(.(example_data1)), get(.(example_data2)))))
    onemap_mis <- filter_missing(comb, 0.25)
    twopts <- rf_2pts(onemap_mis)
    lgs <- group_seq(input.2pts = twopts, seqs = "CHROM", unlink.mks = "all", repeated = T)
    expect_equal(lgs$n.unlinked, n.unlinked)
    expect_equal(lgs$n.repeated, n.repeated)
  }
  
  check_group_seq("vcf_example_out", "onemap_example_out", 0, 0)
  check_group_seq("vcf_example_f2", "onemap_example_f2", 0, 0)
  check_group_seq("vcf_example_riself", "onemap_example_riself",3, 0)
  check_group_seq("vcf_example_bc", "onemap_example_bc",25, 22)
})

test_that("group_upgma", {
  check_group_upgma <- function(example_data, expected.groups, table.groups){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.25)))
    twopts <- rf_2pts(onemap_mis)
    seq1 <- make_seq(twopts, "all")
    lgs <- group_upgma(seq1, expected.groups = expected.groups, comp.mat = T, inter = F)
    expect_equal(as.vector(table(lgs$groups)), table.groups)
  }

  check_group_upgma("vcf_example_out", 2, c(11,13))
  check_group_upgma("vcf_example_f2", 3,  c(9,13,3))
  check_group_upgma(example_data = "vcf_example_riself", 2, c(13, 11))
  check_group_upgma("vcf_example_bc", 2, c(18,7))
})
