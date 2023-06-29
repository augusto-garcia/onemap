context("Export functions")

test_that("Combine and split datasets", {
  
  data("vcf_example_out")
  
  twopts <- rf_2pts(vcf_example_out)
  
  seq_all <- make_seq(twopts, "all")
  lgs <- group(seq_all)
  lg1 <- map(make_seq(lgs,1))
  lg2 <- map(make_seq(lgs,2))
  lg3 <- map(make_seq(lgs,3))
  lg4 <- map(make_seq(lgs,4))
  
  seqs.list <- list(lg1, lg2, lg3, lg4)
  
  viewmap.obj <- export_viewpoly(seqs.list)

  expect_equal(names(viewmap.obj), c("d.p1", "d.p2", "ph.p1", "ph.p2", "maps", "software"))
})
