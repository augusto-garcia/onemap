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
  
  genoprob_list <- lapply(seqs.list, export_mappoly_genoprob)

  expect_equal(sum(genoprob_list[[1]]$map), 66.24, tolerance = 1)
  
  #homoprob <- mappoly::calc_homologprob(genoprob_list)
  #mappoly::plot(homoprob, lg = 2)
  
  # ind.names <- dimnames(genoprob_list[[1]]$probs)[[3]]
  # fake.pheno <- matrix(sample(32:100, length(ind.names)*3, replace = TRUE), nrow=length(ind.names))
  # rownames(fake.pheno) <- ind.names
  # colnames(fake.pheno) <- paste0("pheno", 1:3)
  # 
  # library(qtlpoly)
  # data = read_data(ploidy = 2, geno.prob = genoprob_list, pheno = fake.pheno, step = 1) # fix
  # print(data, detailed = TRUE)
  # 
  # remim.mod = remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 1e-04,
  #                   d.sint = 1.5, n.clusters = 4)
  # print(remim.mod)
  # 
  # 
  # data("maps4x")
  # data("pheno4x")
  # genoprob4x = lapply(maps4x, mappoly::calc_genoprob)
  # data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1) 
  
})
