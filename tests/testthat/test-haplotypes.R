context("test plot haplotypes")

test_that("ordering and HMM parallel", {
  test_haplo <- function(example_data, which.group, sum.counts, parent_haplo = NULL, parent_haplo_ref_alt = NULL){
    onemap_mis <- eval(bquote(filter_missing(.(example_data), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    all_mark <- make_seq(twopt,"all")
    lod_sug <- suggest_lod(all_mark)
    groups <- group_upgma(all_mark, expected.groups = 2, inter = F)
    LG <-  eval(bquote(make_seq(groups, .(which.group))))
    map1 <- onemap::map(LG)
    dist <- cumsum(kosambi(map1$seq.rf))
    expect_equal(dist[length(dist)], 100, tolerance = 5) # The simulated distance of Chr01 is 100
    if(!inherits(onemap_mis, "backcross") & !inherits(onemap_mis, "ri")){
      haplo_default <- parents_haplotypes(map1)
      to_match <- unlist(haplo_default[3, 5:8])
      names(to_match) <- NULL
      eval(bquote(expect_equal(to_match, .(parent_haplo))))
      if(!is.null(onemap_mis$ref_alt_alleles)){
        haplo_ref_alt <- parents_haplotypes(map1, ref_alt_alleles = TRUE)
        to_match <- unlist(haplo_ref_alt[3, 5:8])
        names(to_match) <- NULL
        eval(bquote(expect_equal(to_match, .(parent_haplo_ref_alt))))
      }
    }
    counts <- progeny_haplotypes_counts(x = progeny_haplotypes(map1, most_likely = T, ind = "all"))
    eval(bquote(expect_equal(sum(counts$counts), .(sum.counts), tolerance = 0.01)))
  }
  
  data("simu_example_bc")
  test_haplo(example_data = simu_example_bc, 1, 126)
  data("simu_example_out")
  test_haplo(example_data = simu_example_out, 1, 347, parent_haplo = c("a", "a", "b", "a"))
  data("simu_example_f2")
  test_haplo(example_data = simu_example_f2, 1, 218, parent_haplo = c("a", "b", "a", "b"))
  
  example_data <- onemap_read_vcfR(vcf = system.file("extdata/simu_cod_out.vcf.gz", package = "onemap"),
                                   parent1 = "P2",
                                   parent2 = "P1",
                                   cross = "outcross", only_biallelic = FALSE)
  
  test_haplo(example_data, which.group = 1, sum.counts = 336, c("c", "c", "a", "b"), c("C", "C", "A", "G"))
  
  example_data <- onemap_read_vcfR(vcf = system.file("extdata/simu_cod_f2.vcf.gz", package = "onemap"),
                                   parent1 = "P1",
                                   parent2 = "P2",
                                   cross = "f2 intercross", only_biallelic = FALSE)
  
  test_haplo(example_data, which.group = 1, sum.counts = 306, c("a", "b", "a", "b"), c("T", "A", "T", "A"))
  
})
