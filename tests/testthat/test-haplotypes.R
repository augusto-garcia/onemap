context("test plot haplotypes")

test_that("ordering and HMM parallel", {
  test_haplo <- function(example_data, which.group, sum.counts){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.15)))
    twopt <- rf_2pts(onemap_mis)
    all_mark <- make_seq(twopt,"all")
    lod_sug <- suggest_lod(all_mark)
    groups <- group_upgma(all_mark, expected.groups = 2, inter = F)
    LG <-  eval(bquote(make_seq(groups, .(which.group))))
    map1 <- onemap::map(LG)
    dist <- cumsum(kosambi(map1$seq.rf))
    expect_equal(dist[length(dist)], 100, tolerance = 5) # The simulated distance of Chr01 is 100
    counts <- progeny_haplotypes_counts(x = progeny_haplotypes(map1, most_likely = T, ind = "all"))
    eval(bquote(expect_equal(sum(counts$counts), .(sum.counts))))
  }
  
  test_haplo("simu_example_bc", 1, 126)
  test_haplo("simu_example_out", 1, 347)
  test_haplo("simu_example_f2", 1, 216)
})
