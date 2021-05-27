context("Filters function")

test_that("number of distorted markers",{
  check_dist <- function(example_data, table.h0){
    eval(bquote(data(.(example_data))))
    segre <- eval(bquote(test_segregation(get(.(example_data)))))
    segre_tab <- print(segre)
    eval(bquote(expect_equal(as.vector(table(segre_tab$H0)), .(table.h0))))
    expect_equal(length(select_segreg(segre, distorted = T, numbers = T)), sum(segre_tab$`p-value` < 0.05/length(segre_tab$Marker)))
    expect_equal(length(select_segreg(segre, distorted = T, threshold = 0.01, numbers = T)), sum(segre_tab$`p-value` < 0.01/length(segre_tab$Marker)))
  }
  
  check_dist("onemap_example_out", c(12,8,8,2))
  check_dist("onemap_example_f2", c(36,30))
  check_dist("onemap_example_bc", c(67))
  check_dist("onemap_example_riself", c(68))
})

test_that("number of bins",{
  check_bins <- function(example_data, n.mar){
    eval(bquote(data(.(example_data))))
    bins <- eval(bquote(find_bins(get(.(example_data)))))
    onemap_bins <- eval(bquote(create_data_bins(input.obj = get(.(example_data)), bins)))
    eval(bquote(expect_equal(check_data(onemap_bins),0)))
    eval(bquote(expect_equal(onemap_bins$n.mar, .(n.mar))))
  }
  check_bins("vcf_example_f2", 24)
  check_bins("vcf_example_out", 23)
  check_bins("vcf_example_bc", 25)
  check_bins("vcf_example_riself",25)
})

test_that("number of missing data",{
  check_missing <- function(example_data, n.mar){
    eval(bquote(data(.(example_data))))
    onemap_mis <- eval(bquote(filter_missing(get(.(example_data)), 0.5)))
    eval(bquote(expect_equal(check_data(onemap_mis), 0)))
    eval(bquote(expect_equal(onemap_mis$n.mar, .(n.mar))))
  }
  check_missing("vcf_example_f2", 25)
  check_missing("onemap_example_riself", 64)
  check_missing("onemap_example_out", 30)
  check_missing("onemap_example_bc", 67)
})

test_that("number of repeated ID markers",{
  check_dupli <- function(example_data, n.mar){
    eval(bquote(data(.(example_data))))
    onemap_dupli <- eval(bquote(rm_dupli_mks(get(.(example_data)))))
    eval(bquote(expect_equal(check_data(onemap_dupli), 0)))
    eval(bquote(expect_equal(onemap_dupli$n.mar, .(n.mar))))
  }
  
  check_dupli("vcf_example_f2", 25)
  check_dupli("onemap_example_riself", 68)
  check_dupli("onemap_example_out", 30)
  check_dupli("onemap_example_bc", 67)
  
})
