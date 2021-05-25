context("testing two-points estimations")

test_that("two-points tests",{
  expect_twopts_phases <- function(example_data, phase1 , phase2, phase3, phase4, sum_all){
    eval(bquote(data(.(example_data))))
    twopts <- eval(bquote(rf_2pts(get(.(example_data)))))
    expect_equal(check_twopts(twopts),0)
    eval(bquote(expect_equal(as.vector(twopts$analysis[[1]][1:4,1]), .(phase1), tolerance = 0.00001)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[2]][1:4,1]), .(phase2), tolerance = 0.00001)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[3]][1:4,1]), .(phase3), tolerance = 0.00001)))
    eval(bquote(expect_equal(as.vector(twopts$analysis[[4]][1:4,1]), .(phase4), tolerance = 0.00001)))
    eval(bquote(expect_equal(sum(twopts$analysis[[1]]), .(sum_all), tolerance = 0.000001)))
  }
  
  expect_twopts_phases("onemap_example_out", 
                       c(0, 0.1818151, 0.2954514, 0.5746628),
                       c(0, 0.8181849, 0.2954514, 0.5089259),
                       c(0, 0.1818151, 0.7045486, 0.4910741),
                       c(0, 0.8181849, 0.7045486, 0.4253372),
                       898.9572)
  
  expect_twopts_phases("onemap_example_f2", 
                       c(0, 0.4647862, 0.1621592, 0.7536202),
                       c(0, 0.5352138, 0.1621592, 0.2463798),
                       c(0, 0.4647862, 0.8378408, 0.7536202),
                       c(0, 0.5352138, 0.8378408, 0.2463798),
                       7324.314)
  
  expect_twopts <- function(example_data, values, sum_all){
    eval(bquote(data(.(example_data))))
    onemap.mis <- eval(bquote(filter_missing(get(.(example_data)), 0.25)))
    twopts <- eval(bquote(rf_2pts(onemap.mis)))
    expect_equal(check_twopts(twopts),0)
    eval(bquote(expect_equal(as.vector(twopts$analysis[1:4,1]), .(values), tolerance = 0.00001)))
    eval(bquote(expect_equal(sum(twopts$analysis), .(sum_all), tolerance = 0.000001)))
  }
  
  expect_twopts("onemap_example_bc",
                c(0.00000000, 0.03603604, 0.50000000, 0.50000000), 6734.126)
  
  expect_twopts("onemap_example_riself",
                c(0.00000000, 0.4974619, 0.5000000, 0.5000000), 5546.198)
  
  data("onemap_example_riself")
  expect_message(rf_2pts(onemap_example_riself), "We could not estimate all recombination fraction. Check if these markers have at least one genotype information or if they have segregation pattern deviation. We suggest filter_missing function to avoid excess of missing data, test_segregation and rm_mks argument.")
  expect_message(rf_2pts(onemap_example_riself, rm_mks = T), 
                 "Recombination fraction for 2 markers could not be estimated. They were removed from analysis. Check if these markers have at least one genotype information or if they have segregation pattern deviation. We suggest filter_missing function to avoid excess of missing data and test_segregation.")  
})
