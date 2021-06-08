context("test plot haplotypes")

data("simu_example_bc")

data("onemap_example_bc")

twopts <- rf_2pts(onemap_example_bc)

seq1 <- make_seq(twopts, "all")
LG1 <- group(seq1)
LG1 <- make_seq(LG1,1)

seq1 <- order_seq(LG1)
map_bc <- make_seq(seq1, "force")
rf_graph_table(map1)

probs1 <- progeny_haplotypes(map_bc)
plot(x = probs1, position = "split")

data("onemap_example_out")
twopts <- rf_2pts(onemap_example_out)

seq1 <- make_seq(twopts, "all")
LG1 <- group(seq1)
LG1 <- make_seq(LG1,1)

seq1 <- order_seq(LG1)
map1_out <- make_seq(seq1, "force")
rf_graph_table(map1)

probs_out <- progeny_haplotypes(map1_out)
plot(x = probs_out)

data("onemap_example_riself")

twopts <- rf_2pts(onemap_example_riself)

seq1 <- make_seq(twopts, "all")
LG1 <- group(seq1)
LG1 <- make_seq(LG1,1)

seq1 <- order_seq(LG1)
map_ri <- make_seq(seq1, "force")
rf_graph_table(map1)

probs1 <- progeny_haplotypes(map_ri, ind = 2)
plot(x = probs1, position = "split")

data("onemap_example_f2")

twopts <- rf_2pts(onemap_example_f2)

seq1 <- make_seq(twopts, "all")
LG1 <- group(seq1)
LG1 <- make_seq(LG1,1)

seq1 <- order_seq(LG1)
map_f2 <- make_seq(seq1, "force")
rf_graph_table(map1)

probs1 <- progeny_haplotypes(map_f2, ind = 2)
plot(x = probs1, position = "split")
