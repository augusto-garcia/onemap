# ── adjust_rf_ril ──────────────────────────────────────────────────────────────

test_that("adjust_rf_ril: riself expand/collapse are inverse operations", {
  r_obs <- c(0, 0.1, 0.25, 0.49)
  r_exp <- onemap:::adjust_rf_ril(r_obs, type = "riself", expand = TRUE)
  r_back <- onemap:::adjust_rf_ril(r_exp, type = "riself", expand = FALSE)
  expect_equal(r_back, r_obs, tolerance = 1e-10)
})

test_that("adjust_rf_ril: risib expand/collapse are inverse operations", {
  r_obs <- c(0, 0.05, 0.15, 0.3)
  r_exp <- onemap:::adjust_rf_ril(r_obs, type = "risib", expand = TRUE)
  r_back <- onemap:::adjust_rf_ril(r_exp, type = "risib", expand = FALSE)
  expect_equal(r_back, r_obs, tolerance = 1e-10)
})

test_that("adjust_rf_ril: r=0 stays 0 for both types and directions", {
  expect_equal(onemap:::adjust_rf_ril(0, "riself", expand = TRUE),  0)
  expect_equal(onemap:::adjust_rf_ril(0, "riself", expand = FALSE), 0)
  expect_equal(onemap:::adjust_rf_ril(0, "risib",  expand = TRUE),  0)
  expect_equal(onemap:::adjust_rf_ril(0, "risib",  expand = FALSE), 0)
})

test_that("adjust_rf_ril: expand always returns value larger than input for r > 0", {
  r_obs <- c(0.05, 0.1, 0.2, 0.4)
  expect_true(all(onemap:::adjust_rf_ril(r_obs, "riself", expand = TRUE) > r_obs))
  expect_true(all(onemap:::adjust_rf_ril(r_obs, "risib",  expand = TRUE) > r_obs))
})

test_that("adjust_rf_ril: unknown type raises error", {
  expect_error(onemap:::adjust_rf_ril(0.1, type = "unknown"))
})

test_that("adjust_rf_ril: known formula values for riself", {
  # r_expanded = r*2/(1+2*r)
  r <- 0.1
  expect_equal(onemap:::adjust_rf_ril(r, "riself", expand = TRUE),
               r * 2 / (1 + 2 * r))
})

test_that("adjust_rf_ril: known formula values for risib", {
  # r_expanded = r*4/(1+6*r)
  r <- 0.1
  expect_equal(onemap:::adjust_rf_ril(r, "risib", expand = TRUE),
               r * 4 / (1 + 6 * r))
})

# ── compare ───────────────────────────────────────────────────────────────────

test_that("compare: returns correct class for outcross data", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(12, 14, 15))
  cmp   <- onemap::compare(seq1)
  expect_s3_class(cmp, "compare")
})

test_that("compare: best.ord contains the marker indices", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(12, 14, 15))
  cmp   <- onemap::compare(seq1)
  expect_true(all(sort(cmp$best.ord[1, ]) == c(12, 14, 15)))
})

test_that("compare: n.best controls number of stored orders", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(12, 14, 15))
  cmp   <- onemap::compare(seq1, n.best = 5)
  # For 3 markers there are 3!/2 = 3 distinct orders; stored rows <= n.best
  expect_true(nrow(cmp$best.ord) <= 6)
})

test_that("compare: log-likelihoods are in decreasing order", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(12, 14, 15))
  cmp   <- onemap::compare(seq1)
  lod   <- cmp$best.ord.LOD
  lod   <- lod[!is.infinite(lod) & !is.na(lod)]
  expect_true(all(diff(lod) <= 0))
})

test_that("compare: works for inbred (f2) data", {
  data(onemap_example_f2)
  twopt <- rf_2pts(onemap_example_f2)
  seq1  <- make_seq(twopt, c(17, 26, 29))
  cmp   <- onemap::compare(seq1)
  expect_s3_class(cmp, "compare")
  expect_true(all(sort(cmp$best.ord[1, ]) == c(17, 26, 29)))
})

test_that("compare: print.compare runs without error", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(12, 14, 15))
  cmp   <- onemap::compare(seq1)
  expect_output(print(cmp))
})

# ── drop_marker ───────────────────────────────────────────────────────────────

test_that("drop_marker: reduces sequence length correctly", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  lgs   <- group(make_seq(twopt, "all"))
  lg1   <- make_seq(lgs, 1)
  orig_len <- length(lg1$seq.num)
  to_drop  <- lg1$seq.num[1:2]
  new_seq  <- drop_marker(lg1, to_drop)
  expect_equal(length(new_seq$seq.num), orig_len - 2)
})

test_that("drop_marker: dropped markers are absent from result", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  lgs   <- group(make_seq(twopt, "all"))
  lg1   <- make_seq(lgs, 1)
  to_drop <- lg1$seq.num[1:2]
  new_seq <- drop_marker(lg1, to_drop)
  expect_false(any(to_drop %in% new_seq$seq.num))
})

test_that("drop_marker: warns when marker is not in sequence", {
  data(onemap_example_out)
  twopt   <- rf_2pts(onemap_example_out)
  lgs     <- group(make_seq(twopt, "all"))
  lg1     <- make_seq(lgs, 1)
  absent  <- max(lg1$seq.num) + 999
  expect_warning(drop_marker(lg1, absent))
})

test_that("drop_marker: returns sequence class", {
  data(onemap_example_out)
  twopt   <- rf_2pts(onemap_example_out)
  lgs     <- group(make_seq(twopt, "all"))
  lg1     <- make_seq(lgs, 1)
  new_seq <- drop_marker(lg1, lg1$seq.num[1])
  expect_s3_class(new_seq, "sequence")
})

test_that("drop_marker: dropping all but one leaves single-marker sequence", {
  data(onemap_example_out)
  twopt   <- rf_2pts(onemap_example_out)
  lgs     <- group(make_seq(twopt, "all"))
  lg1     <- make_seq(lgs, 1)
  keep    <- lg1$seq.num[1]
  to_drop <- lg1$seq.num[-1]
  new_seq <- drop_marker(lg1, to_drop)
  expect_equal(new_seq$seq.num, keep)
})

test_that("drop_marker: error on non-sequence input", {
  expect_error(drop_marker("not_a_sequence", 1))
})

# ── filter_missing ────────────────────────────────────────────────────────────

test_that("filter_missing: returns onemap object", {
  data(onemap_example_out)
  filt <- filter_missing(onemap_example_out, threshold = 0.25)
  expect_s3_class(filt, "onemap")
})

test_that("filter_missing: no marker has missing rate above threshold (by markers)", {
  data(onemap_example_out)
  thresh <- 0.25
  filt   <- filter_missing(onemap_example_out, threshold = thresh, by = "markers")
  mis    <- apply(filt$geno, 2, function(x) sum(x == 0) / length(x))
  expect_true(all(mis <= thresh))
})

test_that("filter_missing: n.mar updated correctly", {
  data(onemap_example_out)
  filt <- filter_missing(onemap_example_out, threshold = 0.25)
  expect_equal(filt$n.mar, ncol(filt$geno))
})

test_that("filter_missing: no individual has missing rate above threshold (by individuals)", {
  data(onemap_example_out)
  thresh <- 0.5
  filt   <- filter_missing(onemap_example_out, threshold = thresh, by = "individuals")
  mis    <- apply(filt$geno, 1, function(x) sum(x == 0) / length(x))
  expect_true(all(mis <= thresh))
})

test_that("filter_missing: n.ind updated correctly when filtering by individuals", {
  data(onemap_example_out)
  filt <- filter_missing(onemap_example_out, threshold = 0.5, by = "individuals")
  expect_equal(filt$n.ind, nrow(filt$geno))
})

test_that("filter_missing: threshold=1 keeps all markers", {
  data(onemap_example_out)
  filt <- filter_missing(onemap_example_out, threshold = 1, by = "markers")
  expect_equal(filt$n.mar, onemap_example_out$n.mar)
})

test_that("filter_missing: threshold=0 removes all markers with any missing", {
  data(onemap_example_out)
  filt <- filter_missing(onemap_example_out, threshold = 0, by = "markers")
  mis  <- apply(filt$geno, 2, function(x) sum(x == 0) / length(x))
  expect_true(all(mis == 0))
})

test_that("filter_missing: error on invalid 'by' argument", {
  data(onemap_example_out)
  expect_error(filter_missing(onemap_example_out, by = "invalid"))
})

test_that("filter_missing: error on non-onemap input", {
  expect_error(filter_missing("not_onemap"))
})

# ── filter_prob ───────────────────────────────────────────────────────────────

test_that("filter_prob: returns onemap object", {
  data(onemap_example_out)
  obj  <- create_probs(onemap_example_out, global_error = 0.05)
  filt <- filter_prob(obj, threshold = 0.9)
  expect_s3_class(filt, "onemap")
})

test_that("filter_prob: genotypes below threshold are set to 0 (missing)", {
  data(onemap_example_out)
  obj    <- create_probs(onemap_example_out, global_error = 0.05)
  thresh <- 0.99  # very strict: most genotypes should become missing
  filt   <- filter_prob(obj, threshold = thresh)
  # Every remaining non-zero genotype must have had max prob >= threshold
  max_prob <- apply(filt$error, 1, max)
  non_missing_rows <- which(filt$geno != 0)
  if (length(non_missing_rows) > 0)
    expect_true(all(max_prob[non_missing_rows] >= thresh))
})

test_that("filter_prob: threshold=0 changes nothing", {
  data(onemap_example_out)
  obj  <- create_probs(onemap_example_out, global_error = 0.05)
  filt <- filter_prob(obj, threshold = 0)
  expect_equal(filt$geno, obj$geno)
})

# ── find_bins ─────────────────────────────────────────────────────────────────

test_that("find_bins: returns onemap_bin class (exact=TRUE)", {
  data(onemap_example_out)
  bins <- find_bins(onemap_example_out, exact = TRUE)
  expect_s3_class(bins, "onemap_bin")
})

test_that("find_bins: returns onemap_bin class (exact=FALSE)", {
  data(onemap_example_out)
  bins <- find_bins(onemap_example_out, exact = FALSE)
  expect_s3_class(bins, "onemap_bin")
})

test_that("find_bins: number of bins <= number of markers", {
  data(onemap_example_out)
  bins <- find_bins(onemap_example_out, exact = TRUE)
  expect_lte(length(bins$bins), onemap_example_out$n.mar)
})

test_that("find_bins: total markers across bins equals n.mar", {
  data(onemap_example_out)
  bins      <- find_bins(onemap_example_out, exact = TRUE)
  total_mks <- sum(sapply(bins$bins, nrow))
  expect_equal(total_mks, onemap_example_out$n.mar)
})

test_that("find_bins: exact=FALSE gives <= bins than exact=TRUE", {
  data(onemap_example_out)
  bins_exact    <- find_bins(onemap_example_out, exact = TRUE)
  bins_inexact  <- find_bins(onemap_example_out, exact = FALSE)
  expect_lte(length(bins_inexact$bins), length(bins_exact$bins))
})

test_that("find_bins: error on non-onemap input", {
  expect_error(find_bins("not_onemap"))
})

test_that("find_bins: error when fewer than 2 markers", {
  data(onemap_example_out)
  # create a minimal 1-marker object
  one_mk <- onemap_example_out
  one_mk$geno      <- one_mk$geno[, 1, drop = FALSE]
  one_mk$n.mar     <- 1
  one_mk$segr.type <- one_mk$segr.type[1]
  one_mk$segr.type.num <- one_mk$segr.type.num[1]
  expect_error(find_bins(one_mk))
})

test_that("print.onemap_bin: produces output", {
  data(onemap_example_out)
  bins <- find_bins(onemap_example_out, exact = TRUE)
  expect_output(print(bins))
})

test_that("print.onemap_bin: reports exact search correctly", {
  data(onemap_example_out)
  bins_exact   <- find_bins(onemap_example_out, exact = TRUE)
  bins_inexact <- find_bins(onemap_example_out, exact = FALSE)
  expect_output(print(bins_exact),   "exact")
  expect_output(print(bins_inexact), "non exact")
})

# ── marker_type ───────────────────────────────────────────────────────────────

test_that("marker_type: returns a data.frame", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(3, 6, 8, 12, 16, 25))
  mt    <- marker_type(seq1)
  expect_s3_class(mt, "data.frame")
})

test_that("marker_type: has columns Marker, Marker.name, Type", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  seq1  <- make_seq(twopt, c(3, 6, 8))
  mt    <- marker_type(seq1)
  expect_named(mt, c("Marker", "Marker.name", "Type"))
})

test_that("marker_type: number of rows equals number of markers in sequence", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  mks   <- c(3, 6, 8, 12, 16)
  seq1  <- make_seq(twopt, mks)
  mt    <- marker_type(seq1)
  expect_equal(nrow(mt), length(mks))
})

test_that("marker_type: Marker column matches seq.num", {
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  mks   <- c(3, 6, 8)
  seq1  <- make_seq(twopt, mks)
  mt    <- marker_type(seq1)
  expect_equal(mt$Marker, mks)
})

test_that("marker_type: works for f2 data and shows expected type labels", {
  data(onemap_example_f2)
  twopt <- rf_2pts(onemap_example_f2)
  lgs   <- group(make_seq(twopt, "all"))
  lg1   <- make_seq(lgs, 1)
  mt    <- marker_type(lg1)
  expect_s3_class(mt, "data.frame")
  known_types <- c("A.H.B", "C.A", "D.B")
  expect_true(all(mt$Type %in% c(known_types, "NA")))
})


test_that("marker_type: error on non-sequence input", {
  expect_error(marker_type("not_a_sequence"))
})

