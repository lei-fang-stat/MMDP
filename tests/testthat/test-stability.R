toy_estimates <- function() {
  taxa <- paste0("t", 1:5)
  g <- matrix(0, 5, 5, dimnames = list(taxa, taxa))
  th <- matrix(0, 5, 6, dimnames = list(taxa, c("metab", taxa)))
  diag(g) <- 1
  th[cbind(1:5, 2:6)] <- 1
  th[, 1] <- 0.3
  g[-1, 1] <- 0.5            # t1: selected in every stage 1 fit
  g[-2, 2] <- c(0.2, 0, 0.1, 0.4)  # t2: 3 of 4 in stage 1
  th[-1, 2] <- 0.2           # t1: also in every stage 2 fit
  th[-3, 4] <- 0.7           # t3: stage 2 only
  list(g = g, th = th)
}

test_that("stability_selection classifies taxa into the three roles", {
  e <- toy_estimates()
  r <- stability_selection(e$g, e$th, threshold = 0.8)
  expect_equal(r$both, "t1")
  expect_equal(r$stage1_only, character(0))
  expect_equal(r$stage2_only, "t3")
  expect_equal(r$frequencies$freq_stage1[1:2], c(1, 0.75))
  r2 <- stability_selection(e$g, e$th, threshold = 0.7)
  expect_equal(r2$stage1_only, "t2")
})

test_that("stability_selection matches the original frequency formula", {
  sim <- simulate_cohorts(n_taxa = 8)
  g <- stage1_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts)
  th <- stage2_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts, sim$y)
  r <- stability_selection(g, th, threshold = 0.5)
  a <- g[complete.cases(g), ]
  b <- th[complete.cases(th), ]
  expect_equal(unname(r$frequencies$freq_stage1), unname((colSums(a != 0) - 1) / (nrow(a) - 1)))
  expect_equal(unname(r$frequencies$freq_stage2), unname(((colSums(b != 0) - 1) / (nrow(a) - 1))[-1]))
})

test_that("stability_selection reads written files, including unnamed rows", {
  e <- toy_estimates()
  dir <- tempfile(); dir.create(dir)
  f1 <- file.path(dir, "g.txt"); f2 <- file.path(dir, "th.txt")
  write.table(unname_rows(e$g), f1, sep = "\t"); write.table(unname_rows(e$th), f2, sep = "\t")
  expect_equal(stability_selection(f1, f2)$both, "t1")
})

test_that("NA stage 2 rows count as not selected", {
  e <- toy_estimates()
  e$th[2, ] <- NA
  r <- stability_selection(e$g, e$th, threshold = 0.5)
  expect_equal(r$frequencies$freq_stage2[1], 0.75)
  expect_error(stability_selection(e$g, e$th, threshold = 1), "threshold")
})
