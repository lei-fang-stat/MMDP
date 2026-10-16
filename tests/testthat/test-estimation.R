test_that("stage1_estimation returns one gamma row per reference", {
  sim <- simulate_cohorts(n_taxa = 8)
  g <- stage1_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts,
                         references = c(1, 8))
  expect_equal(dim(g), c(2, 8))
  expect_equal(rownames(g), c("taxon1", "taxon8"))
  expect_equal(unname(c(g[1, 1], g[2, 8])), c(1, 1))
  expect_false(anyNA(g))

  g2 <- stage1_estimation(sim$metabolites[, 1], sim$ref_counts,
                          sim$target_counts[, -5], references = 1)
  expect_true(is.na(g2[1, 5]))
})

test_that("stage2_estimation returns correlations, coefficients and theta", {
  sim <- simulate_cohorts(n_taxa = 8)
  ids <- rownames(sim$target_counts)[1:60]
  th <- stage2_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts,
                          sim$y[1:60], target_ids = ids, references = c(2, 8),
                          m_name = "metab1")
  expect_equal(dim(th), c(2, 9))
  expect_equal(colnames(th)[1], "metab1")
  expect_equal(unname(c(th[1, 3], th[2, 9])), c(1, 1))
  expect_true(all(abs(th[, 1]) <= 1))
  expect_length(attr(th, "theta"), 2)
  expect_error(stage2_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts,
                                 sim$y[1:60]), "one value per target sample")
})

