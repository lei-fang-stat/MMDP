test_that("alr_transform drops the reference and centres", {
  x <- matrix(c(1, 2, 4, 2, 2, 2, 8, 4, 1), 3, dimnames = list(NULL, c("a", "b", "c")))
  out <- alr_transform(x, "b")
  expect_equal(colnames(out), c("a", "c"))
  expect_equal(unname(colMeans(out)), c(0, 0))
  expect_equal(unname(alr_transform(x, 2, center = FALSE)[, 1]), log(c(1, 2, 4) / 2))
  expect_error(alr_transform(replace(x, 1, 0), 1), "strictly positive")
})

test_that("prepare_microbiome aligns taxa across cohorts", {
  sim <- simulate_cohorts(n_taxa = 6)
  tgt <- sim$target_counts[, c(1, 2, 3, 4, 5)]
  prep <- prepare_microbiome(sim$ref_counts, tgt, "taxon1")
  expect_equal(prep$taxa, paste0("taxon", 2:5))
  expect_identical(colnames(prep$X_ref), colnames(prep$X_target))
})

test_that("flhd_test reproduces phd::FLhd with CV ridge", {
  skip_if_not_installed("phd")
  set.seed(3)
  X <- matrix(rnorm(40 * 4), 40)
  y <- X[, 1] + rnorm(40)
  set.seed(10)
  p_ref <- phd::FLhd(y, X, nperm = 300)
  set.seed(10)
  p_new <- flhd_test(y, X, nperm = 300, lambda = "lambda.min")
  expect_equal(as.numeric(p_new), p_ref)
})

test_that("flhd_test uses REML with a single nuisance covariate", {
  set.seed(4)
  X <- matrix(rnorm(50 * 2), 50)
  y <- X[, 1] + X[, 2] + rnorm(50)
  p <- suppressWarnings(flhd_test(y, X, nperm = 200))
  expect_equal(attr(p, "lambda_method"), "reml")
  expect_gt(attr(p, "lambda"), 0)
  expect_lt(as.numeric(p), 0.05)
})

test_that("REML ridge works with several nuisance covariates", {
  set.seed(5)
  X <- matrix(rnorm(50 * 4), 50)
  y <- X[, 1] + rnorm(50)
  p <- suppressWarnings(flhd_test(y, X, nperm = 200, lambda = "reml"))
  expect_true(p >= 0 && p <= 1.01)
})

test_that("perm_cor_test detects correlation", {
  set.seed(6)
  x <- rnorm(60)
  expect_lt(perm_cor_test(x, x + rnorm(60, sd = 0.5), nperm = 500), 0.01)
})

test_that("full pipeline runs and captures errors", {
  sim <- simulate_cohorts(n_taxa = 10)
  prep <- prepare_microbiome(sim$ref_counts, sim$target_counts, 10)
  res <- suppressWarnings(
    test_metabolites(sim$metabolites, prep$X_ref, prep$X_target, sim$y, nperm = 200)
  )
  expect_equal(nrow(res), 3)
  expect_true(all(is.na(res$error)))
  expect_lt(res$p_value[1], 0.05)

  bad <- cbind(sim$metabolites, broken = c(NA, sim$metabolites[-1, 1]))
  res2 <- suppressWarnings(
    test_metabolites(bad, prep$X_ref, prep$X_target, sim$y, nperm = 200)
  )
  expect_match(res2$error[4], "missing")
})

test_that("reference_sensitivity loops over references", {
  sim <- simulate_cohorts(n_taxa = 6, n_metab = 1)
  res <- suppressWarnings(reference_sensitivity(
    sim$ref_counts, sim$target_counts, sim$metabolites, sim$y,
    references = c(1, 6), nperm = 100
  ))
  expect_equal(res$reference, c("taxon1", "taxon6"))
})

test_that("quantile_aggregate", {
  expect_equal(quantile_aggregate(c(0.01, 0.02, 0.5, 0.9)), min(1, quantile(c(0.01, 0.02, 0.5, 0.9) / 0.5, 0.5, names = FALSE)))
  expect_lte(quantile_aggregate(runif(20), gamma = c(0.25, 0.5)), 1)
})
