#' Test one predicted metabolite for association with an outcome
#'
#' Runs the full procedure for one metabolite:
#' 1. Standardise the metabolite and fit a cross-validated lasso of it on the
#'    reference microbiome (no intercept).
#' 2. Predict the metabolite in the target cohort from the target microbiome.
#' 3. Fit a cross-validated lasso of the outcome on the predicted metabolite and
#'    the target microbiome, after centring both.
#' 4. Test the predicted metabolite while adjusting for the taxa that step 3
#'    selected. With one taxon the test is [flhd_test()] with a REML ridge
#'    penalty. With two or more it is [flhd_test()] with a CV ridge penalty.
#'    With none it is [perm_cor_test()].
#'
#' @param m Numeric vector of metabolite levels in the reference cohort.
#' @param X_ref ALR microbiome matrix of the reference cohort (rows match `m`).
#' @param X_target ALR microbiome matrix of the target cohort (rows match `y`),
#'   with the same taxa in the same order as `X_ref`.
#' @param y Numeric outcome in the target cohort (e.g. centred log outcome).
#' @param nperm Number of permutations for [flhd_test()].
#' @param nperm_cor Number of permutations for [perm_cor_test()].
#' @param lambda,flip,nfolds,statistic Passed to [flhd_test()].
#' @param penalty_metab Lasso penalty factor for the predicted metabolite in
#'   step 3. The default 1 penalises it like every taxon. Use 0 to always keep
#'   it in the model.
#' @param seed If not `NULL`, `set.seed(seed)` is called before step 1 so each
#'   metabolite's result can be reproduced. This changes the global RNG state.
#'
#' @return An object of class `"metab_test"`: a list with `p_value`, `method`,
#'   `theta` (lasso coefficient of the predicted metabolite), `gamma` (lasso
#'   coefficients of the metabolite on the taxa), `selected_taxa` (taxa that the
#'   outcome lasso selected) and `lambda_ridge`.
#' @export
test_metabolite <- function(m, X_ref, X_target, y, nperm = 20000, nperm_cor = 5000,
                            lambda = "auto", flip = FALSE, nfolds = 10,
                            statistic = "partialcor", penalty_metab = 1, seed = 1234) {
  X_ref <- as.matrix(X_ref)
  X_target <- as.matrix(X_target)
  y <- as.numeric(y)
  if (ncol(X_ref) != ncol(X_target)) stop("X_ref and X_target must have the same taxa")
  if (!is.null(colnames(X_ref)) && !is.null(colnames(X_target)) &&
      !identical(colnames(X_ref), colnames(X_target))) {
    stop("X_ref and X_target columns are not aligned; see prepare_microbiome()")
  }
  if (length(m) != nrow(X_ref)) stop("length(m) must equal nrow(X_ref)")
  if (length(y) != nrow(X_target)) stop("length(y) must equal nrow(X_target)")
  if (anyNA(m) || anyNA(y)) stop("m and y must not contain missing values")

  if (!is.null(seed)) set.seed(seed)

  # Step 1: metabolite ~ microbiome in the reference cohort
  M <- scale(as.matrix(as.numeric(m)))
  cv_m <- glmnet::cv.glmnet(X_ref, M, intercept = FALSE, family = "gaussian")
  fit_m <- glmnet::glmnet(X_ref, M, lambda = cv_m$lambda.min, intercept = FALSE,
                          family = "gaussian")
  gamma <- fit_m$beta[, 1]
  gamma[abs(gamma) < 1e-10] <- 0
  names(gamma) <- colnames(X_ref)

  if (all(gamma == 0)) {
    return(new_metab_test(NA_real_, "none: metabolite not predictable", NA_real_,
                          gamma, character(0), NA_real_))
  }

  # Step 2: predicted metabolite in the target cohort
  m_hat <- X_target %*% gamma
  X_comb <- scale(cbind(m_hat, X_target), scale = FALSE)

  # Step 3: outcome ~ predicted metabolite + microbiome
  pf <- c(penalty_metab, rep(1, ncol(X_target)))
  cv_y <- glmnet::cv.glmnet(X_comb, y, intercept = FALSE, penalty.factor = pf,
                            family = "gaussian")
  fit_y <- glmnet::glmnet(X_comb, y, intercept = FALSE, lambda = cv_y$lambda.min,
                          penalty.factor = pf, family = "gaussian")
  theta <- fit_y$beta[1]
  sel <- which(fit_y$beta[-1] != 0)
  sel_taxa <- if (is.null(colnames(X_target))) as.character(sel) else colnames(X_target)[sel]

  # Step 4: test the predicted metabolite, adjusting for the selected taxa
  if (length(sel) == 0) {
    p <- perm_cor_test(X_comb[, 1], y, nperm = nperm_cor)
    method <- "perm_cor"
    lam <- NA_real_
  } else {
    p <- flhd_test(y, X_comb[, c(1, 1 + sel), drop = FALSE], nperm = nperm,
                   lambda = lambda, flip = flip, nfolds = nfolds, statistic = statistic)
    method <- paste0("flhd_", attr(p, "lambda_method"))
    lam <- attr(p, "lambda")
    p <- as.numeric(p)
  }
  new_metab_test(p, method, theta, gamma, sel_taxa, lam)
}

new_metab_test <- function(p, method, theta, gamma, selected_taxa, lambda_ridge) {
  structure(list(p_value = p, method = method, theta = theta, gamma = gamma,
                 selected_taxa = selected_taxa, lambda_ridge = lambda_ridge),
            class = "metab_test")
}

#' @export
print.metab_test <- function(x, ...) {
  cat("Predicted-metabolite association test\n")
  cat("  p-value:          ", format(x$p_value, digits = 4), "\n")
  cat("  method:           ", x$method, "\n")
  cat("  theta (lasso):    ", format(x$theta, digits = 4), "\n")
  cat("  taxa predicting M:", sum(x$gamma != 0), "\n")
  cat("  taxa adjusted for:", length(x$selected_taxa), "\n")
  invisible(x)
}

#' Test every metabolite in a matrix
#'
#' Calls [test_metabolite()] on each column of `metabolites`. If a metabolite
#' raises an error, it is recorded in the `error` column and the remaining
#' metabolites still run.
#'
#' @param metabolites Matrix or data frame of metabolites in the reference
#'   cohort, with one column per metabolite and rows matching `X_ref`.
#' @param X_ref,X_target,y See [test_metabolite()].
#' @param ... Further arguments passed to [test_metabolite()].
#' @param verbose Logical; print progress?
#'
#' @return A data frame with one row per metabolite and the columns
#'   `metabolite`, `p_value`, `method`, `theta`, `n_taxa_metabolite`,
#'   `n_taxa_outcome`, `lambda_ridge` and `error`.
#' @export
test_metabolites <- function(metabolites, X_ref, X_target, y, ..., verbose = FALSE) {
  metabolites <- as.data.frame(metabolites)
  ids <- colnames(metabolites)
  if (is.null(ids)) ids <- paste0("M", seq_along(metabolites))

  rows <- lapply(seq_along(metabolites), function(i) {
    if (verbose) message("[", i, "/", length(ids), "] ", ids[i])
    res <- tryCatch(test_metabolite(metabolites[[i]], X_ref, X_target, y, ...),
                    error = function(e) e)
    if (inherits(res, "error")) {
      return(data.frame(metabolite = ids[i], p_value = NA_real_, method = NA_character_,
                        theta = NA_real_, n_taxa_metabolite = NA_integer_,
                        n_taxa_outcome = NA_integer_, lambda_ridge = NA_real_,
                        error = conditionMessage(res)))
    }
    data.frame(metabolite = ids[i], p_value = res$p_value, method = res$method,
               theta = res$theta, n_taxa_metabolite = sum(res$gamma != 0),
               n_taxa_outcome = length(res$selected_taxa),
               lambda_ridge = res$lambda_ridge, error = NA_character_)
  })
  do.call(rbind, rows)
}

#' Repeat the analysis for each choice of ALR reference taxon
#'
#' For each reference taxon, prepares both cohorts with [prepare_microbiome()]
#' and runs [test_metabolites()]. Runs can be spread across cores with
#' [parallel::mclapply()] (not on Windows). The seed is reset for every
#' metabolite, so results do not depend on `mc.cores`.
#'
#' @param ref_counts,target_counts Abundance matrices of the reference and
#'   target cohorts (strictly positive). The rows of `target_counts` must match
#'   `y`.
#' @param metabolites Metabolite matrix of the reference cohort.
#' @param y Outcome in the target cohort.
#' @param references Reference taxa to try, as indices or names. Defaults to
#'   every taxon.
#' @param mc.cores Number of cores.
#' @param out_dir If not `NULL`, each reference's results are also written to
#'   `out_dir/reference<r>.pvalue.txt`.
#' @param ... Further arguments passed to [test_metabolite()].
#'
#' @return A data frame of [test_metabolites()] results for every reference,
#'   with an added `reference` column.
#' @export
reference_sensitivity <- function(ref_counts, target_counts, metabolites, y,
                                  references = seq_len(ncol(ref_counts)),
                                  mc.cores = 1L, out_dir = NULL, ...) {
  run_one <- function(r) {
    prep <- prepare_microbiome(ref_counts, target_counts, r)
    res <- test_metabolites(metabolites, prep$X_ref, prep$X_target, y, ...)
    res$reference <- prep$reference
    if (!is.null(out_dir)) {
      utils::write.table(res, file.path(out_dir, paste0("reference", r, ".pvalue.txt")),
                         sep = "\t")
    }
    res
  }
  out <- if (mc.cores > 1) {
    parallel::mclapply(references, run_one, mc.cores = mc.cores)
  } else {
    lapply(references, run_one)
  }
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    warning("reference(s) ", paste(references[failed], collapse = ", "), " failed: ",
            paste(unique(vapply(out[failed], as.character, "")), collapse = "; "))
  }
  do.call(rbind, out[!failed])
}

#' Quantile aggregation of p-values
#'
#' Combines p-values from repeated analyses, such as several reference taxa or
#' several sample splits, using the method of Meinshausen, Meier and Buhlmann
#' (2009). With one `gamma` it returns `min(1, quantile(p / gamma, gamma))`.
#' With several `gamma` values it takes the minimum over them and multiplies
#' by the correction factor `1 - log(min(gamma))`.
#'
#' @param p Numeric vector of p-values (`NA`s are dropped).
#' @param gamma Quantile level(s) in (0, 1).
#' @return A single aggregated p-value.
#' @references Meinshausen, N., Meier, L. and Buhlmann, P. (2009). P-values for
#'   high-dimensional regression. *JASA*, 104(488), 1671-1681.
#' @export
quantile_aggregate <- function(p, gamma = 0.5) {
  p <- p[!is.na(p)]
  if (length(p) == 0) return(NA_real_)
  q <- vapply(gamma, function(g) min(1, stats::quantile(p / g, g, names = FALSE)), numeric(1))
  if (length(gamma) == 1) return(q)
  min(1, (1 - log(min(gamma))) * min(q))
}
