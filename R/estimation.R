#' Stage 1 estimation: lasso coefficients of a metabolite on the microbiome
#'
#' For each choice of ALR reference taxon, standardizes the ALR-transformed
#' microbiome of the reference cohort (centered and scaled) and fits a
#' cross-validated lasso of the standardized metabolite on it (no intercept,
#' `lambda.min`). The fitted coefficients (`gamma`) show which taxa predict the
#' metabolite and how stable that selection is across reference taxa.
#'
#' @param m Numeric vector of metabolite levels in the reference cohort.
#' @param ref_counts Abundance matrix of the reference cohort (strictly
#'   positive, taxon names as column names, rows matching `m`).
#' @param target_counts Optional abundance matrix of the target cohort. If it
#'   is given, only taxa present in both cohorts are used, as in the original
#'   analysis.
#' @param references Reference taxa to use, as indices or names. Defaults to
#'   every taxon.
#' @param seed If not `NULL`, `set.seed(seed)` is called before each lasso fit.
#'
#' @return A matrix with one row per reference taxon and one column per taxon
#'   in `ref_counts`, holding the lasso coefficients. In each row the reference
#'   taxon's own cell is set to 1, as a marker. Taxa not shared with
#'   `target_counts` are `NA`.
#' @examples
#' sim <- simulate_cohorts(n_taxa = 8)
#' stage1_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts,
#'                   references = 1:3)
#' @export
stage1_estimation <- function(m, ref_counts, target_counts = NULL,
                              references = seq_len(ncol(ref_counts)), seed = 1234) {
  ref_counts <- as.matrix(ref_counts)
  if (is.null(colnames(ref_counts))) stop("ref_counts needs taxon names as column names")
  if (length(m) != nrow(ref_counts)) stop("length(m) must equal nrow(ref_counts)")
  if (!is.null(target_counts)) target_counts <- as.matrix(target_counts)

  gamma_all <- matrix(NA_real_, nrow = length(references), ncol = ncol(ref_counts),
                      dimnames = list(NULL, colnames(ref_counts)))
  ref_names <- character(length(references))

  for (k in seq_along(references)) {
    ref_idx <- resolve_reference(ref_counts, references[k])
    ref_names[k] <- colnames(ref_counts)[ref_idx]
    X_ref <- scale(alr_transform(ref_counts, ref_idx, center = FALSE))
    taxa <- colnames(X_ref)
    if (!is.null(target_counts)) {
      tgt_idx <- resolve_reference(target_counts, references[k])
      taxa <- intersect(colnames(target_counts)[-tgt_idx], taxa)
    }
    X_ref <- X_ref[, match(taxa, colnames(X_ref)), drop = FALSE]

    if (!is.null(seed)) set.seed(seed)
    M <- scale(as.matrix(as.numeric(m)))
    cv_m <- glmnet::cv.glmnet(X_ref, M, alpha = 1, intercept = FALSE, family = "gaussian")
    fit_m <- glmnet::glmnet(X_ref, M, lambda = cv_m$lambda.min, intercept = FALSE,
                            family = "gaussian")
    gamma <- fit_m$beta[, 1]
    gamma[abs(gamma) < 1e-10] <- 0

    gamma_all[k, taxa] <- gamma
    gamma_all[k, ref_idx] <- 1
  }
  rownames(gamma_all) <- ref_names
  gamma_all
}

#' Stage 2 estimation: lasso of the outcome on the predicted metabolite
#'
#' For each choice of ALR reference taxon:
#' 1. Centres the ALR-transformed microbiome of both cohorts.
#' 2. Fits a cross-validated lasso of the standardised metabolite on the
#'    reference microbiome, then predicts the metabolite in the target cohort.
#' 3. Fits a cross-validated lasso of the outcome on the predicted metabolite
#'    and the target microbiome (no intercept, `lambda.min`).
#'
#' If the stage 1 lasso selects no taxa, the metabolite cannot be predicted and
#' that row is `NA`. An error for one reference taxon is reported with
#' [message()], leaves that row `NA`, and does not stop the other references.
#'
#' @param m Numeric vector of metabolite levels in the reference cohort.
#' @param ref_counts Abundance matrix of the reference cohort (strictly
#'   positive, taxon names as column names, rows matching `m`).
#' @param target_counts Abundance matrix of the target cohort.
#' @param y Numeric outcome in the target cohort, e.g. the centred log outcome.
#' @param target_ids Optional row names of `target_counts` matching `y`, such as
#'   the subjects with an observed outcome. Centring uses every row of
#'   `target_counts`, and the rows are then matched to `target_ids`, as in the
#'   original analysis. If `NULL`, the rows of `target_counts` must already
#'   match `y`.
#' @param references Reference taxa to use, as indices or names. Defaults to
#'   every taxon.
#' @param m_name Name of the first output column.
#' @param seed If not `NULL`, `set.seed(seed)` is called before each fit.
#'
#' @return A matrix with one row per reference taxon. The first column (named
#'   `m_name`) is the correlation between the predicted metabolite and `y`.
#'   The remaining columns, one per taxon in `ref_counts`, hold the outcome
#'   lasso coefficients. In each row the reference taxon's own cell is set to
#'   1, as a marker. The lasso coefficient of the predicted metabolite itself
#'   is returned in the `"theta"` attribute, one value per reference.
#' @examples
#' sim <- simulate_cohorts(n_taxa = 8)
#' stage2_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts,
#'                   sim$y, references = 1:3)
#' @export
stage2_estimation <- function(m, ref_counts, target_counts, y, target_ids = NULL,
                              references = seq_len(ncol(ref_counts)),
                              m_name = "metabolite", seed = 1234) {
  ref_counts <- as.matrix(ref_counts)
  target_counts <- as.matrix(target_counts)
  y <- as.numeric(y)
  if (length(m) != nrow(ref_counts)) stop("length(m) must equal nrow(ref_counts)")
  rows <- seq_len(nrow(target_counts))
  if (!is.null(target_ids)) {
    rows <- match(target_ids, rownames(target_counts))
    if (anyNA(rows)) stop("some target_ids are not in rownames(target_counts)")
  }
  if (length(rows) != length(y)) stop("y must have one value per target sample")

  theta_all <- matrix(NA_real_, nrow = length(references), ncol = ncol(ref_counts) + 1,
                      dimnames = list(NULL, c(m_name, colnames(ref_counts))))
  theta <- rep(NA_real_, length(references))
  ref_names <- character(length(references))

  for (k in seq_along(references)) {
    ref_names[k] <- colnames(ref_counts)[resolve_reference(ref_counts, references[k])]
    tryCatch({
      prep <- prepare_microbiome(ref_counts, target_counts, references[k])
      X_ref <- prep$X_ref
      X_target <- prep$X_target[rows, , drop = FALSE]

      if (!is.null(seed)) set.seed(seed)
      M <- scale(as.matrix(as.numeric(m)))
      cv_m <- glmnet::cv.glmnet(X_ref, M, intercept = FALSE, family = "gaussian")
      fit_m <- glmnet::glmnet(X_ref, M, lambda = cv_m$lambda.min, intercept = FALSE,
                              family = "gaussian")
      gamma <- fit_m$beta[, 1]
      selected <- any(gamma != 0)
      gamma[abs(gamma) < 1e-10] <- 0

      if (selected) {
        X_comb <- cbind(X_target %*% gamma, X_target)
        penalty <- rep(1, ncol(X_comb))
        cv_y <- glmnet::cv.glmnet(X_comb, y, intercept = FALSE, penalty.factor = penalty,
                                  family = "gaussian")
        fit_y <- glmnet::glmnet(X_comb, y, intercept = FALSE, lambda = cv_y$lambda.min,
                                penalty.factor = penalty, family = "gaussian")
        beta <- fit_y$beta[, 1]
        theta[k] <- beta[1]
        theta_all[k, prep$taxa] <- beta[-1]
        theta_all[k, 1 + resolve_reference(ref_counts, references[k])] <- 1
        theta_all[k, 1] <- cor(X_comb[, 1], y)
      }
    }, error = function(e) {
      message("reference ", ref_names[k], " failed: ", conditionMessage(e))
    })
  }
  rownames(theta_all) <- ref_names
  names(theta) <- ref_names
  attr(theta_all, "theta") <- theta
  theta_all
}
