#' Stability selection of microbes across stage 1 and stage 2
#'
#' Uses the stage 1 coefficients (`gamma`) and the stage 2 taxon coefficients
#' (`theta[-1]`, i.e. without the first column) estimated for every choice of
#' ALR reference taxon. For each taxon it computes the selection frequency: the
#' share of reference taxa (other than the taxon itself) for which its
#' coefficient is nonzero. Taxa whose frequency exceeds `threshold` count as
#' selected in that stage, and each taxon is then classified as selected
#' in stage 1 only, in both stages, or in stage 2 only.
#'
#' Both stages use the same denominator: the number of stage 1 fits where the
#' taxon was a candidate. A reference taxon for which stage 2 is `NA` (the
#' metabolite could not be predicted) therefore counts as "not selected" in
#' stage 2. With complete estimates this gives the same frequencies as
#' `(colSums(a != 0) - 1) / (nrow(a) - 1)` in the original analysis.
#'
#' @param gamma_all Stage 1 estimates: the matrix returned by
#'   [stage1_estimation()], or the path to a `.gamma.all.txt` file written from
#'   it with `write.table(sep = "\t")`.
#' @param theta_all Stage 2 estimates: the matrix returned by
#'   [stage2_estimation()], or the path to a `.theta.all.txt` file. Its first
#'   column (the correlation) is dropped.
#' @param threshold Selection-frequency threshold in (0, 1). A taxon is
#'   selected in a stage when its frequency is strictly greater than it.
#'
#' @return An object of class `"microbe_roles"`: a list with
#'   * `stage1_only`, `both`, `stage2_only`: taxon names in each category;
#'   * `frequencies`: a data frame with `taxon`, `freq_stage1`, `freq_stage2`
#'     and `role` (`"stage1_only"`, `"both"`, `"stage2_only"` or `"none"`);
#'   * `threshold`.
#'
#' Each row of the estimates must say which reference taxon it used, either
#' through row names (as returned by the estimation functions) or, for files
#' from the original scripts, by row `k` using reference taxon `k`.
#' @examples
#' sim <- simulate_cohorts(n_taxa = 8)
#' g <- stage1_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts)
#' th <- stage2_estimation(sim$metabolites[, 1], sim$ref_counts, sim$target_counts, sim$y)
#' stability_selection(g, th, threshold = 0.8)
#' @export
stability_selection <- function(gamma_all, theta_all, threshold = 0.8) {
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
    stop("threshold must be a single number between 0 and 1")
  }
  gamma_all <- read_estimates(gamma_all)
  theta_all <- read_estimates(theta_all)
  theta_all <- theta_all[, -1, drop = FALSE]
  if (!identical(colnames(gamma_all), colnames(theta_all))) {
    stop("gamma_all and theta_all[, -1] must have the same taxa in the same order")
  }
  if (nrow(gamma_all) != nrow(theta_all) ||
      (!is.null(rownames(gamma_all)) && !identical(rownames(gamma_all), rownames(theta_all)))) {
    stop("gamma_all and theta_all must have one row per reference taxon, in the same order")
  }

  # candidate[r, j]: taxon j was in the stage 1 fit for reference r
  ref_col <- reference_columns(gamma_all)
  is_ref <- matrix(FALSE, nrow(gamma_all), ncol(gamma_all))
  is_ref[cbind(seq_len(nrow(gamma_all)), ref_col)] <- TRUE
  candidate <- !is.na(gamma_all) & !is_ref
  n_candidate <- colSums(candidate)

  freq1 <- colSums(candidate & gamma_all != 0) / n_candidate
  freq2 <- colSums(candidate & !is.na(theta_all) & theta_all != 0) / n_candidate
  s1 <- !is.na(freq1) & freq1 > threshold
  s2 <- !is.na(freq2) & freq2 > threshold

  taxa <- colnames(gamma_all)
  role <- ifelse(s1 & s2, "both", ifelse(s1, "stage1_only", ifelse(s2, "stage2_only", "none")))
  structure(list(
    stage1_only = taxa[s1 & !s2],
    both = taxa[s1 & s2],
    stage2_only = taxa[s2 & !s1],
    frequencies = data.frame(taxon = taxa, freq_stage1 = unname(freq1),
                             freq_stage2 = unname(freq2), role = role),
    threshold = threshold
  ), class = "microbe_roles")
}

#' @export
print.microbe_roles <- function(x, ...) {
  show <- function(label, taxa) {
    cat(label, " (", length(taxa), "): ",
        if (length(taxa)) paste(taxa, collapse = ", ") else "-", "\n", sep = "")
  }
  cat("Stability selection, threshold", x$threshold, "\n")
  show("  Stage 1 only", x$stage1_only)
  show("  Both stages ", x$both)
  show("  Stage 2 only", x$stage2_only)
  invisible(x)
}

read_estimates <- function(x) {
  if (is.character(x) && length(x) == 1) {
    x <- utils::read.table(x, sep = "\t", check.names = FALSE)
  }
  x <- as.matrix(x)
  if (is.null(colnames(x))) stop("the estimates need taxon names as column names")
  x
}

# Column index of the reference taxon used in each row: matched by row name,
# or row k = taxon k for the original scripts' files (rows named 1, 2, ...).
reference_columns <- function(est) {
  idx <- match(rownames(est), colnames(est))
  if (!is.null(rownames(est)) && !anyNA(idx)) return(idx)
  if (nrow(est) == ncol(est)) return(seq_len(nrow(est)))
  stop("cannot tell which reference taxon each row used; ",
       "name the rows after their reference taxa")
}
