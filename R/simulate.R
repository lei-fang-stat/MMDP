#' Simulate a reference and a target cohort
#'
#' Produces toy data with the structure that [reference_sensitivity()] and
#' [prepare_microbiome()] expect. Abundances are strictly positive. Metabolites
#' depend on a few taxa, and the outcome depends on the first metabolite plus
#' one taxon.
#'
#' @param n_ref,n_target Sample sizes of the reference and target cohorts.
#' @param n_taxa Number of taxa.
#' @param n_metab Number of metabolites.
#' @param effect Effect of metabolite 1 on the outcome.
#' @param seed Random seed.
#' @return A list with `ref_counts`, `target_counts`, `metabolites` (reference
#'   cohort) and `y` (target cohort, log scale and centred).
#' @examples
#' sim <- simulate_cohorts()
#' str(sim, max.level = 1)
#' @export
simulate_cohorts <- function(n_ref = 100, n_target = 80, n_taxa = 20, n_metab = 3,
                             effect = 1, seed = 1) {
  set.seed(seed)
  taxa <- paste0("taxon", seq_len(n_taxa))
  gen_counts <- function(n) {
    z <- matrix(stats::rnorm(n * n_taxa), n, n_taxa)
    x <- exp(z) / rowSums(exp(z))
    dimnames(x) <- list(paste0("s", seq_len(n)), taxa)
    x
  }
  ref_counts <- gen_counts(n_ref)
  target_counts <- gen_counts(n_target)
  B <- matrix(0, n_taxa - 1, n_metab)
  for (k in seq_len(n_metab)) B[sample(n_taxa - 1, 3), k] <- stats::runif(3, 0.5, 1.5)
  alr_r <- log(ref_counts[, -n_taxa] / ref_counts[, n_taxa])
  alr_t <- log(target_counts[, -n_taxa] / target_counts[, n_taxa])
  metabolites <- alr_r %*% B + matrix(stats::rnorm(n_ref * n_metab), n_ref, n_metab)
  colnames(metabolites) <- paste0("metab", seq_len(n_metab))
  m1_target <- alr_t %*% B[, 1]
  y <- effect * m1_target + 0.5 * alr_t[, 2] + stats::rnorm(n_target)
  list(ref_counts = ref_counts, target_counts = target_counts,
       metabolites = metabolites, y = as.numeric(y - mean(y)))
}
