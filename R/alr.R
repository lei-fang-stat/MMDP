#' Additive log-ratio (ALR) transform
#'
#' Divides every taxon by a reference taxon, takes logs and drops the
#' reference column. By default the result is column-centred (not scaled).
#'
#' @param x Numeric matrix or data frame of strictly positive abundances
#'   (samples in rows, taxa in columns). Zeros must be replaced beforehand,
#'   e.g. with a pseudocount.
#' @param reference Reference taxon, given as a column index or a column name.
#' @param center Logical; centre each column after the transform?
#'
#' @return A numeric matrix with `ncol(x) - 1` columns.
#' @examples
#' x <- matrix(rexp(20), 5, 4, dimnames = list(NULL, letters[1:4]))
#' alr_transform(x, reference = "b")
#' @export
alr_transform <- function(x, reference, center = TRUE) {
  x <- as.matrix(x)
  ref <- resolve_reference(x, reference)
  if (any(x <= 0, na.rm = TRUE)) {
    stop("alr_transform() needs strictly positive abundances; replace zeros first")
  }
  out <- log(x[, -ref, drop = FALSE] / x[, ref])
  if (center) {
    out <- scale(out, scale = FALSE)
    attr(out, "scaled:center") <- NULL
  }
  out
}

#' ALR-transform a reference and a target cohort and align their taxa
#'
#' Applies [alr_transform()] to both cohorts using the same reference taxon,
#' then keeps only the taxa present in both. The column order follows
#' `target_counts`.
#'
#' @param ref_counts Abundance matrix of the reference cohort, which has the
#'   metabolite profiles.
#' @param target_counts Abundance matrix of the target cohort, which has the
#'   outcome.
#' @param reference Reference taxon. A column name is looked up in each matrix
#'   separately. A column index is used as-is in both, so the two matrices must
#'   list their taxa in the same order.
#' @param center Logical; passed to [alr_transform()].
#'
#' @return A list with `X_ref`, `X_target` (aligned ALR matrices), `taxa` (the
#'   shared taxa) and `reference` (the reference taxon's name, or its index if
#'   the columns have no names).
#' @export
prepare_microbiome <- function(ref_counts, target_counts, reference, center = TRUE) {
  ref_counts <- as.matrix(ref_counts)
  target_counts <- as.matrix(target_counts)
  if (is.null(colnames(ref_counts)) || is.null(colnames(target_counts))) {
    stop("both abundance matrices need taxon names as column names")
  }
  ref_idx <- resolve_reference(ref_counts, reference)
  tgt_idx <- resolve_reference(target_counts, reference)
  ref_name <- colnames(ref_counts)[ref_idx]
  if (ref_name != colnames(target_counts)[tgt_idx]) {
    warning("reference column ", reference, " is '", ref_name, "' in ref_counts but '",
            colnames(target_counts)[tgt_idx], "' in target_counts")
  }

  X_ref <- alr_transform(ref_counts, ref_idx, center = center)
  X_target <- alr_transform(target_counts, tgt_idx, center = center)
  taxa <- intersect(colnames(X_target), colnames(X_ref))
  if (length(taxa) == 0) stop("the two cohorts share no taxa")

  list(
    X_ref = X_ref[, match(taxa, colnames(X_ref)), drop = FALSE],
    X_target = X_target[, match(taxa, colnames(X_target)), drop = FALSE],
    taxa = taxa,
    reference = ref_name
  )
}

resolve_reference <- function(x, reference) {
  if (length(reference) != 1) stop("reference must be a single column index or name")
  if (is.character(reference)) {
    idx <- match(reference, colnames(x))
    if (is.na(idx)) stop("reference taxon '", reference, "' not found in column names")
    return(idx)
  }
  idx <- as.integer(reference)
  if (is.na(idx) || idx < 1 || idx > ncol(x)) stop("reference index out of range")
  idx
}
