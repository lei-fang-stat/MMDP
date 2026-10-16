#' High-dimensional permutation test for one covariate (FLhd)
#'
#' Tests whether the first column of `X` is associated with `y` after
#' adjusting for the remaining columns, using the residual-permutation
#' procedure of Hemerik, Thoresen and Finos (2021). The nuisance columns are
#' handled by a ridge fit. Its penalty comes either from cross-validation, as in
#' `phd::FLhd()`, or from REML: the ridge is refitted as a random-effects model
#' in 'glmmTMB' and the penalty is taken as `sigma^2 / tau^2`.
#'
#' With a fixed seed and the same `lambda` setting, results match
#' `phd::FLhd()` (`lambda = "lambda.min"` / `"lambda.1se"`) and the original
#' `FLHD2()` script (`lambda = "reml"`, one nuisance covariate) exactly.
#'
#' @param y Numeric response vector (or one-column matrix).
#' @param X Numeric matrix. Column 1 is the covariate of interest. The other
#'   columns are nuisance covariates.
#' @param X1 Optional covariate of interest. If supplied it is placed in front
#'   of `X`, so every column of `X` is treated as a nuisance covariate.
#' @param nperm Number of permutations, including the identity.
#' @param lambda Ridge penalty for the nuisance fit. Use `"lambda.min"` or
#'   `"lambda.1se"` for cross-validation (needs at least two nuisance
#'   covariates), `"reml"` for REML, or a positive number. The default,
#'   `"auto"`, uses REML when there is one nuisance covariate and
#'   `"lambda.min"` otherwise. This is the rule used in the real-data
#'   analysis.
#' @param flip Logical; use sign-flipping instead of permutation.
#' @param nfolds Number of cross-validation folds (CV penalties only).
#' @param statistic `"partialcor"` (default) or `"semipartialcor"`.
#'
#' @return A two-sided p-value. The attributes `lambda` (penalty used) and
#'   `lambda_method` record how the penalty was chosen.
#' @references Hemerik, J., Thoresen, M. and Finos, L. (2021). Permutation
#'   testing in high-dimensional linear models: an empirical investigation.
#'   *Journal of Statistical Computation and Simulation*, 91(5), 897-914.
#' @export
flhd_test <- function(y, X, X1, nperm = 20000, lambda = "auto", flip = FALSE,
                      nfolds = 10, statistic = c("partialcor", "semipartialcor")) {
  statistic <- match.arg(statistic)
  flip <- as.logical(flip)
  if (length(flip) != 1 || is.na(flip)) stop("argument flip should be TRUE or FALSE")
  X <- as.matrix(X)
  if (!missing(X1)) {
    if (length(X1) != nrow(X)) stop("X and X1 should both have n rows")
    X <- cbind(X1, X)
  }
  y <- as.numeric(y)
  if (length(y) != nrow(X)) stop("the length of y should equal the nr of rows of X")
  if (is.numeric(lambda) && lambda <= 0) stop("lambda should be larger than 0")
  if (nperm < 100) message("taking the nr of perms very small leads to inaccurate p-values")

  w <- nperm
  n <- nrow(X)
  q <- ncol(X) - 1
  if (q < 1) stop("X needs the covariate of interest plus at least one nuisance covariate")

  X <- apply(X, 2, function(x) (x - mean(x)) / sd(x))
  y <- (y - mean(y)) / sd(y)
  X0 <- X[, -1, drop = FALSE]

  if (identical(lambda, "auto")) lambda <- if (q == 1) "reml" else "lambda.min"
  if (is.numeric(lambda)) {
    lam <- lambda
    lambda_method <- "fixed"
  } else if (identical(lambda, "reml")) {
    lam <- reml_ridge_lambda(y, X0)
    lambda_method <- "reml"
  } else if (lambda %in% c("lambda.min", "lambda.1se")) {
    if (q < 2) stop("cross-validated ridge needs at least 2 nuisance covariates; use lambda = \"reml\"")
    range <- 10^seq(5, -5, -0.1)
    modelzy <- glmnet::cv.glmnet(x = X0, y = y, family = "gaussian", alpha = 0,
                                 intercept = FALSE, lambda = range, nfolds = nfolds)
    # phd::FLhd() also cross-validates a ridge fit of X[, 1] but then uses the
    # y-based penalty for both projections. The fit is kept only so the RNG
    # stream, and therefore the p-value, matches phd::FLhd().
    glmnet::cv.glmnet(x = X0, y = X[, 1], family = "gaussian", alpha = 0,
                      intercept = FALSE, lambda = range, nfolds = nfolds)
    lam <- n * modelzy[[lambda]]
    lambda_method <- paste0("cv_", lambda)
  } else {
    stop("lambda should be \"auto\", \"reml\", \"lambda.min\", \"lambda.1se\" or a positive number")
  }

  Xhatmat <- X0 %*% solve(t(X0) %*% X0 + lam * diag(q)) %*% t(X0)
  yhat <- Xhatmat %*% y
  resx <- X[, 1] - Xhatmat %*% X[, 1]
  resids <- y - yhat
  target <- if (statistic == "partialcor") resx else X[, 1]

  T <- numeric(w)
  T[1] <- cor(resids, target)
  for (j in 2:w) {
    Ystar <- yhat + sample(resids, size = n, replace = FALSE)
    # Unused draw, kept so the RNG stream matches phd::FLhd()
    sample(resids, size = n, replace = FALSE)
    if (flip) {
      flips <- 2 * stats::rbinom(n, 1, 0.5) - 1
      Ystar <- yhat + flips * resids
    }
    T[j] <- cor(Ystar - Xhatmat %*% Ystar, target)
  }
  p <- 2 * min(sum(T >= T[1]) / w, sum(T <= T[1]) / w)
  structure(p, lambda = lam, lambda_method = lambda_method)
}

# Ridge penalty sigma^2 / tau^2 from a REML fit of y ~ 0 + (0 + X0 | grp).
# This matches FLHD2() for one covariate. With several covariates a homdiag()
# term gives them one shared variance, which is a single ridge penalty.
reml_ridge_lambda <- function(y, X0) {
  q <- ncol(X0)
  df <- data.frame(y = y, X0, grp = factor(rep(1, length(y))))
  xnames <- paste0("x", seq_len(q))
  names(df) <- c("y", xnames, "grp")
  re <- paste0("(0 + ", paste(xnames, collapse = " + "), " | grp)")
  if (q > 1) re <- paste0("homdiag", re)
  form <- stats::as.formula(paste("y ~ 0 +", re))
  fit <- glmmTMB::glmmTMB(form, data = df, REML = TRUE)
  sigma_sq <- stats::sigma(fit)^2
  tau_sq <- diag(as.matrix(glmmTMB::VarCorr(fit)$cond$grp))[1]
  unname(sigma_sq / tau_sq)
}

#' Permutation test of a marginal correlation
#'
#' Two-sided permutation test of `cor(x, y)`, used when the outcome lasso
#' selects no microbiome taxa.
#'
#' @param x,y Numeric vectors of equal length.
#' @param nperm Number of permutations.
#' @return The share of permuted `|cor|` values that are at least as large as the
#'   observed `|cor|`.
#' @export
perm_cor_test <- function(x, y, nperm = 5000) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("x and y must have the same length")
  obs <- cor(x, y)
  perm <- vapply(seq_len(nperm), function(k) cor(x, sample(y)), numeric(1))
  mean(abs(perm) >= abs(obs))
}
