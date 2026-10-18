`%||%` <- function(x, y) if (is.null(x)) y else x

as_data_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("All columns of `%s` must be numeric.", name), call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (length(dim(x)) != 2) stop(sprintf("`%s` must be a vector or a matrix.", name), call. = FALSE)
  if (!all(is.finite(x))) {
    stop(sprintf("`%s` must not contain missing or infinite values.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  dimnames(x) <- NULL
  x
}

check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < min) {
    stop(sprintf("`%s` must be an integer of at least %d.", name, min), call. = FALSE)
  }
  as.integer(x)
}

check_grid <- function(grid, n, d) {
  g <- as_data_matrix(grid, "grid")
  if (nrow(g) != n || ncol(g) != d) {
    stop(sprintf("`grid` must be a %d x %d matrix.", n, d), call. = FALSE)
  }
  if (any(rowSums(g^2) >= 1)) stop("Grid points must lie in the open unit ball.", call. = FALSE)
  for (a in c("n_R", "n_S", "n_0", "type")) attr(g, a) <- attr(grid, a)
  g
}

check_sample_sizes <- function(x, y) {
  n <- nrow(x)
  if (nrow(y) != n) stop("`x` and `y` must have the same number of observations.", call. = FALSE)
  if (n < 4) stop("At least 4 observations are needed.", call. = FALSE)
  n
}

# U-centered distance matrix (Szekely and Rizzo, 2014).
ucenter_dist <- function(z) {
  D <- as.matrix(stats::dist(z))
  n <- nrow(D)
  r <- rowSums(D)
  A <- D - outer(r, r, "+") / (n - 2) + sum(r) / ((n - 1) * (n - 2))
  diag(A) <- 0
  A
}

# n times the unbiased (U-statistic) estimator of the squared distance
# covariance, from U-centered distance matrices.
n_dcov_u <- function(A, B) sum(A * B) / (nrow(A) - 3)
