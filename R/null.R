#' Null distribution of the center-outward rank distance covariance
#'
#' Computes the null distribution of the statistic of [hrank_dcov()].  Under
#' independence of two absolutely continuous random vectors, the
#' center-outward ranks of the two samples are independent and each is a
#' uniformly distributed permutation of its grid.  Hence the null
#' distribution of the statistic depends only on `n`, the dimensions, the
#' grids and the score function, and not on the distributions of the data.
#' It can be computed once and reused for any number of tests.
#'
#' Two methods are available.
#' * `"permutation"`: `B` independent uniform random permutations of the
#'   second grid are drawn and the statistic is evaluated on the (scored)
#'   grids.  The resulting Monte Carlo p-value
#'   \eqn{(1 + \#\{b: T_b \ge T\}) / (B + 1)} is exact, i.e. the test has
#'   level at most \eqn{\alpha} for every sample size.
#' * `"asymptotic"`: the limiting null distribution
#'   \eqn{\sum_{k, l} \lambda_k \mu_l (\xi_{kl}^2 - 1)}, with
#'   \eqn{\xi_{kl}} independent standard normal (Shi, Drton and Han, 2022,
#'   Theorem 3.1; Shi, Hallin, Drton and Han, 2022, Corollary 5.1).  The
#'   weights \eqn{\lambda_k} (resp. \eqn{\mu_l}) are the eigenvalues of the
#'   doubly centered kernel \eqn{\|J(R)S - J(R')S'\|} on the spherical
#'   uniform distribution in dimension `p` (resp. `q`), i.e. the limit of
#'   the (scored) grids.  They are computed to high accuracy by exploiting
#'   rotational invariance: by the Funk--Hecke formula the kernel
#'   decomposes over spherical harmonics, so that the eigenvalues of each
#'   harmonic degree are those of a one-dimensional radial integral
#'   equation, which is solved by Gauss--Legendre quadrature.  This refines
#'   the grid-based approximation of Shi, Drton and Han (2022, Section 5.2)
#'   and reproduces their tabulated critical values; for `p = 1` the
#'   eigenvalues are \eqn{4 / (\pi^2 k^2)}.  The largest 1000 products
#'   \eqn{\lambda_k \mu_l} are kept, the remaining ones are replaced by a
#'   single scaled chi-square term with the exact mean and variance, and
#'   tail probabilities are computed by numerical inversion of the
#'   characteristic function (Imhof, 1961).  No resampling is involved.
#'   The eigenvalues for the built-in scores and dimensions up to 10 are
#'   precomputed and stored in the package; other cases take about a
#'   second per dimension.  For moderate `n` the resulting test tends to be
#'   slightly conservative.
#'
#' @param n Sample size, an integer of at least 4.
#' @param p,q Dimensions of the two samples.
#' @param method `"permutation"` (default) or `"asymptotic"`; see Details.
#' @param B Number of Monte Carlo permutations for `method = "permutation"`.
#' @inheritParams hrank_dcov
#'
#' @return An object of class `"hrank_dcov_null"`, which can be passed to
#'   [hrank_dcov_test()] through its `null` argument and has `print()` and
#'   `quantile()` methods.  For `method = "permutation"`, the component
#'   `statistics` holds the `B` simulated values of the statistic.
#'
#' @references
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in
#' normal variables. *Biometrika*, 48(3/4), 419--426.
#' \doi{10.1093/biomet/48.3-4.419}
#'
#' Shi, H., Drton, M. and Han, F. (2022). Distribution-free consistent
#' independence tests via center-outward ranks and signs. *Journal of the
#' American Statistical Association*, 117(537), 395--410.
#' \doi{10.1080/01621459.2020.1782223}
#'
#' Shi, H., Hallin, M., Drton, M. and Han, F. (2022). On universally
#' consistent and fully distribution-free rank tests of vector independence.
#' *The Annals of Statistics*, 50(4), 1933--1959. \doi{10.1214/21-AOS2151}
#'
#' @seealso [hrank_dcov_test()]
#' @export
#' @examples
#' set.seed(1)
#' null <- hrank_dcov_null(200, 2, 3, B = 499)
#' null
#' quantile(null, 0.95)
#'
#' # Critical values from the limiting distribution
#' quantile(hrank_dcov_null(200, 2, 3, method = "asymptotic"), 0.95)
hrank_dcov_null <- function(n, p, q, score = c("rank", "normal", "sign"),
                            method = c("permutation", "asymptotic"), B = 999,
                            n_R = NULL, type = c("regular", "random")) {
  n <- check_count(n, "n", min = 4)
  p <- check_count(p, "p")
  q <- check_count(q, "q")
  if (!is.function(score)) score <- match.arg(score)
  method <- match.arg(method)
  type <- match.arg(type)
  n_R <- split_n_R(n_R)

  grid_x <- co_grid(n, p, n_R = n_R[[1]], type = type)
  grid_y <- co_grid(n, q, n_R = n_R[[2]], type = type)

  out <- list(n = n, p = p, q = q, score = score, method = method, type = type,
              grid_x = grid_x, grid_y = grid_y)
  if (method == "permutation") {
    B <- check_count(B, "B")
    sx <- co_scores(grid_x, score)
    sy <- co_scores(grid_y, score)
    perms <- vapply(seq_len(B), function(b) sample.int(n), integer(n))
    perms <- matrix(perms, nrow = n)
    out$B <- B
    out$statistics <- .Call(C_perm_dcov, ucenter_dist(sx), ucenter_dist(sy), perms) / (n - 3)
  } else {
    out$weights <- chisq_weights(limit_eigenvalues_for(score, p),
                                 limit_eigenvalues_for(score, q))
  }
  structure(out, class = "hrank_dcov_null")
}

#' @export
print.hrank_dcov_null <- function(x, ...) {
  cat("Null distribution of the center-outward rank distance covariance\n")
  cat(sprintf("n = %d, p = %d, q = %d, %s, %s grid\n",
              x$n, x$p, x$q, score_label(x$score), x$type))
  cat(sprintf("Method: %s\n", null_label(x)))
  cat("Upper quantiles of n dCov^2:\n")
  print(signif(stats::quantile(x, c(0.9, 0.95, 0.99)), 4))
  invisible(x)
}

#' Quantiles of the null distribution
#'
#' Critical values of the center-outward rank distance covariance test.
#'
#' @param x An object returned by [hrank_dcov_null()].
#' @param probs Numeric vector of probabilities.
#' @param ... Passed to [stats::quantile()] for permutation null
#'   distributions.
#' @return A named numeric vector of quantiles.
#' @export
quantile.hrank_dcov_null <- function(x, probs = c(0.9, 0.95, 0.99), ...) {
  if (x$method == "permutation") return(stats::quantile(x$statistics, probs, ...))
  w <- x$weights
  sd <- sqrt(2 * w$sum_w2)
  q <- vapply(probs, function(prob) {
    if (prob <= 0) return(-w$sum_w)
    if (prob >= 1) return(Inf)
    stats::uniroot(function(t) wchisq_upper(t, w) - (1 - prob),
                   lower = -w$sum_w, upper = 5 * sd, extendInt = "downX",
                   tol = 1e-8 * sd)$root
  }, numeric(1))
  names(q) <- paste0(format(100 * probs, trim = TRUE), "%")
  q
}

null_label <- function(null) {
  if (null$method == "permutation") {
    sprintf("permutation (Monte Carlo) null distribution with %d replicates", null$B)
  } else {
    "asymptotic (weighted chi-square) null distribution"
  }
}

null_pvalue <- function(null, stat) {
  if (null$method == "permutation") {
    s <- null$statistics
    (1 + sum(s >= stat - 1e-10 * max(1, abs(stat)))) / (length(s) + 1)
  } else {
    wchisq_upper(stat, null$weights)
  }
}

# Weights of the limiting distribution sum_{k,l} lambda_k mu_l (xi_kl^2 - 1)
# from the eigenvalue lists of the two dimensions (see limit_eigen()).  The
# `keep` largest products are kept; the remaining ones are summarized by a
# scaled chi-square variable with matching mean and variance, using the
# exact sums sum(lambda) sum(mu) and sum(lambda^2) sum(mu^2).
chisq_weights <- function(ex, ey, keep = 1000) {
  w <- sort(outer(ex$values, ey$values), decreasing = TRUE)
  w <- w[seq_len(min(keep, length(w)))]
  sum_w <- ex$sum_w * ey$sum_w
  sum_w2 <- ex$sum_w2 * ey$sum_w2
  df <- rep(1, length(w))
  rest <- sum_w - sum(w)
  rest2 <- sum_w2 - sum(w^2)
  if (rest > 1e-12 * sum_w && rest2 > 0) {
    w <- c(w, rest2 / rest)
    df <- c(df, rest^2 / rest2)
  }
  list(w = w, df = df, sum_w = sum_w, sum_w2 = sum_w2)
}

# P(sum_k w_k (chi^2_{df_k} - df_k) > t) by Imhof's (1961) method.
wchisq_upper <- function(t, weights) {
  w <- weights$w
  h <- weights$df
  x <- t + sum(w * h)
  integrand <- function(u) {
    uw <- outer(u, w)
    theta <- 0.5 * drop(atan(uw) %*% h) - 0.5 * x * u
    log_rho <- 0.25 * drop(log1p(uw^2) %*% h)
    sin(theta) / (u * exp(log_rho))
  }
  res <- stats::integrate(integrand, 0, Inf, subdivisions = 10000L,
                          rel.tol = 1e-10, abs.tol = 1e-12, stop.on.error = FALSE)
  min(1, max(0, 0.5 + res$value / pi))
}
