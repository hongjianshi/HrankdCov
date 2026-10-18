#' Center-outward ranks and signs
#'
#' Computes the empirical center-outward distribution function of Hallin,
#' del Barrio, Cuesta-Albertos and Matrán (2021): the sample is matched to the
#' grid points of [co_grid()] by the assignment that minimizes the total
#' squared Euclidean distance,
#' \deqn{\hat\sigma = \arg\min_{\sigma} \sum_{i=1}^n \| X_i - g_{\sigma(i)} \|^2,}
#' and \eqn{F_\pm^{(n)}(X_i) = g_{\hat\sigma(i)}}.  The center-outward rank of
#' \eqn{X_i} is \eqn{(n_R + 1)\|F_\pm^{(n)}(X_i)\|} and its center-outward sign
#' is \eqn{F_\pm^{(n)}(X_i) / \|F_\pm^{(n)}(X_i)\|}.
#'
#' The assignment problem is solved exactly with a shortest augmenting path
#' (Hungarian) algorithm implemented in C, in \eqn{O(n^3)} time and
#' \eqn{O(n^2)} memory; for `d = 1` it reduces to sorting.  Because the
#' optimal assignment is unchanged by translations and by positive rescaling
#' of the data, the data are centered and scaled before the assignment for
#' numerical stability.
#'
#' Ties among the observations (which have probability zero for absolutely
#' continuous distributions) are broken deterministically: by the order of
#' the observations in dimension one, and by the assignment algorithm
#' otherwise.
#'
#' @param x A numeric vector, matrix or data frame with observations in rows.
#' @param n_R,type Passed to [co_grid()].
#' @param grid Optionally, a numeric matrix with `nrow(x)` rows and
#'   `ncol(x)` columns of grid points in the open unit ball, used instead of
#'   [co_grid()].  Attributes `n_R`, `n_S` and `n_0` are used when present.
#'
#' @return An object of class `"co_ranks"`, a list with components
#' \describe{
#'   \item{`F`}{the `n` by `d` matrix of grid points assigned to the
#'     observations, i.e. the empirical center-outward distribution function
#'     evaluated at the data;}
#'   \item{`rank`}{the center-outward ranks \eqn{(n_R + 1)\|F\|}: integers in
#'     \eqn{1, \dots, n_R}, plus 0 for the origin and 1/2 for the extra
#'     points near the origin when \eqn{n_0 > 1};}
#'   \item{`sign`}{the `n` by `d` matrix of center-outward signs (unit
#'     vectors, or zero for the origin);}
#'   \item{`index`}{the row of `grid` assigned to each observation;}
#'   \item{`grid`, `n_R`, `n_S`, `n_0`}{the grid and its factorization;}
#'   \item{`data`}{the data as a matrix.}
#' }
#'
#' @references
#' Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021).
#' Distribution and quantile functions, ranks and signs in dimension d: A
#' measure transportation approach. *The Annals of Statistics*, 49(2),
#' 1139--1165. \doi{10.1214/20-AOS1996}
#'
#' @seealso [co_grid()], [co_scores()]
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400), ncol = 2)
#' r <- co_ranks(x)
#' r
#' table(r$rank)
#' plot(r)
#'
#' # In dimension one, center-outward ranks are a monotone transformation of
#' # the usual ranks.
#' z <- rnorm(11)
#' all.equal(order(co_ranks(z)$F), order(z))
co_ranks <- function(x, n_R = NULL, type = c("regular", "random"), grid = NULL) {
  x <- as_data_matrix(x, "x")
  n <- nrow(x)
  d <- ncol(x)
  g <- if (is.null(grid)) co_grid(n, d, n_R = n_R, type = match.arg(type)) else check_grid(grid, n, d)

  index <- co_assign(x, g)
  F <- g[index, , drop = FALSE]
  norm <- sqrt(rowSums(F^2))
  sign <- F / ifelse(norm > 0, norm, 1)

  n_R <- attr(g, "n_R")
  rank <- if (is.null(n_R)) rep(NA_real_, n) else round((n_R + 1) * norm, 10)

  structure(
    list(F = F, rank = rank, sign = sign, index = index, grid = g,
         n_R = n_R, n_S = attr(g, "n_S"), n_0 = attr(g, "n_0"), data = x),
    class = "co_ranks"
  )
}

# Row i of x is matched to row index[i] of the grid g.
co_assign <- function(x, g) {
  n <- nrow(x)
  if (ncol(x) == 1) {
    index <- integer(n)
    index[order(x[, 1])] <- order(g[, 1])
    return(index)
  }
  # Minimizing sum ||x_i - g_sigma(i)||^2 is equivalent to maximizing
  # sum <x_i, g_sigma(i)>, which is invariant to centering and scaling of x.
  x <- sweep(x, 2, colMeans(x))
  s <- sqrt(mean(x^2))
  if (s > 0) x <- x / s
  .Call(C_lsap, -tcrossprod(x, g))
}

#' @export
print.co_ranks <- function(x, ...) {
  n <- nrow(x$F)
  d <- ncol(x$F)
  cat(sprintf("Center-outward ranks of %d observations in dimension %d\n", n, d))
  if (!is.null(x$n_R)) {
    cat(sprintf("Grid: n_R = %d radii x n_S = %d directions + n_0 = %d (%s)\n",
                x$n_R, x$n_S, x$n_0, attr(x$grid, "type") %||% "user-supplied"))
    cat("Rank counts:\n")
    print(table(x$rank, dnn = NULL))
  } else {
    cat("Grid: user-supplied\n")
  }
  invisible(x)
}

#' Plot center-outward ranks
#'
#' For bivariate data, plots the observations colored by their center-outward
#' rank, so that the empirical center-outward quantile contours become
#' visible.  For univariate data, plots the empirical center-outward
#' distribution function.
#'
#' @param x An object returned by [co_ranks()].
#' @param col Colors for the ranks \eqn{0, \dots, n_R} (bivariate case).
#' @param pch,... Further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.co_ranks <- function(x, col = NULL, pch = 19, ...) {
  d <- ncol(x$data)
  if (d == 1) {
    graphics::plot(x$data[, 1], x$F[, 1], pch = pch,
                   xlab = "x", ylab = "center-outward distribution function", ...)
  } else if (d == 2) {
    level <- ceiling(x$rank)
    if (anyNA(level)) level <- rep(1, length(level))
    levels <- sort(unique(level))
    if (is.null(col)) col <- grDevices::hcl.colors(length(levels), "viridis", rev = TRUE)
    graphics::plot(x$data, col = col[match(level, levels)], pch = pch,
                   xlab = "x1", ylab = "x2", asp = 1, ...)
  } else {
    stop("Plotting is only available in dimensions 1 and 2.", call. = FALSE)
  }
  invisible(x)
}

#' Scored center-outward ranks
#'
#' Applies a score function to center-outward ranks, giving the vectors
#' \deqn{J\big(\|F_\pm^{(n)}(X_i)\|\big) \, S_\pm^{(n)}(X_i),}
#' where \eqn{S_\pm^{(n)}} are the center-outward signs (Shi, Hallin, Drton
#' and Han, 2022).
#'
#' The built-in score functions are
#' * `"rank"` (Wilcoxon or Spearman scores): \eqn{J(u) = u}, so that the
#'   scored ranks are the grid points themselves;
#' * `"normal"` (van der Waerden scores):
#'   \eqn{J(u) = \sqrt{F^{-1}_{\chi^2_d}(u)}}, which maps the uniform
#'   distribution on the ball to the standard normal distribution on
#'   \eqn{R^d};
#' * `"sign"`: \eqn{J(u) = 1}, so that only the center-outward signs are
#'   used.
#'
#' Alternatively, `score` can be a vectorized function \eqn{J} from
#' \eqn{[0, 1)} to \eqn{[0, \infty)}.  The origin is always mapped to the
#' origin.
#'
#' @param x An object returned by [co_ranks()], or a numeric matrix of
#'   points in the open unit ball (for example, a grid from [co_grid()]).
#' @param score One of `"rank"`, `"normal"` or `"sign"`, or a score function;
#'   see Details.
#'
#' @return A numeric matrix with the same dimensions as `x$F` (or `x`).
#'
#' @references
#' Shi, H., Hallin, M., Drton, M. and Han, F. (2022). On universally
#' consistent and fully distribution-free rank tests of vector independence.
#' *The Annals of Statistics*, 50(4), 1933--1959. \doi{10.1214/21-AOS2151}
#'
#' @export
#' @examples
#' set.seed(1)
#' r <- co_ranks(matrix(rnorm(600), ncol = 3))
#' z <- co_scores(r, "normal")
#' colMeans(z)
#' round(cov(z), 2)
co_scores <- function(x, score = c("rank", "normal", "sign")) {
  F <- if (inherits(x, "co_ranks")) x$F else as_data_matrix(x, "x")
  J <- score_function(score, ncol(F))
  norm <- sqrt(rowSums(F^2))
  if (any(norm >= 1)) stop("Points must lie in the open unit ball.", call. = FALSE)
  pos <- norm > 0
  out <- matrix(0, nrow(F), ncol(F))
  if (any(pos)) {
    j <- J(norm[pos])
    if (length(j) != sum(pos) || anyNA(j) || any(j < 0)) {
      stop("The score function must return nonnegative values.", call. = FALSE)
    }
    out[pos, ] <- F[pos, , drop = FALSE] * (j / norm[pos])
  }
  out
}

score_function <- function(score, d) {
  if (is.function(score)) return(score)
  switch(match.arg(score, c("rank", "normal", "sign")),
    rank = function(u) u,
    normal = function(u) sqrt(stats::qchisq(u, df = d)),
    sign = function(u) rep(1, length(u))
  )
}

score_label <- function(score) {
  if (is.function(score)) return("user-supplied scores")
  switch(match.arg(score, c("rank", "normal", "sign")),
    rank = "rank (Wilcoxon) scores",
    normal = "normal (van der Waerden) scores",
    sign = "sign scores"
  )
}
