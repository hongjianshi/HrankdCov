#' Distance covariance of center-outward ranks
#'
#' Computes the test statistic
#' \deqn{n \, \mathrm{dCov}^2_n\Big(\big[J(F_\pm^{(n)}(X_i))\big]_{i=1}^n,
#'   \big[J(G_\pm^{(n)}(Y_i))\big]_{i=1}^n\Big)}
#' of Shi, Drton and Han (2022) and Shi, Hallin, Drton and Han (2022), where
#' \eqn{F_\pm^{(n)}} and \eqn{G_\pm^{(n)}} are the empirical center-outward
#' distribution functions of the two samples (see [co_ranks()]), \eqn{J}
#' is a score function (see [co_scores()]) and \eqn{\mathrm{dCov}^2_n} is
#' the unbiased U-statistic estimator of the squared distance covariance of
#' Székely, Rizzo and Bakirov (2007), computed by U-centering (Székely and
#' Rizzo, 2014).
#'
#' Under independence, and for absolutely continuous distributions, the
#' distribution of the statistic does not depend on the distributions of
#' `x` and `y`.  Use [hrank_dcov_test()] for p-values.
#'
#' @param x,y Numeric vectors, matrices or data frames with the same number
#'   of observations (rows).  Dimensions may differ.
#' @param score One of `"rank"` (default), `"normal"` or `"sign"`, or a score
#'   function; see [co_scores()].  `"rank"` gives the statistic of Shi,
#'   Drton and Han (2022).
#' @param n_R Number of radii of the grids (see [co_grid()]): `NULL` for the
#'   default, a single number used for both samples, or a vector of length 2.
#' @param type Grid type, `"regular"` or `"random"`; see [co_grid()].
#'
#' @return The value of the statistic, a single number.
#'
#' @references
#' Shi, H., Drton, M. and Han, F. (2022). Distribution-free consistent
#' independence tests via center-outward ranks and signs. *Journal of the
#' American Statistical Association*, 117(537), 395--410.
#' \doi{10.1080/01621459.2020.1782223}
#'
#' Shi, H., Hallin, M., Drton, M. and Han, F. (2022). On universally
#' consistent and fully distribution-free rank tests of vector independence.
#' *The Annals of Statistics*, 50(4), 1933--1959. \doi{10.1214/21-AOS2151}
#'
#' Székely, G. J., Rizzo, M. L. and Bakirov, N. K. (2007). Measuring and
#' testing dependence by correlation of distances. *The Annals of
#' Statistics*, 35(6), 2769--2794. \doi{10.1214/009053607000000505}
#'
#' Székely, G. J. and Rizzo, M. L. (2014). Partial distance correlation with
#' methods for dissimilarities. *The Annals of Statistics*, 42(6),
#' 2382--2412. \doi{10.1214/14-AOS1255}
#'
#' @seealso [hrank_dcov_test()]
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400), ncol = 2)
#' y <- x^2 + matrix(rnorm(400), ncol = 2)
#' hrank_dcov(x, y)
#' hrank_dcov(x, y, score = "normal")
hrank_dcov <- function(x, y, score = c("rank", "normal", "sign"), n_R = NULL,
                       type = c("regular", "random")) {
  x <- as_data_matrix(x, "x")
  y <- as_data_matrix(y, "y")
  check_sample_sizes(x, y)
  if (!is.function(score)) score <- match.arg(score)
  type <- match.arg(type)
  n_R <- split_n_R(n_R)

  rx <- co_ranks(x, n_R = n_R[[1]], type = type)
  ry <- co_ranks(y, n_R = n_R[[2]], type = type)
  n_dcov_u(ucenter_dist(co_scores(rx, score)), ucenter_dist(co_scores(ry, score)))
}

split_n_R <- function(n_R) {
  if (is.null(n_R)) return(list(NULL, NULL))
  if (!length(n_R) %in% 1:2) stop("`n_R` must have length 1 or 2.", call. = FALSE)
  list(n_R[[1]], n_R[[length(n_R)]])
}
