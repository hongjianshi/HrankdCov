#' Distribution-free test of independence based on center-outward ranks
#'
#' Tests the null hypothesis that two random vectors `X` and `Y` are
#' independent, using the distance covariance of their (scored)
#' center-outward ranks; see [hrank_dcov()].  The test is
#' * distribution-free: for absolutely continuous `X` and `Y`, the null
#'   distribution of the statistic does not depend on their distributions,
#'   so that critical values are universal;
#' * consistent: it has power tending to one against every fixed
#'   alternative, for distributions with nonvanishing Lebesgue densities
#'   (Shi, Drton and Han, 2022, Theorem 3.2) and, more generally, for the
#'   class of absolutely continuous distributions for which the empirical
#'   center-outward distribution function is uniformly consistent (Shi,
#'   Hallin, Drton and Han, 2022, Theorem 4.1 and Proposition 4.2).
#'
#' @inheritParams hrank_dcov
#' @inheritParams hrank_dcov_null
#' @param null Optionally, a null distribution computed by
#'   [hrank_dcov_null()] for the same sample size and dimensions.  It is
#'   then used for the p-value, and `score`, `method`, `B`, `n_R` and
#'   `type` are taken from it.  This saves time when many tests are run with
#'   the same `n`, `p` and `q`, as in simulation studies.
#'
#' @return A list of class `"htest"` with components `statistic` (the value
#'   of \eqn{n\,\mathrm{dCov}^2_n}), `estimate` (\eqn{\mathrm{dCov}^2_n}),
#'   `p.value`, `method`, `data.name`, and `null`, the
#'   `"hrank_dcov_null"` object used, which can be reused through the `null`
#'   argument.
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
#' @seealso [hrank_dcov()], [hrank_dcov_null()], [co_ranks()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- matrix(rnorm(2 * n), ncol = 2)
#'
#' # Nonlinear, non-monotone dependence
#' y <- cbind(x[, 1]^2, x[, 2]) + matrix(rnorm(2 * n), ncol = 2)
#' hrank_dcov_test(x, y, B = 199)
#'
#' # Independence, heavy tails: the test is valid regardless
#' y0 <- matrix(rt(3 * n, df = 1), ncol = 3)
#' hrank_dcov_test(x, y0, score = "normal", method = "asymptotic")
hrank_dcov_test <- function(x, y, score = c("rank", "normal", "sign"),
                            method = c("permutation", "asymptotic"), B = 999,
                            n_R = NULL, type = c("regular", "random"),
                            null = NULL) {
  data_name <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  x <- as_data_matrix(x, "x")
  y <- as_data_matrix(y, "y")
  n <- check_sample_sizes(x, y)

  if (is.null(null)) {
    null <- hrank_dcov_null(n, ncol(x), ncol(y), score = score, method = method,
                            B = B, n_R = n_R, type = type)
  } else {
    if (!inherits(null, "hrank_dcov_null")) {
      stop("`null` must be an object returned by hrank_dcov_null().", call. = FALSE)
    }
    if (null$n != n || null$p != ncol(x) || null$q != ncol(y)) {
      stop(sprintf("`null` was computed for n = %d, p = %d, q = %d, but the data have n = %d, p = %d, q = %d.",
                   null$n, null$p, null$q, n, ncol(x), ncol(y)), call. = FALSE)
    }
  }

  rx <- co_ranks(x, grid = null$grid_x)
  ry <- co_ranks(y, grid = null$grid_y)
  stat <- n_dcov_u(ucenter_dist(co_scores(rx, null$score)),
                   ucenter_dist(co_scores(ry, null$score)))

  structure(
    list(
      statistic = c("n dCov^2" = stat),
      estimate = c("dCov^2" = stat / n),
      p.value = null_pvalue(null, stat),
      method = sprintf("Center-outward rank distance covariance test of independence (%s; %s)",
                       score_label(null$score), null_label(null)),
      data.name = data_name,
      null = null
    ),
    class = "htest"
  )
}
