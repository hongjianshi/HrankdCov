#' Center-outward reference grid
#'
#' Builds the grid of `n` points in the open unit ball of \eqn{R^d} onto which
#' the data are optimally transported to obtain center-outward ranks and
#' signs.  The sample size is factorized as \eqn{n = n_R n_S + n_0} with
#' \eqn{n_S} even and \eqn{0 \le n_0 < 2 n_R}; the grid consists of the
#' \eqn{n_R n_S} points
#' \deqn{\frac{r}{n_R + 1} s_k, \qquad r = 1, \dots, n_R,\ k = 1, \dots, n_S,}
#' where \eqn{s_1, \dots, s_{n_S}} are unit vectors, plus \eqn{n_0} extra
#' points close to the origin.
#'
#' The grid is symmetric: whenever \eqn{u} is a grid point, so is \eqn{-u}.
#' To this end the directions consist of \eqn{n_S / 2} unit vectors together
#' with their antipodes, and the \eqn{n_0} extra points are
#' \eqn{\lfloor n_0 / 2 \rfloor} antipodal pairs \eqn{\pm s_k / (2(n_R + 1))}
#' on distinct directions, plus the origin when \eqn{n_0} is odd.  This
#' refines the constructions of Shi, Drton and Han (2022), who use
#' \eqn{n_0} copies of the origin, and of Shi, Hallin, Drton and Han (2022),
#' who place the \eqn{n_0} points on \eqn{n_0} random directions at radius
#' \eqn{1 / (2(n_R + 1))}.  Symmetry ensures that the center-outward signs
#' sum to zero, which is needed for the sign-based statistics of Shi, Drton,
#' Hallin and Han (2025, Remark 2.1).
#'
#' In dimension `d = 1` the grid is
#' \eqn{\{\pm r/(n_R + 1): r = 1, \dots, n_R\}} with \eqn{n_S = 2} and
#' \eqn{n_R = \lfloor n/2 \rfloor}, plus the origin when `n` is odd; `n_R`
#' and `type` are then ignored.
#'
#' The directions are chosen as follows.
#' * `type = "regular"` (default): for `d = 2`, the \eqn{n_S} directions
#'   split the unit circle into arcs of equal length; for `d >= 3`, half of
#'   them are obtained by mapping a Kronecker low-discrepancy sequence
#'   (based on the generalized golden ratio) through the standard normal
#'   quantile function and projecting onto the sphere, and the other half
#'   are their antipodes.  The grid is deterministic, so test statistics
#'   and null distributions are reproducible without a seed.  This
#'   construction differs from the one in Appendix B of Shi, Drton and Han
#'   (2022), which uses a product grid in spherical coordinates, but like it
#'   satisfies the only requirement of the theory, namely that the
#'   empirical distribution of the directions converges to the uniform
#'   distribution on the sphere.
#' * `type = "random"`: half of the directions are drawn independently and
#'   uniformly on the unit sphere (the other half being their antipodes),
#'   and the directions of the extra points are sampled at random.
#'
#' Either choice is independent of the data, which is all that is needed for
#' the tests in this package to be distribution-free.
#'
#' @param n Sample size, an integer of at least 2.
#' @param d Dimension, a positive integer.
#' @param n_R Number of radii (concentric spheres) for `d >= 2`.  Defaults to
#'   \eqn{\lfloor n^{1/4} \rfloor}, the choice made in the simulation code
#'   of Shi, Drton and Han (2022); \eqn{n_S = 2 \lfloor n / (2 n_R) \rfloor}
#'   and \eqn{n_0 = n - n_R n_S} then follow.  The theory only requires
#'   \eqn{n_R, n_S \to \infty}.
#' @param type Either `"regular"` or `"random"`; see Details.
#'
#' @return An `n` by `d` numeric matrix of grid points with attributes `n_R`,
#'   `n_S`, `n_0` and `type`.
#'
#' @references
#' Hallin, M., del Barrio, E., Cuesta-Albertos, J. and Matrán, C. (2021).
#' Distribution and quantile functions, ranks and signs in dimension d: A
#' measure transportation approach. *The Annals of Statistics*, 49(2),
#' 1139--1165. \doi{10.1214/20-AOS1996}
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
#' Shi, H., Drton, M., Hallin, M. and Han, F. (2025). Distribution-free tests
#' of multivariate independence based on center-outward quadrant, Spearman,
#' Kendall, and van der Waerden statistics. *Bernoulli*, 31(1), 106--129.
#' \doi{10.3150/24-BEJ1721}
#'
#' @seealso [co_ranks()]
#' @export
#' @examples
#' g <- co_grid(100, 2)
#' attributes(g)[c("n_R", "n_S", "n_0")]
#' plot(g, asp = 1, pch = 19, cex = 0.5)
co_grid <- function(n, d, n_R = NULL, type = c("regular", "random")) {
  n <- check_count(n, "n", min = 2)
  d <- check_count(d, "d", min = 1)
  type <- match.arg(type)

  if (d == 1) {
    n_R <- n %/% 2
    n_0 <- n %% 2
    r <- seq_len(n_R) / (n_R + 1)
    g <- matrix(c(-rev(r), rep(0, n_0), r), ncol = 1)
    return(new_grid(g, n_R = n_R, n_S = 2L, n_0 = n_0, type = type))
  }

  n_R <- if (is.null(n_R)) max(1L, floor(n^(1 / 4))) else check_count(n_R, "n_R", min = 1)
  n_S <- 2L * (n %/% (2L * n_R))
  n_0 <- n - n_R * n_S
  if (n_S < 2L || n_0 >= n_S) {
    stop(sprintf("`n_R = %d` is too large for `n = %d`: need n_0 < n_S.", n_R, n),
         call. = FALSE)
  }

  half <- if (type == "regular") half_directions(n_S / 2, d) else random_directions(n_S / 2, d)
  s <- rbind(half, -half)
  radii <- seq_len(n_R) / (n_R + 1)
  g <- radii %x% s

  n_pairs <- n_0 %/% 2
  if (n_pairs > 0) {
    idx <- if (type == "regular") {
      floor((seq_len(n_pairs) - 1) * (n_S / 2) / n_pairs) + 1
    } else {
      sample.int(n_S / 2, n_pairs)
    }
    extra <- half[idx, , drop = FALSE] / (2 * (n_R + 1))
    g <- rbind(g, extra, -extra)
  }
  if (n_0 %% 2 == 1) g <- rbind(g, matrix(0, 1, d))
  new_grid(g, n_R = n_R, n_S = n_S, n_0 = n_0, type = type)
}

new_grid <- function(g, n_R, n_S, n_0, type) {
  dimnames(g) <- NULL
  structure(g, n_R = as.integer(n_R), n_S = as.integer(n_S),
            n_0 = as.integer(n_0), type = type)
}

# m deterministic unit vectors in dimension d >= 2 which, together with their
# antipodes, are roughly uniform on the sphere.
half_directions <- function(m, d) {
  if (d == 2) {
    angle <- pi * (seq_len(m) - 1) / m
    return(cbind(cos(angle), sin(angle)))
  }
  # Kronecker sequence with the generalized golden ratio phi_d, the positive
  # root of x^(d + 1) = x + 1.
  phi <- 2
  for (i in 1:100) phi <- (1 + phi)^(1 / (d + 1))
  alpha <- (1 / phi)^seq_len(d)
  u <- (0.5 + outer(seq_len(m), alpha)) %% 1
  normalize_rows(stats::qnorm(u))
}

random_directions <- function(m, d) {
  normalize_rows(matrix(stats::rnorm(m * d), m, d))
}

normalize_rows <- function(z) z / sqrt(rowSums(z^2))
