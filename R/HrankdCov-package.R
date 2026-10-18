#' HrankdCov: Distribution-free tests of multivariate independence
#'
#' Center-outward ranks and signs (Hallin et al., 2021) extend univariate
#' ranks to random vectors through optimal transport of the data onto a
#' regular grid in the unit ball.  This package computes them and uses them
#' in the distance covariance to obtain tests of independence between two
#' random vectors that are both distribution-free and consistent (Shi,
#' Drton and Han, 2022; Shi, Hallin, Drton and Han, 2022).
#'
#' Main functions:
#' * [hrank_dcov_test()]: the test of independence;
#' * [hrank_dcov()]: the test statistic;
#' * [hrank_dcov_null()]: its distribution-free null distribution, for reuse
#'   across many tests;
#' * [co_ranks()], [co_grid()], [co_scores()]: center-outward ranks, signs,
#'   grids and scores.
#'
#' See `vignette("HrankdCov")` for an introduction.
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
#' @keywords internal
#' @useDynLib HrankdCov, .registration = TRUE
#' @importFrom stats quantile
"_PACKAGE"
