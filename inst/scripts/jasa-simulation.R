# Simulation study of Section 6 of
#
#   Shi, H., Drton, M. and Han, F. (2022). Distribution-free consistent
#   independence tests via center-outward ranks and signs. Journal of the
#   American Statistical Association, 117(537), 395-410.
#
# Compares the empirical rejection rates at level alpha = 0.05 of
#   (i)   the center-outward rank distance covariance test (HrankdCov),
#   (ii)  distance covariance of coordinatewise ranks (Lin, 2017), and
#   (iii) distance covariance (Szekely and Rizzo, 2013),
# for Gaussian data (Example 6.1) and data with Cauchy margins (Example 6.2).
# Tests (ii) and (iii) use n permutations, as in the paper.
#
# Requires the packages mvtnorm and energy.  The full design takes a long
# time to run; reduce `n_rep` or the parameter grids for a quick look.

library(HrankdCov)
library(mvtnorm)
library(energy)

alpha <- 0.05
n_rep <- 1000
n_list <- c(216, 432, 864, 1728)
p_list <- c(2, 3, 5, 7)
tau_list <- c(0, 0.5, 0.9)
rho_list <- seq(0, 0.15, by = 0.005)
tests <- c("HrankdCov", "rank dCov", "dCov")
examples <- c("Gaussian", "Cauchy")

rejections <- array(
  NA_real_,
  dim = c(length(n_list), length(p_list), length(tau_list), length(rho_list),
          length(examples), length(tests)),
  dimnames = list(n = n_list, p = p_list, tau = tau_list, rho = rho_list,
                  example = examples, test = tests)
)

for (n in n_list) {
  for (p in p_list) {
    # The null distribution is the same for all data-generating distributions,
    # so it is computed once for each (n, p).
    null <- hrank_dcov_null(n, p, p, method = "asymptotic")
    for (tau in tau_list) {
      for (rho in rho_list) {
        sigma <- diag(2 * p)
        sigma[1, 2] <- sigma[2, 1] <- tau
        sigma[1, p + 1] <- sigma[p + 1, 1] <- rho

        reject <- replicate(n_rep, {
          z <- rmvnorm(n, sigma = sigma)
          sapply(list(Gaussian = z, Cauchy = qt(pnorm(z), df = 1)), function(w) {
            x <- w[, 1:p]
            y <- w[, (p + 1):(2 * p)]
            c(hrank_dcov_test(x, y, null = null)$p.value,
              dcov.test(apply(x, 2, rank), apply(y, 2, rank), R = n)$p.value,
              dcov.test(x, y, R = n)$p.value) <= alpha
          })
        })
        # `reject` is a 3 (tests) x 2 (examples) x n_rep array.
        rejections[as.character(n), as.character(p), as.character(tau),
                   as.character(rho), , ] <- t(apply(reject, 1:2, mean))
      }
    }
  }
}

saveRDS(rejections, "jasa-simulation.rds")

# Empirical sizes (rho = 0), as in Table 6.1 of the paper
round(rejections[, , , "0", , ], 3)
