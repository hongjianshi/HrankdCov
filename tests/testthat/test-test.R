test_that("the test returns an htest object", {
  set.seed(1)
  x <- matrix(rnorm(100), ncol = 2)
  y <- matrix(rnorm(150), ncol = 3)
  res <- hrank_dcov_test(x, y, B = 99)
  expect_s3_class(res, "htest")
  expect_named(res$statistic, "n dCov^2")
  expect_equal(unname(res$statistic), hrank_dcov(x, y))
  expect_true(res$p.value > 0 && res$p.value <= 1)
  expect_equal(res$data.name, "x and y")
  expect_s3_class(res$null, "hrank_dcov_null")
  expect_output(print(res), "Center-outward rank distance covariance test")
})

test_that("strong dependence is detected", {
  set.seed(2)
  x <- matrix(rnorm(400), ncol = 2)
  y <- cbind(x[, 1]^2, sin(3 * x[, 2])) + 0.3 * matrix(rnorm(400), ncol = 2)
  expect_lt(hrank_dcov_test(x, y, B = 199)$p.value, 0.01)
  expect_lt(hrank_dcov_test(x, y, method = "asymptotic")$p.value, 0.001)
  expect_lt(hrank_dcov_test(x, y, score = "normal", B = 199)$p.value, 0.01)
})

test_that("a precomputed null distribution can be reused", {
  set.seed(3)
  null <- hrank_dcov_null(60, 2, 1, B = 199)
  x <- matrix(rnorm(120), ncol = 2)
  y <- rnorm(60)
  res <- hrank_dcov_test(x, y, null = null)
  expect_equal(res$p.value,
               (1 + sum(null$statistics >= res$statistic - 1e-10)) / 200)
  expect_identical(hrank_dcov_test(x, y, null = res$null)$p.value, res$p.value)
  expect_error(hrank_dcov_test(x, cbind(y, y), null = null), "was computed for")
  expect_error(hrank_dcov_test(x, y, null = list()), "hrank_dcov_null")
})

test_that("random grids give valid tests", {
  set.seed(4)
  res <- hrank_dcov_test(matrix(rnorm(150), ncol = 3), rnorm(50), type = "random", B = 99)
  expect_true(res$p.value > 0 && res$p.value <= 1)
})

test_that("the permutation null distribution is centered", {
  set.seed(5)
  null <- hrank_dcov_null(80, 2, 2, B = 2000)
  expect_lt(abs(mean(null$statistics)), 4 * sd(null$statistics) / sqrt(2000))
})

test_that("permutation and asymptotic null distributions agree", {
  set.seed(6)
  perm <- hrank_dcov_null(300, 2, 3, B = 3000)
  asym <- hrank_dcov_null(300, 2, 3, method = "asymptotic")
  qp <- quantile(perm, c(0.5, 0.9, 0.95))
  qa <- quantile(asym, c(0.5, 0.9, 0.95))
  expect_equal(unname(qa), unname(qp), tolerance = 0.1)
})

test_that("Imhof's method matches simulation", {
  set.seed(7)
  w <- list(w = c(1, 0.5, 0.2, 0.1), df = c(1, 1, 1, 3))
  sims <- colSums(w$w * matrix(rchisq(4e5, df = rep(w$df, 1e5)), nrow = 4)) - sum(w$w * w$df)
  for (t in c(-0.5, 0, 1, 3)) {
    expect_equal(wchisq_upper(t, w), mean(sims > t), tolerance = 0.01)
  }
})

test_that("asymptotic quantiles invert the tail probabilities", {
  asym <- hrank_dcov_null(100, 2, 2, method = "asymptotic")
  q <- quantile(asym, c(0.9, 0.99))
  expect_equal(wchisq_upper(q[[1]], asym$weights), 0.1, tolerance = 1e-6)
  expect_equal(wchisq_upper(q[[2]], asym$weights), 0.01, tolerance = 1e-6)
  expect_output(print(asym), "asymptotic")
})
