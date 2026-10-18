test_that("the statistic is n times the U-statistic distance covariance of the ranks", {
  skip_if_not_installed("energy")
  set.seed(1)
  x <- matrix(rnorm(240), ncol = 3)
  y <- cbind(x[, 1]^2, rnorm(80))
  for (score in c("rank", "normal", "sign")) {
    sx <- co_scores(co_ranks(x), score)
    sy <- co_scores(co_ranks(y), score)
    expect_equal(hrank_dcov(x, y, score = score), 80 * unname(energy::dcovU(sx, sy)))
  }
})

test_that("the statistic is distribution-free under monotone transformations in d = 1", {
  set.seed(2)
  x <- rnorm(60)
  y <- x + rnorm(60)
  expect_equal(hrank_dcov(x, y), hrank_dcov(exp(x), y^3))
})

test_that("n_R can be set separately for x and y", {
  set.seed(3)
  x <- matrix(rnorm(200), ncol = 2)
  y <- matrix(rnorm(200), ncol = 2)
  expect_equal(hrank_dcov(x, y, n_R = c(2, 2)), hrank_dcov(x, y, n_R = 2))
  expect_false(isTRUE(all.equal(hrank_dcov(x, y, n_R = c(2, 5)), hrank_dcov(x, y, n_R = 2))))
  expect_error(hrank_dcov(x, y, n_R = 1:3), "length 1 or 2")
})

test_that("sample sizes must agree", {
  expect_error(hrank_dcov(rnorm(10), rnorm(11)), "same number")
  expect_error(hrank_dcov(rnorm(3), rnorm(3)), "At least 4")
})
