total_cost <- function(x, g, index) sum((x - g[index, , drop = FALSE])^2)

test_that("ranks are a permutation of the grid", {
  set.seed(1)
  x <- matrix(rnorm(300), ncol = 3)
  r <- co_ranks(x)
  expect_s3_class(r, "co_ranks")
  expect_equal(sort(r$index), 1:100)
  expect_equal(r$F, r$grid[r$index, ])
  expect_equal(sort(r$rank), sort(round((attr(r$grid, "n_R") + 1) * sqrt(rowSums(r$grid^2)), 10)))
  norms <- sqrt(rowSums(r$sign^2))
  expect_true(all(abs(norms - 1) < 1e-12 | norms == 0))
})

test_that("the assignment is optimal", {
  skip_if_not_installed("clue")
  set.seed(2)
  for (d in 2:4) {
    x <- matrix(rt(120 * d, df = 2), ncol = d)
    r <- co_ranks(x)
    D <- outer(rowSums(x^2), rowSums(r$grid^2), "+") - 2 * tcrossprod(x, r$grid)
    best <- as.vector(clue::solve_LSAP(D - min(D)))
    expect_equal(total_cost(x, r$grid, r$index), total_cost(x, r$grid, best))
  }
})

test_that("univariate center-outward ranks preserve the order of the data", {
  set.seed(3)
  z <- rnorm(31)
  r <- co_ranks(z)
  expect_equal(order(r$F[, 1]), order(z))
  expect_equal(r$rank, abs(rank(z) - 16))
})

test_that("ranks are invariant to translation and positive scaling", {
  set.seed(4)
  x <- matrix(rexp(200), ncol = 2)
  expect_equal(co_ranks(x)$index, co_ranks(3 * x + 7)$index)
})

test_that("ranks are equivariant under permutations of the observations", {
  set.seed(5)
  x <- matrix(rnorm(150), ncol = 3)
  p <- sample(50)
  expect_equal(co_ranks(x[p, ])$F, co_ranks(x)$F[p, ])
})

test_that("a user-supplied grid is used", {
  set.seed(6)
  x <- matrix(rnorm(80), ncol = 2)
  g <- co_grid(40, 2, type = "random")
  r <- co_ranks(x, grid = g)
  expect_equal(r$grid, g, ignore_attr = TRUE)
  expect_error(co_ranks(x, grid = g[-1, ]), "40 x 2")
  expect_error(co_ranks(x, grid = 2 * g), "unit ball")
})

test_that("data frames are accepted and bad input is rejected", {
  set.seed(7)
  df <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_equal(co_ranks(df)$index, co_ranks(as.matrix(df))$index)
  expect_error(co_ranks(data.frame(a = letters)), "numeric")
  expect_error(co_ranks(c(1, NA, 3)), "missing")
})

test_that("print and plot methods work", {
  set.seed(8)
  r <- co_ranks(matrix(rnorm(100), ncol = 2))
  expect_output(print(r), "Center-outward ranks of 50 observations")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(r))
  expect_invisible(plot(co_ranks(rnorm(10))))
  expect_error(plot(co_ranks(matrix(rnorm(30), ncol = 3))), "dimensions 1 and 2")
})

test_that("scores have the right norms", {
  set.seed(9)
  r <- co_ranks(matrix(rnorm(303), ncol = 3))
  u <- sqrt(rowSums(r$F^2))
  expect_equal(co_scores(r, "rank"), r$F)
  expect_equal(sqrt(rowSums(co_scores(r, "normal")^2)), sqrt(qchisq(u, 3)))
  expect_equal(sqrt(rowSums(co_scores(r, "sign")^2)), as.numeric(u > 0))
  expect_equal(co_scores(r, function(u) u^2), r$F * u)
  expect_error(co_scores(r, function(u) -u), "nonnegative")
  expect_error(co_scores(r, "spearman"))
})
