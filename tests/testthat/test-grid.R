grid_is_symmetric <- function(g) {
  key <- function(m) apply(round(m, 10), 1, paste, collapse = ",")
  setequal(key(g), key(-g)) && all(table(key(g)) == table(key(-g))[names(table(key(g)))])
}

test_that("grid factorization satisfies n = n_R n_S + n_0 with n_S even and n_0 < 2 n_R", {
  for (n in c(4, 17, 100, 257)) {
    for (d in 1:4) {
      g <- co_grid(n, d)
      expect_equal(dim(g), c(n, d))
      n_R <- attr(g, "n_R")
      n_S <- attr(g, "n_S")
      n_0 <- attr(g, "n_0")
      expect_equal(n_R * n_S + n_0, n)
      expect_equal(n_S %% 2, 0)
      expect_lt(n_0, if (d == 1) 2 else 2 * n_R)
      expect_true(all(rowSums(g^2) < 1))
      expect_true(grid_is_symmetric(g))
    }
  }
})

test_that("univariate grid is symmetric and evenly spaced", {
  g <- co_grid(7, 1)
  expect_equal(g[, 1], c(-3:3) / 4)
  expect_equal(co_grid(6, 1)[, 1], c(-3:-1, 1:3) / 4)
})

test_that("radii and directions are as specified", {
  g <- co_grid(100, 2, n_R = 4)
  expect_equal(attr(g, "n_S"), 24)
  expect_equal(attr(g, "n_0"), 4)
  expect_equal(sort(unique(round(sqrt(rowSums(g[1:96, ]^2)), 12))), round((1:4) / 5, 12))
  angles <- sort(atan2(g[1:24, 2], g[1:24, 1]) %% (2 * pi))
  expect_equal(diff(angles), rep(2 * pi / 24, 23))
})

test_that("extra points come in antipodal pairs plus the origin when n_0 is odd", {
  g4 <- co_grid(100, 3, n_R = 4)   # n_S = 24, n_0 = 4: two pairs
  r <- sqrt(rowSums(g4^2))
  expect_equal(sum(abs(r - 1 / 10) < 1e-12), 4)
  expect_equal(sum(r == 0), 0)
  expect_true(grid_is_symmetric(g4))

  g5 <- co_grid(101, 3, n_R = 4)   # n_0 = 5: two pairs and the origin
  r <- sqrt(rowSums(g5^2))
  expect_equal(sum(abs(r - 1 / 10) < 1e-12), 4)
  expect_equal(sum(r == 0), 1)
  expect_true(grid_is_symmetric(g5))

  g1 <- co_grid(97, 3, n_R = 4)    # n_0 = 1: the origin only
  expect_equal(sum(rowSums(g1^2) == 0), 1)
  expect_equal(anyDuplicated(round(g5, 12)), 0)
})

test_that("regular grids are deterministic and random grids are not", {
  expect_identical(co_grid(50, 3), co_grid(50, 3))
  set.seed(1)
  g1 <- co_grid(50, 3, type = "random")
  g2 <- co_grid(50, 3, type = "random")
  expect_false(isTRUE(all.equal(g1, g2)))
  expect_equal(attr(g1, "n_R"), 2)
  expect_equal(attr(g1, "n_S"), 24)
  expect_equal(sqrt(rowSums(g1^2))[1:24], rep(1 / 3, 24))
  expect_true(grid_is_symmetric(g1))
})

test_that("regular directions are roughly uniform on the sphere", {
  s <- 2 * co_grid(4000, 3, n_R = 1)
  expect_equal(rowSums(s^2), rep(1, 4000))
  expect_equal(colMeans(s), rep(0, 3))
  expect_lt(max(abs(crossprod(s) / 4000 - diag(3) / 3)), 0.01)
})

test_that("invalid arguments are rejected", {
  expect_error(co_grid(1, 2), "`n`")
  expect_error(co_grid(10, 0), "`d`")
  expect_error(co_grid(10, 2, n_R = 4), "too large")
  expect_error(co_grid(10, 2, n_R = 1.5), "`n_R`")
})
