# Checks of the limiting null distribution against Shi, Drton and Han (2022).

paper_critical_values <- list(
  # (p, q) = c(alpha = 0.1, 0.05, 0.01), Table C.1 of the supplement,
  # tabulated with accuracy 5e-3.
  "1,1" = c(0.306, 0.490, 0.945),
  "1,2" = c(0.215, 0.320, 0.563),
  "2,2" = c(0.145, 0.205, 0.338),
  "2,3" = c(0.114, 0.159, 0.255),
  "3,3" = c(0.090, 0.124, 0.194),
  "5,5" = c(0.052, 0.070, 0.105),
  "7,7" = c(0.037, 0.049, 0.073),
  "10,10" = c(0.027, 0.035, 0.052)
)

test_that("asymptotic critical values agree with Table C.1 of Shi, Drton and Han (2022)", {
  for (key in names(paper_critical_values)) {
    pq <- as.integer(strsplit(key, ",")[[1]])
    null <- hrank_dcov_null(100, pq[1], pq[2], method = "asymptotic")
    q <- unname(quantile(null, c(0.9, 0.95, 0.99)))
    expect_lt(max(abs(q - paper_critical_values[[key]])), 0.003,
              label = sprintf("max |difference| for (p, q) = (%d, %d)", pq[1], pq[2]))
  }
})

test_that("the p = 1 eigenvalues match the closed form 4 / (pi^2 j^2)", {
  e <- limit_eigen(1, function(u) u)
  expect_equal(e$values[1:10], 4 / (pi^2 * (1:10)^2), tolerance = 1e-3)
  expect_equal(e$sum_w, 2 / 3, tolerance = 1e-5)      # E|U - U'| for U uniform on [-1, 1]
  expect_equal(e$sum_w2, 16 / 90, tolerance = 1e-4)
  expect_equal(limit_eigenvalues[["rank 1"]]$values[1:10], e$values[1:10])
})

test_that("limiting eigenvalues have the multiplicities of spherical harmonics", {
  e3 <- limit_eigenvalues[["rank 3"]]
  expect_equal(e3$values[1:3], rep(e3$values[1], 3))    # degree 1 in R^3: multiplicity 3
  expect_lt(e3$values[4], e3$values[1])
  e5 <- limit_eigenvalues[["normal 5"]]
  expect_equal(e5$values[1:5], rep(e5$values[1], 5))
})

test_that("the eigenvalue sums are the moments of the limiting distribution", {
  # sum of eigenvalues = E||U - U'|| and sum of squares from the kernel
  set.seed(1)
  u <- 2 * co_grid(20000, 3, n_R = 1, type = "random") * runif(20000)   # spherical uniform
  u <- u[sample(20000), ]
  e3 <- limit_eigenvalues[["rank 3"]]
  mc <- mean(sqrt(rowSums((u[1:10000, ] - u[10001:20000, ])^2)))
  expect_equal(e3$sum_w, mc, tolerance = 0.02)
  expect_lt(sum(e3$values), e3$sum_w)
  expect_gt(sum(e3$values) / e3$sum_w, 0.9)
  expect_equal(sum(e3$values^2), e3$sum_w2, tolerance = 1e-4)
})

test_that("limiting eigenvalues are stable under refinement of the quadrature", {
  for (J in list(function(u) u, function(u) sqrt(qchisq(u, 2)))) {
    a <- limit_eigen(2, J)
    b <- limit_eigen(2, J, n_r = 400, n_theta = 800)
    expect_equal(a$values[1:20], b$values[1:20], tolerance = 1e-4)
    expect_equal(a$sum_w, b$sum_w, tolerance = 1e-5)
  }
})

test_that("user-supplied score functions and large dimensions are computed on the fly", {
  e <- limit_eigenvalues_for(function(u) u^2, 2)
  expect_true(all(diff(e$values) <= 0))
  e12 <- limit_eigenvalues_for("rank", 12)
  expect_identical(limit_eigenvalues_for("rank", 12), e12)   # cached
  expect_error(limit_eigen(2, function(u) -u), "nonnegative")
})
