# Eigenvalues of the limiting null distribution.
#
# Under independence, n dCov^2_n of the scored center-outward ranks converges
# to sum_{k, l} lambda_k mu_l (xi_kl^2 - 1), where lambda_k (resp. mu_l) are
# the eigenvalues of the integral operator with the doubly centered kernel
#   k(u, u') = || J(r) s - J(r') s' ||,   u = r s,  u' = r' s',
# on L^2 of the spherical uniform distribution U_d (R uniform on [0, 1],
# S uniform on the unit sphere, independent).
#
# The kernel is rotation invariant, so by the Funk-Hecke formula it acts on
# functions of the form f(r) Y_l(s), Y_l a spherical harmonic of degree l, as
# the one-dimensional operator with kernel
#   kappa_l(r, r') = c_d int_0^pi k(r, r', cos theta) P_l(cos theta) sin^{d-2}(theta) dtheta,
# P_l the Legendre polynomial in dimension d.  Centering only affects l = 0.
# The eigenvalues of kappa_l (each with multiplicity N(d, l), the dimension of
# the degree-l harmonics) are computed by Nystrom's method with
# Gauss-Legendre quadrature in r and theta.  The sum and the sum of squares of
# all eigenvalues are available in closed form from the centered kernel and
# are used for the tail of the eigenvalue sequence.

gauss_legendre <- function(n, a = 0, b = 1) {
  i <- seq_len(n - 1)
  beta <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- beta
  J[cbind(i + 1, i)] <- beta
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(nodes = (b - a) / 2 * e$values[o] + (a + b) / 2,
       weights = (b - a) * e$vectors[1, o]^2)
}

# Legendre polynomials in dimension d >= 2 (normalized Gegenbauer polynomials
# with parameter (d - 2) / 2), P_l(t) for l = 0, ..., L, as columns.
dim_legendre <- function(t, L, d) {
  P <- matrix(0, length(t), L + 1)
  P[, 1] <- 1
  if (L == 0) return(P)
  if (d == 2) {
    theta <- acos(pmin(1, pmax(-1, t)))
    for (l in seq_len(L)) P[, l + 1] <- cos(l * theta)
    return(P)
  }
  a <- (d - 2) / 2
  P[, 2] <- t
  if (L >= 2) {
    for (l in 2:L) {
      P[, l + 1] <- (2 * (l + a - 1) * t * P[, l] - (l - 1) * P[, l - 1]) / (2 * a + l - 1)
    }
  }
  P
}

# Dimension of the space of spherical harmonics of degree l in R^d.
harmonic_mult <- function(d, l) {
  if (d == 1) return(rep(1, length(l)))
  if (d == 2) return(ifelse(l == 0, 1, 2))
  (2 * l + d - 2) / (d - 2) * choose(l + d - 3, l)
}

# Returns a list with the largest `keep` eigenvalues (positive convention,
# decreasing, with multiplicities), and the exact sum and sum of squares of
# all eigenvalues.
limit_eigen <- function(d, J, n_r = 200, n_theta = 400, L = 100, keep = 1000,
                        tol = 1e-8) {
  gr <- gauss_legendre(n_r, 0, 1)
  r <- gr$nodes
  w <- gr$weights
  Jr <- J(r)
  if (length(Jr) != n_r || anyNA(Jr) || any(Jr < 0) || any(!is.finite(Jr))) {
    stop("The score function must return finite nonnegative values on (0, 1).", call. = FALSE)
  }
  kernel <- function(t) sqrt(pmax(outer(Jr^2, Jr^2, "+") - 2 * t * outer(Jr, Jr), 0))

  if (d == 1) {
    # Directions +1 and -1; "harmonics" are the even (l = 0) and odd (l = 1) parts.
    L <- 1
    kp <- kernel(1)
    km <- kernel(-1)
    kappa <- array(c((kp + km) / 2, (kp - km) / 2), c(n_r, n_r, 2))
  } else {
    gt <- gauss_legendre(n_theta, 0, pi)
    theta <- gt$nodes
    c_d <- gamma(d / 2) / (sqrt(pi) * gamma((d - 1) / 2))
    C <- dim_legendre(cos(theta), L, d) * (c_d * gt$weights * sin(theta)^(d - 2))
    kappa <- array(0, c(n_r, n_r, L + 1))
    chunk <- 50L
    for (start in seq(1L, n_theta, by = chunk)) {
      idx <- start:min(start + chunk - 1L, n_theta)
      K <- vapply(idx, function(i) as.vector(kernel(cos(theta[i]))), numeric(n_r^2))
      kappa <- kappa + array(K %*% C[idx, , drop = FALSE], c(n_r, n_r, L + 1))
    }
  }

  # Center the degree-zero block: k - m(r) - m(r') + c.
  m <- as.vector(kappa[, , 1] %*% w)
  cc <- sum(w * m)
  kappa[, , 1] <- kappa[, , 1] - outer(m, m, "+") + cc
  sum_w <- cc
  sum_w2 <- 2 * sum(w * Jr^2) - 2 * sum(w * m^2) + cc^2

  sw <- sqrt(w)
  values <- numeric(0)
  mult <- numeric(0)
  for (l in 0:L) {
    A <- sw * kappa[, , l + 1] * rep(sw, each = n_r)
    ev <- -eigen(A, symmetric = TRUE, only.values = TRUE)$values
    ev <- ev[ev > tol * sum_w]
    if (length(ev) > 0) {
      values <- c(values, ev)
      mult <- c(mult, rep(harmonic_mult(d, l), length(ev)))
    }
    # The degree-zero block can be empty (sign scores); higher degrees decay.
    if (l > 0 && (length(ev) == 0 || max(ev) < 1e-6 * sum_w)) break
  }
  o <- order(values, decreasing = TRUE)
  values <- rep(values[o], mult[o])
  list(values = values[seq_len(min(keep, length(values)))], sum_w = sum_w, sum_w2 = sum_w2)
}

limit_cache <- new.env(parent = emptyenv())

# Eigenvalues for a built-in or user-supplied score in dimension d, using the
# precomputed values in R/sysdata.rda when available.
limit_eigenvalues_for <- function(score, d) {
  key <- if (is.function(score)) NULL else paste(score, d)
  if (!is.null(key)) {
    stored <- limit_eigenvalues[[key]]
    if (!is.null(stored)) return(stored)
    cached <- limit_cache[[key]]
    if (!is.null(cached)) return(cached)
  }
  e <- limit_eigen(d, score_function(score, d))
  if (!is.null(key)) assign(key, e, envir = limit_cache)
  e
}
