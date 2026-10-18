#include <R.h>
#include <Rinternals.h>

/*
 * Permutation distribution of the (unnormalized) U-statistic
 *   sum_{i != j} A[i, j] * B[pi(i), pi(j)]
 * for symmetric n x n matrices A and B with zero diagonals (U-centered
 * distance matrices) and the permutations stored in the columns of the
 * integer matrix `perms` (1-based).
 */
SEXP C_perm_dcov(SEXP A_, SEXP B_, SEXP perms_) {
  const int n = nrows(A_);
  if (nrows(B_) != n || ncols(A_) != n || ncols(B_) != n || nrows(perms_) != n)
    error("dimension mismatch");
  const int nperm = ncols(perms_);
  const double *A = REAL(A_);
  const double *B = REAL(B_);
  const int *P = INTEGER(perms_);

  SEXP ans = PROTECT(allocVector(REALSXP, nperm));
  double *out = REAL(ans);
  int *pi = (int *) R_alloc(n, sizeof(int));

  for (int b = 0; b < nperm; b++) {
    if ((b & 15) == 0) R_CheckUserInterrupt();
    for (int i = 0; i < n; i++) pi[i] = P[(size_t) b * n + i] - 1;
    double s = 0.0;
    for (int j = 0; j < n; j++) {
      const double *Aj = A + (size_t) j * n;
      const double *Bj = B + (size_t) pi[j] * n;
      for (int i = j + 1; i < n; i++) s += Aj[i] * Bj[pi[i]];
    }
    out[b] = 2.0 * s;
  }

  UNPROTECT(1);
  return ans;
}
