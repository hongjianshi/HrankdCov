#include <R.h>
#include <Rinternals.h>
#include <float.h>

/*
 * Linear sum assignment: find the permutation sigma minimizing
 * sum_i cost[i, sigma(i)] for a dense n x n cost matrix.
 *
 * Shortest augmenting path version of the Hungarian algorithm with dual
 * potentials (Kuhn-Munkres with Dijkstra-type updates), O(n^3) in the worst
 * case.  Rows are added one at a time; arrays are 1-based with index 0 used
 * as a sentinel column.
 *
 * Returns an integer vector `sigma` (1-based) with row i matched to column
 * sigma[i].
 */
SEXP C_lsap(SEXP cost_) {
  if (!isReal(cost_) || !isMatrix(cost_) || nrows(cost_) != ncols(cost_))
    error("'cost' must be a square double matrix");
  const int n = nrows(cost_);
  const double *c0 = REAL(cost_);

  /* Row-major copy so that the inner loop reads contiguous memory. */
  double *c = (double *) R_alloc((size_t) n * n, sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      c[(size_t) i * n + j] = c0[(size_t) j * n + i];

  double *u = (double *) R_alloc(n + 1, sizeof(double));
  double *v = (double *) R_alloc(n + 1, sizeof(double));
  double *minv = (double *) R_alloc(n + 1, sizeof(double));
  int *match = (int *) R_alloc(n + 1, sizeof(int));  /* column -> row */
  int *way = (int *) R_alloc(n + 1, sizeof(int));
  char *used = (char *) R_alloc(n + 1, sizeof(char));

  for (int j = 0; j <= n; j++) {
    u[j] = 0.0;
    v[j] = 0.0;
    match[j] = 0;
    way[j] = 0;
  }

  for (int i = 1; i <= n; i++) {
    if ((i & 127) == 0) R_CheckUserInterrupt();
    match[0] = i;
    int j0 = 0;
    for (int j = 0; j <= n; j++) {
      minv[j] = DBL_MAX;
      used[j] = 0;
    }
    do {
      used[j0] = 1;
      const int i0 = match[j0];
      const double *row = c + (size_t) (i0 - 1) * n - 1;
      double delta = DBL_MAX;
      int j1 = 0;
      for (int j = 1; j <= n; j++) {
        if (used[j]) continue;
        const double cur = row[j] - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; j++) {
        if (used[j]) {
          u[match[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (match[j0] != 0);
    do {
      const int j1 = way[j0];
      match[j0] = match[j1];
      j0 = j1;
    } while (j0 != 0);
  }

  SEXP ans = PROTECT(allocVector(INTSXP, n));
  int *sigma = INTEGER(ans);
  for (int j = 1; j <= n; j++) sigma[match[j] - 1] = j;
  UNPROTECT(1);
  return ans;
}
