# HrankdCov 0.1.0

* First release as an R package, replacing the original scripts
  `discov-Hallin-stat.R` and `simulation.R`.  `discov.Hallin.stat(X, Y)`
  corresponds to `hrank_dcov(X, Y, type = "random")`, up to the randomly
  drawn directions and the placement of the extra grid points near the
  origin.
* New functions:
  * `co_grid()`, `co_ranks()` and `co_scores()` for center-outward grids,
    ranks, signs and scores, with `print()` and `plot()` methods;
  * `hrank_dcov()` for the test statistic, now with rank, normal (van der
    Waerden), sign and user-supplied score functions;
  * `hrank_dcov_test()`, returning an `"htest"` object with p-values from the
    exact permutation null distribution or from the limiting weighted
    chi-square distribution;
  * `hrank_dcov_null()` for computing the distribution-free null
    distribution once and reusing it across tests.
* The grid is symmetric about the origin: the directions come in antipodal
  pairs (so `n_S` is even), and the `n_0` extra points are antipodal pairs
  at radius `1 / (2 (n_R + 1))` plus the origin when `n_0` is odd.  The
  grid is deterministic by default (`type = "regular"`).
* The eigenvalues of the limiting null distribution are computed by
  reducing the rotation-invariant kernel to one-dimensional radial integral
  equations via spherical harmonics, instead of from a finite grid.  They
  reproduce the critical values tabulated in Shi, Drton and Han (2022) to
  within about 1e-3, and are precomputed for the built-in scores and
  dimensions up to 10.
* The optimal assignment is now computed by a C implementation of the
  Hungarian algorithm, replacing `clue` and `adagio`.
* The simulation script of Shi, Drton and Han (2022) is in
  `inst/scripts/jasa-simulation.R`.
