# HrankdCov

<!-- badges: start -->
[![R-CMD-check](https://github.com/hongjianshi/HrankdCov/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hongjianshi/HrankdCov/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**HrankdCov** provides distribution-free and consistent tests of independence
between two random vectors of arbitrary dimensions.  It combines the
distance covariance with the *center-outward ranks and signs* of Hallin et
al. (2021), which extend univariate ranks to random vectors via optimal
transport.  The methodology is developed in

* Shi, H., Drton, M. and Han, F. (2022). Distribution-free consistent
  independence tests via center-outward ranks and signs. *Journal of the
  American Statistical Association*, 117(537), 395–410.
  [doi:10.1080/01621459.2020.1782223](https://doi.org/10.1080/01621459.2020.1782223)
* Shi, H., Hallin, M., Drton, M. and Han, F. (2022). On universally
  consistent and fully distribution-free rank tests of vector independence.
  *The Annals of Statistics*, 50(4), 1933–1959.
  [doi:10.1214/21-AOS2151](https://doi.org/10.1214/21-AOS2151)

Under independence, the null distribution of the test statistic is the same
for all absolutely continuous distributions.  The test therefore needs no
moment assumptions, is robust to heavy tails and outliers, and its critical
values can be computed once and reused.

## Installation

```r
# install.packages("remotes")
remotes::install_github("hongjianshi/HrankdCov")
```

The package contains C code, so installing from source requires a compiler
(Rtools on Windows, Xcode command line tools on macOS).

## Example

```r
library(HrankdCov)

set.seed(1)
n <- 300
x <- matrix(rnorm(2 * n), ncol = 2)
y <- cbind(cos(2 * x[, 1]), x[, 2]^2) + matrix(rnorm(2 * n, sd = 0.5), ncol = 2)

hrank_dcov_test(x, y)
#>
#> 	Center-outward rank distance covariance test of independence (rank
#> 	(Wilcoxon) scores; permutation (Monte Carlo) null distribution with
#> 	999 replicates)
#>
#> data:  x and y
#> n dCov^2 = 2.4356, p-value = 0.001
#> sample estimates:
#>     dCov^2
#> 0.008118526
```

Normal (van der Waerden) or sign scores, and p-values from the limiting
distribution instead of permutations:

```r
hrank_dcov_test(x, y, score = "normal", method = "asymptotic")
```

Because the null distribution is distribution-free, it can be computed once
and reused, for example in a simulation study:

```r
null <- hrank_dcov_null(n, p = 2, q = 2, method = "asymptotic")
x0 <- matrix(rt(2 * n, df = 1), ncol = 2)   # Cauchy data
y0 <- matrix(rt(2 * n, df = 1), ncol = 2)
hrank_dcov_test(x0, y0, null = null)$p.value
```

## Functions

| Function | Purpose |
|---|---|
| `hrank_dcov_test()` | Test of independence (returns an `htest` object) |
| `hrank_dcov()` | Test statistic: *n* times the distance covariance of scored center-outward ranks |
| `hrank_dcov_null()` | Distribution-free null distribution (permutation or asymptotic), with a `quantile()` method for critical values |
| `co_ranks()` | Center-outward ranks and signs, with `print()` and `plot()` methods |
| `co_grid()` | The reference grid in the unit ball |
| `co_scores()` | Rank, normal (van der Waerden), sign or user-supplied scores |

See `vignette("HrankdCov")` for an introduction.  The simulation study of
Shi, Drton and Han (2022) can be reproduced with
[`inst/scripts/jasa-simulation.R`](inst/scripts/jasa-simulation.R).

## Citation

Run `citation("HrankdCov")` for BibTeX entries.

## License

GPL (>= 3)
