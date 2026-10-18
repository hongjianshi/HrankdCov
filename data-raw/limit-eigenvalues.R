# Precomputes the eigenvalues of the limiting null distribution for the
# built-in score functions and dimensions 1 to 10 (see R/limit.R), and stores
# them in R/sysdata.rda.  Rerun after any change to limit_eigen().
#
#   Rscript data-raw/limit-eigenvalues.R

pkgload::load_all(".", quiet = TRUE)

limit_eigenvalues <- list()
for (score in c("rank", "normal", "sign")) {
  for (d in 1:10) {
    key <- paste(score, d)
    t0 <- proc.time()[[3]]
    limit_eigenvalues[[key]] <- limit_eigen(d, score_function(score, d))
    cat(sprintf("%-10s largest %.5f  sum %.5f  (%.1fs)\n", key,
                limit_eigenvalues[[key]]$values[1], limit_eigenvalues[[key]]$sum_w,
                proc.time()[[3]] - t0))
  }
}

usethis::use_data(limit_eigenvalues, internal = TRUE, overwrite = TRUE, compress = "xz")
