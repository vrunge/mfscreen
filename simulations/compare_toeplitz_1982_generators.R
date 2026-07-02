# Compare two generators of an n x 1982 Gaussian Toeplitz block
#
# Old implementation:
#   dense Toeplitz covariance + mvtnorm::rmvnorm()
#
# New implementation:
#   stationary AR(1) recursion
#
# Both target:
#   Cov(X_j, X_k) = sigma^2 * rho^|j-k|.

# -------------------------------------------------------------------------
# 1. Parameters
# -------------------------------------------------------------------------
n <- 500L
d <- 1982L
rho <- 0.70
sigma <- 1
B <- 3L  # Number of timing repetitions. Increase for a more stable benchmark.

# Required only for the old dense implementation:
# install.packages("mvtnorm")
stopifnot(requireNamespace("mvtnorm", quietly = TRUE))

# -------------------------------------------------------------------------
# 2. Old dense Toeplitz implementation, extracted from the original Rmd
# -------------------------------------------------------------------------
simulate_toeplitz_dense_old <- function(
    n,
    d,
    mu = rep(0, d),
    rho = 0.70,
    sigma = 1
) {
  create_cor <- function(i) rho^(i - 1L) * sigma^2

  covariance_vector <- sapply(seq_len(d), create_cor)
  Sigma <- toeplitz(covariance_vector)

  # The original code used:
  # Sigma <- lqmm::make.positive.definite(toeplitz(covariance_vector))
  #
  # For abs(rho) < 1, this AR(1) Toeplitz covariance is already positive
  # definite. Therefore, no correction is required here.

  mvtnorm::rmvnorm(
    n = n,
    mean = mu,
    sigma = Sigma
  )
}

# -------------------------------------------------------------------------
# 3. New exact stationary-AR(1) implementation
# -------------------------------------------------------------------------
simulate_toeplitz_ar1_fast <- function(n, d, rho, sigma = 1) {
  stopifnot(d >= 1L, abs(rho) < 1, sigma > 0)

  X <- matrix(0, nrow = n, ncol = d)

  # Stationary initial distribution: Var(X_1) = sigma^2.
  X[, 1L] <- rnorm(n, mean = 0, sd = sigma)

  if (d > 1L) {
    innovation_sd <- sigma * sqrt(1 - rho^2)

    innovations <- matrix(
      rnorm(n * (d - 1L), mean = 0, sd = innovation_sd),
      nrow = n,
      ncol = d - 1L
    )

    for (j in 2:d) {
      X[, j] <- rho * X[, j - 1L] + innovations[, j - 1L]
    }
  }

  X
}

# -------------------------------------------------------------------------
# 4. Check that both generators target the same covariance structure
# -------------------------------------------------------------------------
set.seed(1)
X_old <- simulate_toeplitz_dense_old(n, d, rho = rho, sigma = sigma)

set.seed(1)
X_fast <- simulate_toeplitz_ar1_fast(n, d, rho = rho, sigma = sigma)

lags <- 0:5

covariance_check <- data.frame(
  lag = lags,
  theoretical_covariance = sigma^2 * rho^lags,
  old_dense_empirical_covariance = vapply(
    lags,
    function(h) {
      if (h == 0L) {
        var(X_old[, 1L])
      } else {
        cov(X_old[, 1L], X_old[, 1L + h])
      }
    },
    numeric(1)
  ),
  ar1_empirical_covariance = vapply(
    lags,
    function(h) {
      if (h == 0L) {
        var(X_fast[, 1L])
      } else {
        cov(X_fast[, 1L], X_fast[, 1L + h])
      }
    },
    numeric(1)
  )
)

print(covariance_check)

# -------------------------------------------------------------------------
# 5. Timing benchmark
# -------------------------------------------------------------------------
# Garbage collection before each timing is useful because these are large
# temporary objects.

benchmark_one <- function(generator, B) {
  elapsed <- numeric(B)

  for (b in seq_len(B)) {
    gc()
    set.seed(1000L + b)

    elapsed[b] <- system.time({
      X <- generator(n = n, d = d, rho = rho, sigma = sigma)
    })[["elapsed"]]

    rm(X)
  }

  elapsed
}

old_elapsed <- benchmark_one(simulate_toeplitz_dense_old, B = B)
fast_elapsed <- benchmark_one(simulate_toeplitz_ar1_fast, B = B)

benchmark_summary <- data.frame(
  method = c("Dense Toeplitz + rmvnorm", "Stationary AR(1) recursion"),
  median_elapsed_seconds = c(median(old_elapsed), median(fast_elapsed)),
  mean_elapsed_seconds = c(mean(old_elapsed), mean(fast_elapsed)),
  min_elapsed_seconds = c(min(old_elapsed), min(fast_elapsed)),
  max_elapsed_seconds = c(max(old_elapsed), max(fast_elapsed))
)

benchmark_summary$speedup_vs_dense <-
  benchmark_summary$median_elapsed_seconds[1L] /
  benchmark_summary$median_elapsed_seconds

print(benchmark_summary, row.names = FALSE)

# -------------------------------------------------------------------------
# 6. Memory estimates
# -------------------------------------------------------------------------
bytes_per_double <- 8
dense_covariance_mb <- d^2 * bytes_per_double / 1024^2
output_matrix_mb <- n * d * bytes_per_double / 1024^2
innovation_matrix_mb <- n * (d - 1L) * bytes_per_double / 1024^2

memory_estimate <- data.frame(
  object = c(
    "Dense 1982 x 1982 covariance matrix",
    "Generated n x 1982 output matrix",
    "AR(1) innovation matrix"
  ),
  approximate_megabytes = c(
    dense_covariance_mb,
    output_matrix_mb,
    innovation_matrix_mb
  )
)

print(memory_estimate, row.names = FALSE)

cat(
  "\nComplexity comparison:\n",
  "- Dense method: O(d^3) for a Cholesky factorization and O(n d^2) for sampling.\n",
  "- AR(1) method: O(n d), with no d x d covariance matrix or factorization.\n",
  sep = ""
)

