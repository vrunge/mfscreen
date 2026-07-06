# Compare screening_test_matrix_R() and screening_test_matrix()
# for several n / p ratios, with repeated timings and a plot.

set.seed(123)

q <- 0.15
n_rep <- 10

# Keep p fixed and vary n so that n / p changes.
p <- 10000
n_values <- 10*c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000)

results <- data.frame(
  n = integer(),
  p = integer(),
  ratio_np = numeric(),
  replication = integer(),
  method = character(),
  time_seconds = numeric()
)

for (n in n_values) {

  cat("Running n =", n, ", p =", p, "\n")

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- 3 * X[, 1] - 2 * X[, 2] + rnorm(n)

  # Warm-up calls: exclude one-time overhead from timings.
  #invisible(screening_test_matrix_R(X[, 1:10], y, q = q))
  #invisible(screening_test_matrix(X[, 1:10], y, q = q))

  for (r in seq_len(n_rep)) {

    time_R <- system.time({
      result_R <- screening_test_matrix_R(X, y, q = q)
    })["elapsed"]

    time_cpp <- system.time({
      result_cpp <- screening_test_matrix(X, y, q = q)
    })["elapsed"]

    # Optional consistency check on the first replication only.
    if (r == 1L) {
      same_selection <- identical(
        unname(which(result_R$selected)),
        result_cpp$selected_indices
      )

      cat("  Same selected variables:", same_selection, "\n")
    }

    results <- rbind(
      results,
      data.frame(
        n = n,
        p = p,
        ratio_np = n / p,
        replication = r,
        method = "R",
        time_seconds = unname(time_R)
      ),
      data.frame(
        n = n,
        p = p,
        ratio_np = n / p,
        replication = r,
        method = "Rcpp",
        time_seconds = unname(time_cpp)
      )
    )
  }

  rm(X, y)
  gc()
}

# Mean computation time across replications.
summary_results <- aggregate(
  time_seconds ~ n + p + ratio_np + method,
  data = results,
  FUN = mean
)

# Standard deviation across replications.
sd_results <- aggregate(
  time_seconds ~ n + p + ratio_np + method,
  data = results,
  FUN = sd
)

names(sd_results)[names(sd_results) == "time_seconds"] <- "sd_time"

summary_results <- merge(
  summary_results,
  sd_results,
  by = c("n", "p", "ratio_np", "method")
)

print(summary_results)

# Plot average computation time against n / p.
methods <- unique(summary_results$method)

plot(
  NA,
  xlim = range(summary_results$ratio_np),
  ylim = range(summary_results$time_seconds),
  xlab = "n / p",
  ylab = "Mean computation time (seconds)",
  log = "x"
)

for (method_name in methods) {

  tmp <- summary_results[summary_results$method == method_name, ]
  tmp <- tmp[order(tmp$ratio_np), ]

  lines(
    tmp$ratio_np,
    tmp$time_seconds,
    type = "b",
    pch = if (method_name == "R") 16 else 17
  )
}

legend(
  "topleft",
  legend = c("screening_test_matrix_R", "screening_test_matrix"),
  pch = c(16, 17),
  lty = 1,
  bty = "n"
)



df <- summary_results


df[df$method == "R",5]/df[df$method == "Rcpp",5]





##############################################################################



p <- 10000
n <- 100000
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
dim(X)
y <- rnorm(n)
system.time({result_cpp <- screening_test_matrix(X, y, q = q)})





