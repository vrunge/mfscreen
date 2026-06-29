# Compare the R and Rcpp implementations of reciprocal model-free screening.
#
# Required objects in the current R session:
#   screening_test_matrix_R : reference pure-R implementation
#   screening_test_matrix   : compiled Rcpp implementation
#
# Important: the C++ entry point must be named screening_test_matrix().
# If screening_test.cpp currently exports screening_test_matrix_R(), rename the
# exported C++ function before loading it; otherwise it will mask the R version.
#
# Both functions are expected to return:
#   $scores, $threshold, $gamma, $q, $selected, $selected_variables

# --------------------------- Configuration ---------------------------
set.seed(20260629)

n <- 2500            # observations
p <- 10000             # predictors
q <- 0.10              # screening level
n_repetitions <- 20L   # timing repetitions after one warm-up call
tol <- sqrt(.Machine$double.eps)
relative_tolerance <- 1e-10
absolute_tolerance <- 1e-12
write_csv <- FALSE

# ----------------------- Preconditions and data -----------------------
required_functions <- c("screening_test_matrix_R", "screening_test_matrix")
missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]


# A few marginally associated predictors and many null predictors.
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- sprintf("X%04d", seq_len(p))
y <- 1.25 * X[, 1L] - 0.80 * X[, 2L] + 0.45 * X[, 3L] + rnorm(n)


run_R <- function() {
  screening_test_matrix_R(X, y, q = q, tol = tol)
}

run_cpp <- function() {
  screening_test_matrix(X, y, q = q, tol = tol)
}

# -------------------------- Result comparison -------------------------
extract_scores <- function(result, implementation) {
  if (!is.list(result) || is.null(result$scores)) {
    stop(implementation, " did not return a list containing $scores.")
  }
  as.numeric(result$scores)
}

extract_selected <- function(result, implementation) {
  if (!is.list(result) || is.null(result$selected)) {
    stop(implementation, " did not return a list containing $selected.")
  }
  as.logical(result$selected)
}

numerically_equal <- function(a, b, rtol = relative_tolerance,
                              atol = absolute_tolerance) {
  if (length(a) != length(b)) {
    return(rep(FALSE, max(length(a), length(b))))
  }

  same_infinity <- is.infinite(a) & is.infinite(b) & sign(a) == sign(b)
  same_nan <- is.nan(a) & is.nan(b)
  both_finite <- is.finite(a) & is.finite(b)

  out <- same_infinity | same_nan
  out[both_finite] <- abs(a[both_finite] - b[both_finite]) <=
    atol + rtol * pmax(abs(a[both_finite]), abs(b[both_finite]))
  out[!(same_infinity | same_nan | both_finite)] <- FALSE
  out
}

r_result <- run_R()
cpp_result <- run_cpp()

r_scores <- extract_scores(r_result, "screening_test_matrix_R")
cpp_scores <- extract_scores(cpp_result, "screening_test_matrix")
r_selected <- extract_selected(r_result, "screening_test_matrix_R")
cpp_selected <- extract_selected(cpp_result, "screening_test_matrix")

if (length(r_scores) != p || length(cpp_scores) != p) {
  stop("At least one implementation returned a score vector with the wrong length.")
}

score_equal <- numerically_equal(r_scores, cpp_scores)
selected_equal <- !xor(r_selected, cpp_selected)

score_difference <- r_scores - cpp_scores
score_difference[is.infinite(r_scores) & is.infinite(cpp_scores) &
                   sign(r_scores) == sign(cpp_scores)] <- 0

comparison <- data.frame(
  variable = colnames(X),
  score_R = r_scores,
  score_cpp = cpp_scores,
  absolute_difference = abs(score_difference),
  score_equal = score_equal,
  selected_R = r_selected,
  selected_cpp = cpp_selected,
  selected_equal = selected_equal,
  row.names = NULL,
  check.names = FALSE
)

threshold_equal <- isTRUE(all.equal(
  r_result$threshold, cpp_result$threshold,
  tolerance = relative_tolerance
))
gamma_equal <- isTRUE(all.equal(
  r_result$gamma, cpp_result$gamma,
  tolerance = relative_tolerance
))
selected_variables_equal <- identical(
  as.character(r_result$selected_variables),
  as.character(cpp_result$selected_variables)
)

cat("\n=== Numerical comparison ===\n")
cat("Predictors compared:             ", p, "\n", sep = "")
cat("Matching reciprocal statistics:  ", sum(score_equal), " / ", p, "\n", sep = "")
cat("Matching selection decisions:    ", sum(selected_equal), " / ", p, "\n", sep = "")
cat("Matching threshold:              ", threshold_equal, "\n", sep = "")
cat("Matching gamma:                  ", gamma_equal, "\n", sep = "")
cat("Matching selected-variable list: ", selected_variables_equal, "\n", sep = "")

mismatch <- comparison[!comparison$score_equal | !comparison$selected_equal, ]
if (nrow(mismatch)) {
  cat("\nFirst mismatches:\n")
  print(utils::head(mismatch, 10L), row.names = FALSE)
} else {
  cat("\nAll scores and screening decisions agree within tolerance.\n")
}














# --------------------------- Timing comparison ------------------------
measure_elapsed <- function(fun) {
  unname(system.time({
    ans <- fun()
    invisible(ans)
  })[["elapsed"]])
}

# Warm-up avoids charging one-time Rcpp symbol resolution and allocation costs.
invisible(run_R())
invisible(run_cpp())

r_elapsed <- numeric(n_repetitions)
cpp_elapsed <- numeric(n_repetitions)

# Alternate the order to reduce systematic cache / CPU-frequency bias.
for (i in seq_len(n_repetitions)) {
  if (i %% 2L == 1L) {
    r_elapsed[i] <- measure_elapsed(run_R)
    cpp_elapsed[i] <- measure_elapsed(run_cpp)
  } else {
    cpp_elapsed[i] <- measure_elapsed(run_cpp)
    r_elapsed[i] <- measure_elapsed(run_R)
  }
}

summarise_times <- function(times) {
  c(
    min_seconds = min(times),
    median_seconds = median(times),
    mean_seconds = mean(times),
    max_seconds = max(times)
  )
}

timing <- rbind(
  R = summarise_times(r_elapsed),
  Rcpp = summarise_times(cpp_elapsed)
)
timing <- as.data.frame(timing)
timing$implementation <- rownames(timing)
rownames(timing) <- NULL
timing <- timing[, c("implementation", "min_seconds", "median_seconds",
                     "mean_seconds", "max_seconds")]

median_speedup <- median(r_elapsed) / median(cpp_elapsed)
mean_speedup <- mean(r_elapsed) / mean(cpp_elapsed)

cat("\n=== Execution time (seconds) ===\n")
print(timing, row.names = FALSE, digits = 6)
cat("\nMedian speed-up (R / Rcpp): ", format(median_speedup, digits = 5), "x\n", sep = "")
cat("Mean speed-up   (R / Rcpp): ", format(mean_speedup, digits = 5), "x\n", sep = "")

if (write_csv) {
  utils::write.csv(comparison, "screening_statistic_comparison.csv", row.names = FALSE)
  utils::write.csv(timing, "screening_timing_comparison.csv", row.names = FALSE)
}
