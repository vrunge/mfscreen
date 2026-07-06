
set.seed(123)

n <- 100
p <- 100000
q <- 0.15

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- 3 * X[, 1] - 2 * X[, 2] + rnorm(n)

# Threshold for screening_score_old_R()
gamma <- qnorm(1 - q / 2)

# ------------------------------------------------------------
# 1. screening_score_old_R(), called once per column
# ------------------------------------------------------------
time_old_R <- system.time({
  scores_old_R <- vapply(
    seq_len(p),
    function(j) screening_score_old_R(X[, j], y),
    numeric(1)
  )

  selected_old_R <- which(abs(scores_old_R) >= gamma)
})["elapsed"]

# ------------------------------------------------------------
# 2. Matrix-level reciprocal R implementation
# ------------------------------------------------------------
time_matrix_R <- system.time({
  result_matrix_R <- screening_test_matrix_R(X, y, q = q)
})["elapsed"]

# ------------------------------------------------------------
# 3. Matrix-level Rcpp implementation
# ------------------------------------------------------------
time_cpp <- system.time({
  result_cpp <- screening_test_matrix(X, y, q = q)
})["elapsed"]

# ------------------------------------------------------------
# Results
# ------------------------------------------------------------
cat("\nExecution times\n")
cat("-----------------------------\n")
cat("screening_score_old_R :", time_old_R, "seconds\n")
cat("screening_test_matrix_R:", time_matrix_R, "seconds\n")
cat("screening_test_matrix :", time_cpp, "seconds\n")

cat("\nSpeed-up relative to screening_score_old_R\n")
cat("------------------------------------------\n")

cat(
  "screening_test_matrix_R:",
  round(time_old_R / time_matrix_R, 2),
  "x\n"
)

cat(
  "screening_test_matrix:",
  round(time_old_R / time_cpp, 2),
  "x\n"
)


cat(
  "screening_test_matrix:",
  round(time_old_R / time_cpp, 2),
  "x\n")






################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



cat("\n Selected-variable checks \n ")
cat("----------------------------- \n ")

# These should agree because T_j^2 = 1 / D_j^2,
# apart from rare numerical boundary differences.
cat(
  "old R vs matrix R:",
  identical(
    selected_old_R,
    unname(which(result_matrix_R$selected))),"\n")

cat(
  "old R vs C++:",
  identical(
    selected_old_R,
    result_cpp$selected_indices),"\n")

cat(
  "matrix R vs C++:",
  identical(
    unname(which(result_matrix_R$selected)),
    result_cpp$selected_indices
  ),"\n")

