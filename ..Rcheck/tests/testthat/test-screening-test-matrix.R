test_that("screening_test_matrix returns the expected structure", {
  set.seed(123)
  n <- 120
  p <- 25
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- 2 * X[, 1] - 1.25 * X[, 2] + rnorm(n)

  result <- screening_test_matrix(X, y, q = 0.10)

  expect_named(
    result,
    c("scores", "threshold", "gamma", "q", "selected",
      "selected_indices")
  )
  expect_length(result$scores, p)
  expect_length(result$selected, p)
  expect_type(result$selected, "logical")
  expect_type(result$selected_indices, "integer")
  expect_equal(result$threshold, 1 / result$gamma^2)
  expect_equal(result$q, 0.10)
  expect_equal(result$selected_indices, which(result$selected))
})

test_that("screening_test_matrix agrees with columnwise screening_score", {
  set.seed(456)
  n <- 150
  p <- 40
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, 5] <- 1
  y <- X[, 1] - 2 * X[, 3] + rnorm(n)
  q <- 0.15

  result <- screening_test_matrix(X, y, q = q)
  expected_scores <- vapply(
    seq_len(p),
    function(j) screening_score(X[, j], y),
    numeric(1)
  )
  threshold <- 1 / qnorm(1 - q / 2)^2

  expect_equal(result$scores, expected_scores, tolerance = 1e-9)
  expect_equal(result$threshold, threshold)
  expect_equal(result$selected, is.finite(expected_scores) &
                 expected_scores <= threshold)
  expect_equal(result$selected_indices, which(result$selected))
})

test_that("screening_test_matrix agrees with the R reference selection", {
  set.seed(789)
  n <- 100
  p <- 30
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- 1.2 * X[, 2] - 0.8 * X[, 4] + rnorm(n)

  cpp <- screening_test_matrix(X, y, q = 0.10)
  r_ref <- screening_test_matrix_R(X, y, q = 0.10)

  expect_equal(cpp$scores, unname(r_ref$scores), tolerance = 1e-9)
  expect_equal(cpp$selected_indices, unname(which(r_ref$selected)))
})
