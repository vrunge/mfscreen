test_that("large-offset predictors use the stable path", {
  set.seed(321)
  n <- 120
  p <- 12
  X <- matrix(1e9 + rnorm(n * p), nrow = n, ncol = p)
  y <- rnorm(n)

  result <- screening_test_matrix(X, y, q = 0.10)
  expected_scores <- vapply(
    seq_len(p),
    function(j) screening_score(X[, j], y),
    numeric(1)
  )

  expect_true(all(is.finite(result$scores)))
  expect_equal(result$scores, expected_scores, tolerance = 1e-9)
  expect_equal(result$selected_indices, which(result$selected))
})
