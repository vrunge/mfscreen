test_that("screening_score agrees with the R reference implementation", {
  set.seed(42)
  x <- rnorm(200)
  y <- 1.5 * x + 0.3 * x^2 + rnorm(200)

  expect_equal(
    screening_score(x, y),
    screening_score_R(x, y),
    tolerance = 1e-10
  )
})

test_that("screening_score returns Inf for constant and uncorrelated inputs", {
  y <- as.numeric(seq_len(20))

  expect_equal(screening_score(rep(1.0, 20), y), Inf)
  expect_equal(screening_score(c(-1, 1, -1, 1), rep(2.0, 4)), Inf)
})
