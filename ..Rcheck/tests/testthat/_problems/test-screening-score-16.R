# Extracted from test-screening-score.R:16

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "mfscreen", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
y <- seq_len(20)
expect_equal(screening_score(rep(1, 20), y), Inf)
