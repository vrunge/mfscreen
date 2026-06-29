## =========================================================
## Simple Gaussian simulation for the screening procedure
## =========================================================

##
##
## One simulated dataset
##
##
simulate_simple_gaussian_data <- function(n, p = 200)
{
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  ## True active variables: X1,...,X5
  beta <- c(2.0, -1.5, 1.2, 1.0, -0.8)

  ## Response built only from the first 5 variables
  y <- X[, 1] * beta[1] +
    X[, 2] * beta[2] +
    X[, 3] * beta[3] +
    X[, 4] * beta[4] +
    X[, 5] * beta[5] +
    rnorm(n, mean = 0, sd = 2)

  list(
    X = as.data.frame(X),
    y = y,
    active_set = 1:5
  )
}

##
##
## One Monte Carlo replicate
##
##
one_rep_simple_gaussian <- function(n, p = 200, q = 0.10)
{
  dat <- simulate_simple_gaussian_data(n = n, p = p)
  out <- screening_test_matrixR(dat$X, dat$y, q = q)

  selected <- as.logical(out$results$selected)

  list(
    X = dat$X,
    y = dat$y,
    active_set = dat$active_set,
    selected = selected,
    results = out$results,
    threshold = out$threshold
  )
}

##
##
## Reproduce a simple screening table for one sample size
##
##
reproduce_simple_gaussian <- function(n = 200, p = 200, N = 100, q = 0.10, seed = 123)
{
  set.seed(seed)

  active_set <- 1:5
  inactive_set <- setdiff(1:p, active_set)

  selected_mat <- matrix(FALSE, nrow = N, ncol = p)
  colnames(selected_mat) <- paste0("X", 1:p)

  pb <- txtProgressBar(min = 0, max = N, style = 3)

  for (b in seq_len(N))
  {
    sim_b <- one_rep_simple_gaussian(n = n, p = p, q = q)
    selected_mat[b, ] <- sim_b$selected
    setTxtProgressBar(pb, b)
  }

  close(pb)

  selection_freq <- colMeans(selected_mat)

  per_variable <- data.frame(
    variable = paste0("X", 1:p),
    selection_frequency = selection_freq
  )

  list(
    per_variable = per_variable,
    TPR = mean(selected_mat[, active_set]),
    FPR = mean(selected_mat[, inactive_set]),
    mean_Mhat = mean(rowSums(selected_mat)),
    selected_matrix = selected_mat,
    active_variables = paste0("X", active_set),
    threshold = qnorm(1 - q / 2),
    q = q,
    n = n,
    p = p,
    N = N
  )
}

################################################################################

res <- reproduce_simple_gaussian(n = 2500, p = 2000, N = 100,
                                 q = 0.10, seed = 123)
plot(res$per_variable[,2], type = 'l')
res$per_variable[1:10,2]
res$TPR
res$FPR
res$mean_Mhat

