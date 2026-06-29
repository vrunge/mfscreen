## =========================================================
## Assumes these two functions are already defined:
##   - screening_scoreR(x, y)
##   - screening_test_matrixR(X, y, q = 0.10)
## =========================================================

## Helper: Gaussian sample with covariance R'R = Sigma
rmvnorm_chol <- function(n, R)
{
  Z <- matrix(rnorm(n * ncol(R)), nrow = n, ncol = ncol(R))
  Z %*% R
}

##
##
## Build the explanatory variables exactly as in the paper
##
##

make_X_continuous_paper <- local({
  ## First block: X1,...,X5 with rho = 0.57
  rho_5 <- 0.57
  Sigma_5 <- toeplitz(rho_5^(0:4))
  R_5 <- chol(Sigma_5)

  ## Last block: X19,...,X2000 stationary Gaussian with rho = 0.7
  rho_tail <- 0.70
  Sigma_tail <- toeplitz(rho_tail^(0:1981))
  R_tail <- chol(Sigma_tail)

  function(n) {
    ## X1,...,X5
    X1to5 <- rmvnorm_chol(n, R_5)
    colnames(X1to5) <- paste0("X", 1:5)

    ## X6,...,X10
    X6  <- rbinom(n, size = 1, prob = 0.35)
    X7  <- rchisq(n, df = 2)
    Z1  <- rpois(n, lambda = 2)
    Z2  <- rnorm(n, mean = 1, sd = 1)
    X8  <- X1to5[, 1] * Z1
    X9  <- X1to5[, 2] * Z2
    X10 <- rt(n, df = 14)

    ## Noise variables for X11,...,X18
    ## In the paper N(0, a) is mathematical notation, so I use variance = a
    Z3 <- rnorm(n, mean = 0, sd = sqrt(0.4))
    Z4 <- rnorm(n, mean = 0, sd = sqrt(0.02))
    Z5 <- rnorm(n, mean = 0, sd = sqrt(0.35))
    Z6 <- rnorm(n, mean = 0, sd = sqrt(0.55))

    ## X11,...,X14
    X11 <- X1to5[, 1] + Z3
    X12 <- X1to5[, 3] + Z4
    X13 <- X1to5[, 4] + Z5
    X14 <- X1to5[, 5] + Z6

    ## X15,...,X18
    X15 <- Z3
    X16 <- Z4
    X17 <- Z5
    X18 <- Z6

    ## X19,...,X2000
    X19to2000 <- rmvnorm_chol(n, R_tail)
    colnames(X19to2000) <- paste0("X", 19:2000)

    X <- cbind(
      X1to5,
      X6 = X6, X7 = X7, X8 = X8, X9 = X9, X10 = X10,
      X11 = X11, X12 = X12, X13 = X13, X14 = X14,
      X15 = X15, X16 = X16, X17 = X17, X18 = X18,
      X19to2000
    )

    as.data.frame(X)
  }
})

##
##
## Continuous response from the paper
##
##
make_Y_continuous_paper <- function(X)
{
  theta <- c(
    -1, 3.619350, -3.274923, 2.963273, -2.681280,
    2, 4, 6, 3, 2, 4
  )

  eta <- theta[1] +
    theta[2]  * X$X1  +
    theta[3]  * X$X2  +
    theta[4]  * X$X3  +
    theta[5]  * X$X4  +
    theta[6]  * X$X5  +
    theta[7]  * X$X6  +
    theta[8]  * X$X7  +
    theta[9]  * X$X8  +
    theta[10] * X$X9  +
    theta[11] * X$X10

  eps <- rchisq(nrow(X), df = 3) - 3
  abs(eta)^0.8 + X$X5 * eps / 3
}

##
##
## One replicate of the continuous simulation
##
##
one_rep_continuous <- function(n = 500, q = 0.10)
{
  X <- make_X_continuous_paper(n)
  y <- make_Y_continuous_paper(X)

  out <- screening_test_matrixR(X, y, q = q)
  selected <- as.logical(out$results$selected)

  list(
    X = X,
    y = y,
    selected = selected,
    scores = out$results$score,
    score_abs = out$results$score_abs,
    threshold = out$threshold
  )
}

##
##
## Monte Carlo reproduction of Table 1 for one sample size n
##
##
reproduce_table1_continuous <- function(n = 500, N = 500, q = 0.10, seed = 1)
{
  set.seed(seed)

  p <- 2000
  selected_mat <- matrix(FALSE, nrow = N, ncol = p)
  colnames(selected_mat) <- paste0("X", 1:p)

  for (b in seq_len(N))
  {
    sim_b <- one_rep_continuous(n = n, q = q)
    selected_mat[b, ] <- sim_b$selected
  }

  selection_freq <- colMeans(selected_mat)

  per_variable <- data.frame(
    variable = paste0("X", 1:18),
    selection_frequency = selection_freq[1:18]
  )

  list(
    per_variable = per_variable,
    TPR = mean(selected_mat[, 1:14]),         # variables correlated with Y
    FPR = mean(selected_mat[, 15:2000]),      # variables not correlated with Y
    mean_Mhat = mean(rowSums(selected_mat)),
    selected_matrix = selected_mat
  )
}

###########################################################################

res_small <- reproduce_table1_continuous(n = 500, N = 500, q = 0.10)

res_small$per_variable
res_small$TPR
res_small$FPR
res_small$mean_Mhat



