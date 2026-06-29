
# ------------------------------------------------------------
# Example data for screening_test_matrix_old_R()
# ------------------------------------------------------------

n <- 2500
p <- 2000

# Five correlated Gaussian predictors
rho <- 0.50
Sigma5 <- toeplitz(rho^(0:4))

X_signal <- MASS::mvrnorm(
  n = n,
  mu = rep(0, 5),
  Sigma = Sigma5
)

# Additional signal variables
X6 <- rbinom(n, size = 1, prob = 0.35)
X7 <- rchisq(n, df = 2)
X8 <- X_signal[, 1] * rpois(n, lambda = 2)
X9 <- X_signal[, 2] * rnorm(n, mean = 1, sd = 1)
X10 <- rt(n, df = 14)

# Correlated proxy variables: associated with Y through X1, X3, X4, X5
Z3 <- rnorm(n, mean = 0, sd = 0.04)
Z4 <- rnorm(n, mean = 0, sd = 0.02)
Z5 <- rnorm(n, mean = 0, sd = 0.35)
Z6 <- rnorm(n, mean = 0, sd = 0.55)

X11 <- X_signal[, 1] + Z3
X12 <- X_signal[, 3] + Z4
X13 <- X_signal[, 4] + Z5
X14 <- X_signal[, 5] + Z6

# Noise terms: should not be selected
X15 <- Z3
X16 <- Z4
X17 <- Z5
X18 <- Z6

# Independent null variables
X_noise <- matrix(rnorm(n * (p - 18L)), nrow = n, ncol = p - 18L)


# Complete design matrix
Xcomplet <- cbind(
  X_signal,
  X6, X7, X8, X9, X10,
  X11, X12, X13, X14,
  X15, X16, X17, X18,
  X_noise
)

colnames(Xcomplet) <- c(
  paste0("X", 1:5),
  "X6_bern",
  "X7_chisq2",
  "X8_X1_times_Pois2",
  "X9_X2_times_Norm",
  "X10_t14",
  "X11_X1_plus_Z3",
  "X12_X3_plus_Z4",
  "X13_X4_plus_Z5",
  "X14_X5_plus_Z6",
  "X15_Z3",
  "X16_Z4",
  "X17_Z5",
  "X18_Z6",
  paste0("X_noise_", seq_len(p - 18L))
)

dim(Xcomplet)

################################################################################

# Response: zero-inflated Poisson-like response
theta <- c(
  -1.0,   # intercept
  1.2,  -1.0, 0.9, -0.8, 0.7,
  0.9,   0.7, 0.6, -0.5, 0.7
)

linear_predictor <- as.vector(
  theta[1] + Xcomplet[, 1:10] %*% theta[-1]
)

prob_nonzero <- plogis(linear_predictor)

Y1 <- rbinom(n, size = 1, prob = prob_nonzero)
Y2 <- rpois(n, lambda = abs(linear_predictor))

Y <- Y1 * Y2

################################################################################


# ------------------------------------------------------------
# Run the screening procedure
# ------------------------------------------------------------

q_val <- 0.1

fit <- screening_test_matrix_old_R(
  X = Xcomplet,
  y = as.vector(Y),
  q = q_val
)

fit$scores
fit$threshold
fit$selected_variables



################################################################################

# ------------------------------------------------------------
# Truth definition for this simulation
# ------------------------------------------------------------

# Variables directly used to generate Y
true_active_variables <- colnames(Xcomplet)[1:10]

# Optional: identify proxies separately
proxy_variables <- colnames(Xcomplet)[11:14]

# ------------------------------------------------------------
# Detailed results table
# ------------------------------------------------------------

results_table <- data.frame(
  variable = names(fit$scores),
  score = fit$scores,
  abs_score = abs(fit$scores),
  selected = fit$selected,
  true_active = names(fit$scores) %in% true_active_variables,
  variable_type = ifelse(
    names(fit$scores) %in% true_active_variables,
    "true_signal",
    ifelse(
      names(fit$scores) %in% proxy_variables,
      "proxy_not_directly_active",
      "null"
    )
  ),
  row.names = NULL
)

# Classification of each screening decision
results_table$test_result <- with(
  results_table,
  ifelse(
    selected & true_active, "TP",   # true positive
    ifelse(
      !selected & true_active, "FN", # false negative
      ifelse(
        selected & !true_active, "FP", # false positive
        "TN"                           # true negative
      )
    )
  )
)

# Was the decision correct for this variable?
results_table$correct_decision <- results_table$test_result %in% c("TP", "TN")

# Sort by absolute score, from strongest to weakest
results_table <- results_table[
  order(-results_table$abs_score),
]

print(results_table)


# ------------------------------------------------------------
# Summary: discoveries, FDP, power, error counts
# ------------------------------------------------------------


TP <- sum(results_table$test_result == "TP")
FP <- sum(results_table$test_result == "FP")
TN <- sum(results_table$test_result == "TN")
FN <- sum(results_table$test_result == "FN")

R <- TP + FP

# False positive rate among truly null variables
FPR <- if ((FP + TN) == 0L) NA_real_ else FP / (FP + TN)

power <- TP / length(true_active_variables)

screening_summary <- data.frame(
  q_target = q_val,
  threshold = fit$threshold,
  total_variables = nrow(results_table),
  true_active_variables = length(true_active_variables),
  selected_variables = R,
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  FPR_observed = FPR,
  power = power
)

print(screening_summary)
