mfscreen

mfscreen is an R package, powered by Rcpp, for model-free marginal screening in high-dimensional data.

It identifies predictors associated with a response without requiring a fully specified outcome model, such as linear regression, logistic regression, or a particular response distribution. The procedure is based on marginal moment conditions and can be used with continuous, binary, or discrete responses when its assumptions are appropriate.

What it does

Given a response vector y and a predictor matrix X, mfscreen evaluates each predictor separately.

For each predictor, the package computes a reciprocal screening statistic,

[
D_j^2 = \frac{1}{T_j^2},
]

where (T_j) is the corresponding studentized marginal association statistic.

Smaller values of (D_j^2) indicate stronger marginal association between predictor (X_j) and the response. A predictor is retained when

[
D_j^2 \leq \frac{1}{\gamma^2},
\qquad
\gamma = \Phi^{-1}(1-q/2),
]

where q is the target two-sided false-positive rate.

The package is implemented in C++ through Rcpp and is designed for fast first-stage screening when the number of candidate predictors is large.

Installation

Install the development version from GitHub:

install.packages("remotes")
remotes::install_github("vrunge/mfscreen")

Typical workflow

1. Store the response in a numeric vector y of length n.
2. Store predictors in a numeric matrix X with n rows and one column per candidate variable.
3. Choose a target marginal false-positive rate q.
4. Run screening_test_matrix().
5. Extract the retained variables using selected or selected_indices.

Example

library(mfscreen)
set.seed(123)
n <- 300
p <- 1000
# Candidate predictors
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
# Response associated with the first three predictors
y <- 2 * X[, 1] - 1.5 * X[, 2] + X[, 3] + rnorm(n)
# Model-free marginal screening
result <- screening_test_matrix(
  X = X,
  y = y,
  q = 0.10
)
# One-based indices of selected predictors
result$selected_indices
# Logical selection indicator for each predictor
result$selected
# Reciprocal screening scores D_j^2:
# smaller values indicate stronger marginal association
head(result$scores)
# Selection threshold and corresponding normal threshold
result$threshold
result$gamma

In this example, X1, X2, and X3 will typically have small reciprocal screening scores and are likely to be retained. Independent noise variables will generally have larger scores and will not be selected.

Output

screening_test_matrix() returns a list with:

* scores: reciprocal screening scores (D_j^2), one per predictor;
* threshold: reciprocal selection threshold (1/\gamma^2);
* gamma: normal-theory threshold (\Phi^{-1}(1-q/2));
* q: supplied false-positive rate;
* selected: logical vector indicating retained predictors;
* selected_indices: one-based indices of retained predictors.

Single-predictor score

For a single predictor, use screening_score():

score_x1 <- screening_score(X[, 1], y)
score_x1

This returns the reciprocal score (D^2). Smaller values indicate stronger marginal association.

Input requirements

The C++ functions assume that:

* X and y contain finite numeric values;
* length(y) == nrow(X);
* q is strictly between 0 and 1;
* tol is non-negative.

Predictors with zero empirical variance, or with empirical covariance numerically indistinguishable from zero, receive score Inf and are not selected.

Method

For each predictor, mfscreen centers the predictor and response, computes a marginal association statistic, and applies a heteroskedasticity-robust normalization. The resulting reciprocal score is algebraically equivalent to the inverse squared studentized marginal statistic.

The procedure is intended as a computationally efficient screening stage before more detailed modelling or variable-selection methods.

Reference

Dedecker, J., Taupin, M. L., and Tocquet, A. S. (2025). A model-free Screening procedure.
