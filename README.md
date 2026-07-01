# mfscreen

`mfscreen` is an R package, powered by Rcpp, for **model-free marginal screening** in high-dimensional data.

It evaluates each candidate predictor separately and retains variables showing evidence of marginal association with a response. The procedure does not require fitting a full linear, logistic, or other parametric outcome model.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("vrunge/mfscreen")
```

## Main function

```r
screening_test_matrix(X, y, q = 0.10, tol = sqrt(.Machine$double.eps))
```

where:

- `X` is a numeric matrix with one row per observation and one column per candidate predictor;
- `y` is a numeric response vector with `length(y) == nrow(X)`;
- `q` is the target two-sided marginal false-positive rate;
- `tol` is a numerical tolerance for treating an empirical covariance as zero.

The function returns:

- `scores`: reciprocal screening scores \(D_j^2\);
- `threshold`: the reciprocal selection threshold \(1 / \gamma^2\);
- `gamma`: the normal-theory threshold \(\gamma = \Phi^{-1}(1-q/2)\);
- `q`: the supplied false-positive rate;
- `selected`: logical selection indicators;
- `selected_indices`: one-based indices of selected predictors.

## Screening rule

For predictor \(X_j\), the package computes a reciprocal screening statistic

\[
D_j^2 = \frac{1}{T_j^2},
\]

where \(T_j\) is the corresponding studentized marginal screening statistic.

A predictor is selected when

\[
D_j^2 \leq \frac{1}{\gamma^2},
\qquad
\gamma = \Phi^{-1}(1-q/2).
\]

Therefore:

- **smaller** values of `scores` indicate stronger marginal association;
- a predictor is selected when its score is **below** `threshold`.

Predictors with zero empirical variance, or with empirical covariance numerically indistinguishable from zero, receive score `Inf` and are not selected.

## Example

```r
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

# One-based indices of retained predictors
result$selected_indices

# Selection indicator for every predictor
result$selected

# Reciprocal scores: smaller means stronger association
head(result$scores)

# Selection threshold
result$threshold

# Corresponding normal-theory threshold
result$gamma
```

In this example, predictors `X1`, `X2`, and `X3` will generally have small scores and are likely to be selected. Most independent noise predictors should have larger scores and should not be selected.

## Single-predictor score

Use `screening_score()` to compute the reciprocal statistic for one predictor:

```r
score_x1 <- screening_score(X[, 1], y)

score_x1
```

The return value is \(D^2\). Smaller values indicate stronger marginal association between the predictor and the response.

## Input requirements

The Rcpp functions expect:

- finite numeric values in `X` and `y`;
- `length(y) == nrow(X)`;
- `q` strictly between `0` and `1`;
- `tol` greater than or equal to `0`.

The low-level C++ implementation does not perform full input validation. Validate inputs before calling the exported functions when working outside the usual R workflow.

## Method

For each predictor, `mfscreen` centers the predictor and response, computes marginal empirical moments, and forms a reciprocal statistic that is algebraically equivalent to the inverse squared studentized screening statistic.

The package is designed as a fast first-stage screening method before more detailed modeling, regularization, or variable-selection procedures.

## Reference

Dedecker, J., Taupin, M. L., and Tocquet, A. S. (2025). *A model-free Screening procedure*.
