# mfscreen

`mfscreen` is an R package, powered by Rcpp, for **model-free marginal screening** in high-dimensional data.

The package identifies explanatory variables associated with a response without assuming a particular outcome model, such as a linear or logistic regression, or a specific data distribution. It can therefore be applied to continuous, binary, or discrete responses when its moment conditions are appropriate.

## What it does

Given a response vector `y` and a matrix of predictors `x`, `mfscreen` evaluates each predictor separately using a studentized slope statistic. Predictors are retained when their absolute scores exceed a threshold associated with a user-chosen false-positive rate `q`.

This provides a fast first-stage screening method for data sets with many candidate predictors.

## Installation

The development version can be installed from GitHub once the package source is available:

```r
install.packages("remotes")
remotes::install_github("vrunge/mfscreen")
```

## Typical workflow

1. Store the response in a vector `y` of length `n`.
2. Store predictors in a matrix or data frame `x`, with `n` rows and one column per candidate variable.
3. Choose a target false-positive rate `q`.
4. Run the screening procedure and retain the selected variables.

Function-level examples will be added once the public package interface is finalized.

## Method

For each predictor, the method estimates its marginal association with the response and forms a studentized score. A variable is selected when its absolute score exceeds the normal-theory threshold associated with `q`.

The method is designed to support sure screening and false-positive-rate control under the conditions described in the accompanying paper.

## Reference

Dedecker, J., Taupin, M. L., and Tocquet, A. S. (2025). *A model-free Screening procedure*.
