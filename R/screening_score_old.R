
#' Model-free screening statistic for one predictor
#'
#' Computes the studentized marginal screening statistic for one predictor
#' \code{x} and response vector \code{y}. The statistic is based on the
#' marginal ordinary least-squares slope and a heteroskedasticity-robust
#' variance estimator.
#'
#' Let \eqn{\widehat{\tau}} denote the marginal OLS slope obtained by
#' regressing \eqn{Y} on \eqn{X}. The function returns
#' \deqn{
#'   T =
#'   \frac{\sqrt{n}\widehat{\tau}}{\sqrt{\widehat{v}}},
#' }
#' where
#' \deqn{
#'   \widehat{v} =
#'   \frac{
#'     n^{-1} \sum_{i=1}^n
#'     (X_i - \overline{X})^2 \widehat{\varepsilon}_i^2
#'   }{
#'     \left\{
#'       n^{-1} \sum_{i=1}^n
#'       (X_i - \overline{X})^2
#'     \right\}^2
#'   }.
#' }
#'
#' Here, \eqn{\widehat{\varepsilon}_i} is the residual from the marginal
#' regression of \eqn{Y} on \eqn{X}. Large absolute values of the returned
#' statistic indicate a marginal association between the predictor and the
#' response.
#'
#' @param x Numeric vector containing observations of one predictor.
#' @param y Numeric vector containing observations of the response. Must have
#'   the same length as \code{x}.
#'
#' @return A numeric scalar containing the signed studentized screening score.
#'
#' @details
#' Non-finite values in \code{x} or \code{y} are removed pairwise before the
#' statistic is computed. At least three complete observations are required.
#'
#' The function stops when the inputs are not numeric vectors of equal length,
#' fewer than three complete pairs remain, \code{x} has zero empirical
#' variance, or the robust variance estimate is not strictly positive.
#'
#' @seealso
#' \code{\link{screening_test_matrix_old_R}} for applying the statistic to all
#' columns of a predictor matrix.
#'
#' @examples
#' set.seed(123)
#'
#' x_signal <- rnorm(200)
#' x_null <- rnorm(200)
#' y <- 2 * x_signal + rnorm(200)
#'
#' screening_score_old_R(x_signal, y)
#' screening_score_old_R(x_null, y)
#'
#' @export
screening_score_old_R <- function(x, y)
{
  ##
  #####
  #########################################
  #####
  ##

  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y))
  {
    stop("x and y must be numeric vectors of equal length.")
  }

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  n <- length(x)
  if (n < 3L)
  {
    stop("Need at least 3 complete observations.")
  }

  ##
  #####
  #########################################
  #####
  ##

  x_bar <- mean(x)
  y_bar <- mean(y)

  xc <- x - x_bar
  yc <- y - y_bar

  sxx <- sum(xc^2)
  if (!is.finite(sxx) || sxx <= 0)
  {
    stop("x has zero empirical variance.")
  }

  ##
  #####
  #########################################
  #####
  ##
  # Marginal OLS slope and intercept
  tau_hat <- sum(xc * yc) / sxx

  # Equivalent to y - b_hat - tau_hat * x
  eps_hat <- yc - tau_hat * xc

  # Paper's heteroskedasticity-robust estimator
  v_hat <- mean(xc^2 * eps_hat^2) / mean(xc^2)^2

  if (!is.finite(v_hat) || v_hat <= 0)
  {
    stop("Estimated variance is not positive.")
  }

  score <- sqrt(n) * tau_hat / sqrt(v_hat)

  return(score)
}





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


#' Model-free marginal screening for a predictor matrix
#'
#' Applies \code{screening_score_old_R()} separately to every column of
#' \code{X} and selects predictors whose absolute studentized screening score
#' exceeds a two-sided normal threshold.
#'
#' For predictor \eqn{j}, let \eqn{T_j} be its studentized marginal screening
#' statistic. Predictor \eqn{j} is selected when
#' \deqn{
#'   |T_j| \geq \Phi^{-1}(1-q/2),
#' }
#' where \eqn{q} is the target marginal false-positive rate and \eqn{\Phi}
#' is the standard normal distribution function.
#'
#' @param X Numeric matrix or data frame of predictors. Rows correspond to
#'   observations and columns correspond to candidate predictors. A numeric
#'   vector is accepted and treated as a one-column matrix.
#' @param y Numeric response vector with one entry per row of \code{X}.
#' @param q Numeric scalar in \code{(0, 1)} giving the target two-sided
#'   false-positive rate used to form the normal threshold.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{scores}{Named numeric vector of studentized screening scores, one
#'   per predictor.}
#'   \item{threshold}{Two-sided normal threshold
#'   \eqn{\Phi^{-1}(1-q/2)}.}
#'   \item{q}{The supplied false-positive rate.}
#'   \item{selected}{Named logical vector indicating selected predictors.}
#'   \item{selected_variables}{Character vector containing the names of the
#'   selected predictors.}
#' }
#'
#' @details
#' The score for each predictor is computed independently using pairwise
#' complete observations for that predictor and \code{y}. Consequently, if
#' columns of \code{X} have different missing-value patterns, their scores may
#' be based on different subsets of observations.
#'
#' When \code{X} has no column names, the function creates names
#' \code{"X1"}, \code{"X2"}, and so on.
#'
#' @seealso
#' \code{\link{screening_score_old_R}} for the single-predictor statistic.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 2000
#' p <- 100
#'
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' colnames(X) <- paste0("X", seq_len(p))
#'
#' y <- 2 * X[, 1] - 1.5 * X[, 2] + rnorm(n)
#'
#' result <- screening_test_matrix_old_R(X, y, q = 0.10)
#'
#' result$selected_variables
#' head(result$scores)
#' result$threshold
#'
#' @export
screening_test_matrix_old_R <- function(X, y, q = 0.10)
{
  ##
  #####
  #########################################
  #####
  ##
  # Validate q
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) ||
      q <= 0 || q >= 1)
  {
    stop("q must be a finite number in (0, 1).")
  }

  # Allow one predictor supplied as a vector
  if (is.vector(X))
  {
    X <- matrix(X, ncol = 1L)
  }

  X <- as.data.frame(X, check.names = FALSE)

  # Validate input dimensions and types
  if (!is.numeric(y))
  {
    stop("y must be a numeric vector.")
  }

  if (nrow(X) != length(y))
  {
    stop("X and y must have the same number of observations.")
  }

  if (!all(vapply(X, is.numeric, logical(1))))
  {
    stop("All columns of X must be numeric.")
  }

  # Create names when X has no column names
  if (is.null(names(X)))
  {
    names(X) <- paste0("X", seq_len(ncol(X)))
  }

  ##
  #####
  #########################################
  #####
  ##
  # Two-sided normal threshold used by the screening procedure
  threshold <- qnorm(1 - q / 2)

  # Compute one studentized screening score per predictor
  scores <- vapply(
    X,
    FUN = function(xj) screening_score_old_R(xj, y),
    FUN.VALUE = numeric(1)
  )

  # Screening decision
  selected <- abs(scores) >= threshold
  selected_variables <- names(scores)[selected]

  list(
    scores = scores,
    threshold = threshold,
    q = q,
    selected = selected,
    selected_variables = selected_variables
  )
}

