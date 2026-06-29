
#' Reciprocal model-free screening statistic for one predictor
#'
#' Computes the squared reciprocal screening statistic for one predictor
#' \code{x} and one response vector \code{y}. This is the distance-based
#' formulation of the model-free marginal screening statistic.
#'
#' We center and normalize the predictor,
#' \deqn{
#'   \widetilde X_i =
#'   \frac{X_i - \overline X}
#'        {\left\{\sum_{k=1}^n (X_k - \overline X)^2\right\}^{1/2}}
#' }
#' and also center the response,
#' \deqn{
#'   \widetilde Y_i =
#'   Y_i - \overline Y
#' }
#' We define
#' \deqn{
#'   A = \sum_{i=1}^n \widetilde X_i \widetilde Y_i.
#' }
#'
#' For nonzero \eqn{A}, the function returns
#' \deqn{
#'   D^2 =
#'   \sum_{i=1}^n
#'   \left[
#'     \frac{\widetilde X_i \widetilde Y_i}{A}
#'     - \widetilde X_i^2
#'   \right]^2.
#' }
#'
#' The statistic satisfies
#' \deqn{
#'   D^2 = \frac{1}{T^2}
#'        = \frac{\widehat v}{n\widehat\tau^2},
#' }
#' where \eqn{T} is the studentized screening statistic returned by
#' \code{screening_score_old_R()}.
#'
#' A predictor is selected by the reciprocal rule when
#' \deqn{
#'   D^2 \leq \frac{1}{\gamma^2},
#' }
#' where \eqn{\gamma = \code{qnorm}(1 - q / 2)} for a chosen screening
#' level \eqn{q}.
#'
#' @param x Numeric vector containing observations of one predictor.
#' @param y Numeric vector containing observations of the response.
#'   Must have the same length as \code{x}.
#' @param tol Non-negative numerical tolerance used to determine whether
#'   the normalized empirical covariance \eqn{A} is effectively zero.
#'   Defaults to \code{sqrt(.Machine$double.eps)}.
#'
#' @return A numeric scalar containing the squared reciprocal statistic
#'   \eqn{D^2}. If \eqn{A} is zero, or numerically indistinguishable from
#'   zero according to \code{tol}, the function returns \code{Inf}.
#'
#' @details
#' Non-finite values in \code{x} or \code{y} are removed pairwise before
#' computing the statistic. At least three complete observations are required.
#'
#' The function stops with an error when \code{x} and \code{y} are not
#' numeric vectors of equal length, when fewer than three complete pairs
#' remain, or when \code{x} has zero empirical variance.
#'
#' @seealso
#' \code{\link{screening_score_old_R}} for the corresponding studentized
#' screening statistic \eqn{T}.
#'
#' @examples
#'
#' n <- 2500
#' q <- 0.05
#'
#' # One predictor with a marginal association with y
#' x_signal <- rnorm(n)
#'
#' # One predictor independent of y
#' x_null <- rnorm(n)
#'
#' # Response generated from x_signal only
#' y <- x_signal^2 + 0.5 * x_signal + rnorm(n, sd = 1)
#'
#' # Reciprocal distance statistic
#' D2_signal <- screening_scoreR(x_signal, y)
#' D2_null <- screening_scoreR(x_null, y)
#'
#' # screening decisions
#' gamma <- qnorm(1 - q / 2)
#'
#' selected <- c(
#'   signal = D2_signal <= 1 / gamma^2,
#'   null = D2_null <= 1 / gamma^2
#' )
#'
#' data.frame(
#'   variable = c("signal", "null"),
#'   D2 = c(D2_signal, D2_null),
#'   selected = selected
#' )
#'
#' @export
screening_scoreR <- function(
    x,
    y,
    tol = sqrt(.Machine$double.eps)
)
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

  # Remove missing and non-finite observations pairwise
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

  # Center x and y
  xc <- x - mean(x)
  yc <- y - mean(y)

  # Normalize centered x so that sum(x_tilde^2) = 1
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

  x_tilde <- xc / sqrt(sxx)

  # A_j = sum_i x_tilde_i * y_tilde_i
  A <- sum(x_tilde * yc)

  # By definition: if A_j = 0, D_j^2 = +Inf
  A_scale <- max(1, sqrt(sum(x_tilde^2 * yc^2)))

  if (!is.finite(A))
  {
    return(Inf)
  }

  ##
  #####
  #########################################
  #####
  ##

  # p_j and q_j
  p <- (x_tilde * yc) / A
  q <- x_tilde^2

  # D_j^2 = ||p_j - q_j||_2^2
  D2 <- sum((p - q)^2)

  return(D2)
}



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



screening_test_matrix_R <- function(
    X,
    y,
    q = 0.10,
    tol = sqrt(.Machine$double.eps)
)
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

  # Validate tolerance
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol < 0)
  {
    stop("tol must be a non-negative finite number.")
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
  # Two-sided normal threshold
  gamma <- qnorm(1 - q / 2)

  # Reciprocal screening threshold:
  # D_j^2 <= 1 / gamma^2
  threshold <- 1 / gamma^2

  # Compute one squared reciprocal statistic per predictor
  D2 <- vapply(
    X,
    FUN = function(xj) screening_scoreR(xj, y, tol = tol),
    FUN.VALUE = numeric(1)
  )

  # Screening decision
  selected <- is.finite(D2) & D2 <= threshold
  selected_variables <- names(D2)[selected]

  list(
    scores = D2,
    threshold = threshold,
    gamma = gamma,
    q = q,
    selected = selected,
    selected_variables = selected_variables
  )
}
