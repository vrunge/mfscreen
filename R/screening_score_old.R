
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

