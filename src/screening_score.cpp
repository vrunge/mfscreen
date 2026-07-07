#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "screening_score.h"

using namespace Rcpp;
//' Reciprocal model-free screening statistic for one predictor
//'
//' Computes the squared reciprocal model-free screening statistic for one
//' predictor \code{x} and one response vector \code{y}.
//'
//' Let
//' \deqn{
//'   z_i = x_i - \overline{x},
//'   \qquad
//'   u_i = y_i - \overline{y}.
//' }
//'
//' Define
//' \deqn{
//'   b = \sum_{i=1}^n z_i u_i,
//'   \qquad
//'   s_{xx} = \sum_{i=1}^n z_i^2.
//' }
//'
//' The returned statistic is
//' \deqn{
//'   D^2 =
//'   \sum_{i=1}^n
//'   \left\{
//'     \frac{z_i u_i}{b}
//'     -
//'     \frac{z_i^2}{s_{xx}}
//'   \right\}^2.
//' }
//'
//' It is the squared reciprocal of the studentized marginal screening score:
//' \deqn{
//'   D^2 = \frac{1}{T^2}.
//' }
//'
//' Smaller values of \eqn{D^2} indicate stronger marginal association between
//' \code{x} and \code{y}. A value of \code{Inf} indicates that the empirical
//' covariance between the centered predictor and response is effectively zero.
//'
//' @param x Numeric vector containing one predictor.
//' @param y Numeric response vector with the same length as \code{x}.
//' @param tol Non-negative numerical tolerance used to determine whether the
//'   empirical covariance is effectively zero. Defaults to
//'   \code{sqrt(.Machine$double.eps)}.
//'
//' @return A numeric scalar containing the reciprocal screening statistic
//'   \eqn{D^2}. Returns \code{Inf} when the predictor has zero empirical
//'   variance or its empirical covariance with \code{y} is below the
//'   tolerance threshold.
//'
//' @details
//' This is a low-level Rcpp implementation. It assumes that \code{x} and
//' \code{y} are finite numeric vectors of equal positive length and that
//' \code{tol} is non-negative. No input validation is performed in C++.
//'
//' The covariance is treated as numerically zero when
//' \deqn{
//'   |b| \leq
//'   \mathrm{tol}
//'   \max\left\{
//'     \sqrt{s_{xx}},
//'     \sqrt{\sum_{i=1}^n z_i^2 u_i^2}
//'   \right\}.
//' }
//'
//' @examples
//' set.seed(123)
//'
//' x_signal <- rnorm(200)
//' x_null <- rnorm(200)
//' y <- 2 * x_signal + rnorm(200)
//'
//' screening_score(x_signal, y)
//' screening_score(x_null, y)
//'
//' @export
// [[Rcpp::export]]
double screening_score(SEXP x,
                       SEXP y,
                       const double tol = 1.4901161193847656e-08)
{
  const R_xlen_t n = XLENGTH(x);

  const double* x_data = REAL(x);
  const double* y_data = REAL(y);

  // Pass 1: means.
  double x_sum = 0.0;
  double y_sum = 0.0;

  for (R_xlen_t i = 0; i < n; ++i) {
    x_sum += x_data[i];
    y_sum += y_data[i];
  }

  const double inv_n = 1.0 / static_cast<double>(n);
  const double x_mean = x_sum * inv_n;
  const double y_mean = y_sum * inv_n;

  // Pass 2: all centered moments needed for D^2.
  double sxx = 0.0; // sum z^2
  double b = 0.0;   // sum z * u
  double m22 = 0.0; // sum z^2 * u^2
  double m31 = 0.0; // sum z^3 * u
  double m40 = 0.0; // sum z^4

  for (R_xlen_t i = 0; i < n; ++i) {
    const double z = x_data[i] - x_mean;
    const double u = y_data[i] - y_mean;
    const double z2 = z * z;
    const double u2 = u * u;

    sxx += z2;
    b += z * u;
    m22 += z2 * u2;
    m31 += z2 * z * u;
    m40 += z2 * z2;
  }

  const double tol2 = tol * tol;
  if (sxx <= 0.0 || b * b <= tol2 * std::max(sxx, m22)) {
    return R_PosInf;
  }

  const double inv_b = 1.0 / b;
  const double inv_sxx = 1.0 / sxx;

  double d2 =
    m22 * inv_b * inv_b
  - 2.0 * m31 * inv_b * inv_sxx
  + m40 * inv_sxx * inv_sxx;

  if (d2 < 0.0) {
    d2 = 0.0;
  }

  return d2;
}
