#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "screening_score.h"

using namespace Rcpp;


double screening_score_core(SEXP x,
                            SEXP y,
                            const double tol)
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

  // Pass 2:
  // sxx       = sum (x_i - x_mean)^2
  // b         = sum (x_i - x_mean)(y_i - y_mean)
  // sum_z2_y2 = sum (x_i - x_mean)^2 (y_i - y_mean)^2
  double sxx = 0.0;
  double b = 0.0;
  double sum_z2_y2 = 0.0;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double z = x_data[i] - x_mean;
    const double yc = y_data[i] - y_mean;
    const double z2 = z * z;

    sxx += z2;
    b += z * yc;
    sum_z2_y2 += z2 * yc * yc;
  }

  // Equivalent to:
  // abs(A) <= tol * max(1, sqrt(sum(x_tilde^2 * y_tilde^2))).
  const double b_limit =
    tol * std::max(std::sqrt(sxx), std::sqrt(sum_z2_y2));

  if (std::fabs(b) <= b_limit) {
    return R_PosInf;
  }

  const double inv_b = 1.0 / b;
  const double inv_sxx = 1.0 / sxx;

  // Pass 3:
  // D^2 = sum [z_i * y_centered_i / b - z_i^2 / sxx]^2.
  double d2 = 0.0;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double z = x_data[i] - x_mean;
    const double yc = y_data[i] - y_mean;

    const double difference =
      z * (yc * inv_b - z * inv_sxx);

    d2 += difference * difference;
  }

  return d2;
}
