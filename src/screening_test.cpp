#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Preconditions:
// - X and y are finite doubles.
// - length(y) == nrow(X).
// - q in (0, 1), tol >= 0.
// - No input validation.
//
// [[Rcpp::export]]
List screening_test_matrix(
    const NumericMatrix& X,
    const NumericVector& y,
    const double q = 0.10,
    const double tol = 1.4901161193847656e-08
)
{
  const int n = X.nrow();
  const int p = X.ncol();

  const double* x_data = REAL(X);
  const double* y_data = REAL(y);

  // Center y once.
  double y_mean = 0.0;
  for (int i = 0; i < n; ++i)
  {
    y_mean += y_data[i];
  }
  y_mean /= static_cast<double>(n);

  NumericVector y_centered(n);
  double* yc_data = REAL(y_centered);

  ///////  ///////
  ///////  ///////
  for (int i = 0; i < n; ++i)
  {
    yc_data[i] = y_data[i] - y_mean;
  }
  ///////  ///////
  ///////  ///////


  const double gamma = R::qnorm5(1.0 - q / 2.0, 0.0, 1.0, 1, 0);
  const double threshold = 1.0 / (gamma * gamma);

  NumericVector scores(p);
  LogicalVector selected(p);

  std::vector<int> selected_buffer;
  selected_buffer.reserve(static_cast<std::size_t>(p));

  const double inv_n = 1.0 / static_cast<double>(n);

  /////  /////  /////
  /////  /////  /////
  /////  /////  /////

  for (int j = 0; j < p; ++j)
  {
    const double* xj = x_data + static_cast<R_xlen_t>(j) * n;

    ///////  ///////
    ///////  ///////
    // Pass 1: predictor mean.
    double x_sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
      x_sum += xj[i];
    }

    const double x_mean = x_sum * inv_n;

    // Pass 2: all centered moments needed for D^2.
    double sxx = 0.0;   // sum z^2
    double b = 0.0;     // sum z * u
    double m22 = 0.0;   // sum z^2 * u^2
    double m31 = 0.0;   // sum z^3 * u
    double m40 = 0.0;   // sum z^4

    double z;
    double u;

    double z2;
    double u2;

    for (int i = 0; i < n; ++i)
    {
      z = xj[i] - x_mean;
      u = yc_data[i];
      z2 = z * z;
      u2 = u * u;

      sxx += z2;
      b += z * u;
      m22 += z2 * u2;
      m31 += z2 * z * u;
      m40 += z2 * z2;
    }

    if (sxx <= 0.0)
    {
      scores[j] = R_PosInf;
      selected[j] = false;
      continue;
    }

    // Equivalent to the original A-based tolerance check.
    const double b_limit = tol * std::max(std::sqrt(sxx), std::sqrt(m22));

    if (std::fabs(b) <= b_limit)
    {
      scores[j] = R_PosInf;
      selected[j] = false;
      continue;
    }

    const double inv_b = 1.0 / b;
    const double inv_sxx = 1.0 / sxx;

    double d2 =
      m22 * inv_b * inv_b
    - 2.0 * m31 * inv_b * inv_sxx
    + m40 * inv_sxx * inv_sxx;

    // D^2 is theoretically non-negative; protect against tiny round-off.
    if (d2 < 0.0) {d2 = 0.0;}
    scores[j] = d2;

    const bool keep = d2 <= threshold;
    selected[j] = keep;

    if (keep) {selected_buffer.push_back(j + 1);}
  }

  /////  /////  /////
  /////  /////  /////
  /////  /////  /////

  IntegerVector selected_indices(
      selected_buffer.begin(),
      selected_buffer.end()
  );

  return List::create(
    _["scores"] = scores,
    _["threshold"] = threshold,
    _["gamma"] = gamma,
    _["q"] = q,
    _["selected"] = selected,
    _["selected_indices"] = selected_indices
  );
}
