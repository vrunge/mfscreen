#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

//' Reciprocal model-free screening for a predictor matrix
//'
//' Computes the squared reciprocal model-free screening statistic
//' \eqn{D_j^2} for each column \eqn{j} of \code{X}, then selects predictors
//' whose score is below a reciprocal normal-theory threshold.
//'
//' For predictor \eqn{j}, define
//' \deqn{
//'   z_{ij} = X_{ij} - \overline{X}_j,
//'   \qquad
//'   u_i = y_i - \overline{y}.
//' }
//'
//' The returned score is
//' \deqn{
//'   D_j^2 =
//'   \frac{\sum_{i=1}^n z_{ij}^2 u_i^2}
//'        {\left(\sum_{i=1}^n z_{ij} u_i\right)^2}
//'   -
//'   \frac{
//'     2 \sum_{i=1}^n z_{ij}^3 u_i
//'   }{
//'     \left(\sum_{i=1}^n z_{ij} u_i\right)
//'     \left(\sum_{i=1}^n z_{ij}^2\right)
//'   }
//'   +
//'   \frac{\sum_{i=1}^n z_{ij}^4}
//'        {\left(\sum_{i=1}^n z_{ij}^2\right)^2}.
//' }
//'
//' The statistic is algebraically equivalent to the squared reciprocal of the
//' studentized marginal screening statistic:
//' \deqn{
//'   D_j^2 = \frac{1}{T_j^2}.
//' }
//'
//' Predictor \eqn{j} is selected when
//' \deqn{
//'   D_j^2 \leq \frac{1}{\gamma^2},
//'   \qquad
//'   \gamma = \Phi^{-1}(1 - q/2).
//' }
//'
//' Smaller values of \eqn{D_j^2} indicate stronger marginal association with
//' the response.
//'
//' @param X Numeric predictor matrix. Rows correspond to observations and
//'   columns correspond to candidate predictors.
//' @param y Numeric response vector with length equal to \code{nrow(X)}.
//' @param q Numeric scalar in \code{(0, 1)} giving the target two-sided
//'   false-positive rate. Defaults to \code{0.10}.
//' @param tol Non-negative numerical tolerance used to treat a marginal
//'   empirical covariance as effectively zero. Defaults to
//'   \code{sqrt(.Machine$double.eps)}.
//'
//' @return A list with components:
//' \describe{
//'   \item{scores}{Numeric vector of reciprocal screening scores
//'   \eqn{D_j^2}. Smaller values indicate stronger marginal association.}
//'   \item{threshold}{Selection threshold \eqn{1 / \gamma^2}.}
//'   \item{gamma}{Normal-theory threshold \eqn{\Phi^{-1}(1-q/2)}.}
//'   \item{q}{The supplied target false-positive rate.}
//'   \item{selected}{Logical vector indicating selected predictors.}
//'   \item{selected_indices}{One-based integer indices of selected
//'   predictors.}
//' }
//'
//' @details
//' The response is centered once. Each predictor column is processed
//' independently using centered empirical moments. Columns with zero empirical
//' variance, or with empirical covariance numerically indistinguishable from
//' zero, receive score \code{Inf} and are not selected.
//'
//' This low-level Rcpp implementation performs no input validation. The inputs
//' must satisfy: \code{X} and \code{y} contain finite values,
//' \code{length(y) == nrow(X)}, \code{q} is in \code{(0, 1)}, and
//' \code{tol >= 0}.
//'
//' @seealso
//' \code{\link{screening_score}} for the one-predictor reciprocal statistic.
//'
//' @examples
//' set.seed(123)
//'
//' n <- 300
//' p <- 500
//'
//' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
//' y <- 2 * X[, 1] - 1.5 * X[, 2] + rnorm(n)
//'
//' result <- screening_test_matrix(X, y, q = 0.10)
//'
//' result$selected_indices
//' head(result$scores)
//' result$threshold
//'
//' @export
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

