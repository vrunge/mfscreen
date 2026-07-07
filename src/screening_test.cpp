#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <thread>
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

  std::vector<double> y_centered(static_cast<std::size_t>(n));
  std::vector<double> y_centered2(static_cast<std::size_t>(n));
  double y_centered_sum = 0.0;
  double y_centered2_sum = 0.0;
  for (int i = 0; i < n; ++i)
  {
    const double u = y_data[i] - y_mean;
    y_centered[i] = u;
    y_centered2[i] = u * u;
    y_centered_sum += u;
    y_centered2_sum += u * u;
  }

  const double gamma = R::qnorm5(1.0 - q / 2.0, 0.0, 1.0, 1, 0);
  const double threshold = 1.0 / (gamma * gamma);
  const double tol2 = tol * tol;

  const double inv_n = 1.0 / static_cast<double>(n);
  std::vector<double> score_data(static_cast<std::size_t>(p));
  std::vector<int> selected_data(static_cast<std::size_t>(p));

  const auto process_columns = [&](const int begin, const int end) {
    for (int j = begin; j < end; ++j)
    {
      const double* xj = x_data + static_cast<R_xlen_t>(j) * n;

      double sx = 0.0;
      double sx2 = 0.0;
      double sx3 = 0.0;
      double sx4 = 0.0;
      double sxu = 0.0;
      double sx2u = 0.0;
      double sx3u = 0.0;
      double sxu2 = 0.0;
      double sx2u2 = 0.0;

      for (int i = 0; i < n; ++i)
      {
        const double x = xj[i];
        const double u = y_centered[i];
        const double u2 = y_centered2[i];
        const double x2 = x * x;
        const double x3 = x2 * x;

        sx += x;
        sx2 += x2;
        sx3 += x3;
        sx4 += x2 * x2;
        sxu += x * u;
        sx2u += x2 * u;
        sx3u += x3 * u;
        sxu2 += x * u2;
        sx2u2 += x2 * u2;
      }

      const double x_mean = sx * inv_n;
      const double x_mean2 = x_mean * x_mean;
      const double x_mean3 = x_mean2 * x_mean;
      const double x_mean4 = x_mean2 * x_mean2;

      double sxx = sx2 - sx * sx * inv_n;
      double b = sxu - x_mean * y_centered_sum;
      double m22 = sx2u2 - 2.0 * x_mean * sxu2 +
        x_mean2 * y_centered2_sum;
      double m31 = sx3u - 3.0 * x_mean * sx2u +
        3.0 * x_mean2 * sxu - x_mean3 * y_centered_sum;
      double m40 = sx4 - 4.0 * x_mean * sx3 +
        6.0 * x_mean2 * sx2 - 4.0 * x_mean3 * sx +
        static_cast<double>(n) * x_mean4;

      if (sxx <= std::sqrt(std::numeric_limits<double>::epsilon()) * sx2 ||
          m22 < 0.0 ||
          m40 < 0.0)
      {
        sxx = 0.0;
        b = 0.0;
        m22 = 0.0;
        m31 = 0.0;
        m40 = 0.0;

        for (int i = 0; i < n; ++i)
        {
          const double z = xj[i] - x_mean;
          const double u = y_centered[i];
          const double z2 = z * z;

          sxx += z2;
          b += z * u;
          m22 += z2 * y_centered2[i];
          m31 += z2 * z * u;
          m40 += z2 * z2;
        }
      }

      if (sxx <= 0.0)
      {
        score_data[j] = R_PosInf;
        selected_data[j] = false;
        continue;
      }

      if (b * b <= tol2 * std::max(sxx, m22))
      {
        score_data[j] = R_PosInf;
        selected_data[j] = false;
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
      score_data[j] = d2;
      selected_data[j] = d2 <= threshold;
    }
  };

  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  const int max_threads = hardware_threads == 0 ? 1 :
    static_cast<int>(hardware_threads);
  const int thread_count = std::min(p, max_threads);

  if (thread_count <= 1 || p < 64)
  {
    process_columns(0, p);
  }
  else
  {
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(thread_count));

    const int chunk = (p + thread_count - 1) / thread_count;
    for (int t = 0; t < thread_count; ++t)
    {
      const int begin = t * chunk;
      const int end = std::min(p, begin + chunk);
      if (begin < end)
      {
        threads.emplace_back(process_columns, begin, end);
      }
    }

    for (std::thread& thread : threads)
    {
      thread.join();
    }
  }

  NumericVector scores(p);
  LogicalVector selected(p);
  std::vector<int> selected_buffer;
  selected_buffer.reserve(static_cast<std::size_t>(p));

  double* scores_out = REAL(scores);
  int* selected_out = LOGICAL(selected);
  for (int j = 0; j < p; ++j)
  {
    scores_out[j] = score_data[j];
    selected_out[j] = selected_data[j];
    if (selected_data[j])
    {
      selected_buffer.push_back(j + 1);
    }
  }

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
