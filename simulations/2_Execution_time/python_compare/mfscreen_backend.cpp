#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <thread>
#include <vector>

extern "C" void fill_data(const int n, const int p, double* x, double* y)
{
  for (int j = 0; j < p; ++j) {
    double* xj = x + static_cast<std::size_t>(j) * n;
    const double jj = static_cast<double>(j + 1);

    for (int i = 0; i < n; ++i) {
      const double ii = static_cast<double>(i + 1);
      xj[i] = std::sin(0.00013 * ii * jj) +
        std::cos(0.017 * ii + 0.011 * jj);
    }
  }

  const double* x1 = x;
  const double* x2 = p > 1 ? x + static_cast<std::size_t>(n) : x;
  for (int i = 0; i < n; ++i) {
    const double ii = static_cast<double>(i + 1);
    y[i] = 2.0 * x1[i] - 1.5 * x2[i] + 0.25 * std::sin(0.019 * ii);
  }
}

extern "C" int screening_test_matrix_cpp(
    const int n,
    const int p,
    const double* x_data,
    const double* y_data,
    const double q,
    const double tol,
    double* scores,
    int* selected)
{
  double y_mean = 0.0;
  for (int i = 0; i < n; ++i) {
    y_mean += y_data[i];
  }
  y_mean /= static_cast<double>(n);

  std::vector<double> y_centered(static_cast<std::size_t>(n));
  std::vector<double> y_centered2(static_cast<std::size_t>(n));
  double y_centered_sum = 0.0;
  double y_centered2_sum = 0.0;

  for (int i = 0; i < n; ++i) {
    const double u = y_data[i] - y_mean;
    y_centered[i] = u;
    y_centered2[i] = u * u;
    y_centered_sum += u;
    y_centered2_sum += u * u;
  }

  const double gamma = 1.959963984540054; // qnorm(1 - 0.10 / 2)
  const double threshold = q == 0.10 ? 1.0 / (gamma * gamma) : 0.0;
  const double tol2 = tol * tol;
  const double inv_n = 1.0 / static_cast<double>(n);

  const auto process_columns = [&](const int begin, const int end) {
    for (int j = begin; j < end; ++j) {
      const double* xj = x_data + static_cast<std::size_t>(j) * n;

      double sx = 0.0;
      double sx2 = 0.0;
      double sx3 = 0.0;
      double sx4 = 0.0;
      double sxu = 0.0;
      double sx2u = 0.0;
      double sx3u = 0.0;
      double sxu2 = 0.0;
      double sx2u2 = 0.0;

      for (int i = 0; i < n; ++i) {
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
          m40 < 0.0) {
        sxx = 0.0;
        b = 0.0;
        m22 = 0.0;
        m31 = 0.0;
        m40 = 0.0;

        for (int i = 0; i < n; ++i) {
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

      if (sxx <= 0.0 || b * b <= tol2 * std::max(sxx, m22)) {
        scores[j] = std::numeric_limits<double>::infinity();
        selected[j] = 0;
        continue;
      }

      const double inv_b = 1.0 / b;
      const double inv_sxx = 1.0 / sxx;
      double d2 = m22 * inv_b * inv_b -
        2.0 * m31 * inv_b * inv_sxx +
        m40 * inv_sxx * inv_sxx;

      if (d2 < 0.0) {
        d2 = 0.0;
      }

      scores[j] = d2;
      selected[j] = d2 <= threshold;
    }
  };

  const unsigned int hardware_threads = std::thread::hardware_concurrency();
  const int max_threads = hardware_threads == 0 ? 1 :
    static_cast<int>(hardware_threads);
  const int thread_count = std::min(p, max_threads);

  if (thread_count <= 1 || p < 64) {
    process_columns(0, p);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(thread_count));

    const int chunk = (p + thread_count - 1) / thread_count;
    for (int t = 0; t < thread_count; ++t) {
      const int begin = t * chunk;
      const int end = std::min(p, begin + chunk);
      if (begin < end) {
        threads.emplace_back(process_columns, begin, end);
      }
    }

    for (std::thread& thread : threads) {
      thread.join();
    }
  }

  int selected_count = 0;
  for (int j = 0; j < p; ++j) {
    selected_count += selected[j] != 0;
  }

  return selected_count;
}
