set.seed(20260709)

# FPR-control stress test.
#
# Under the global null, predictors and response have zero mean association,
# but all observations share a random row-level volatility factor. Ordinary
# Pearson p-values treat the observations as homoskedastic and are expected to
# be anti-conservative. Spearman is less sensitive but still affected. The
# robust reciprocal statistic should remain closer to the nominal q level.

q_levels <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30)
n <- 2000L
p <- 1000L
n_rep <- 100L
scale_sd <- 0.5

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  getwd()
}
script_dir <- dirname(normalizePath(script_path))
package_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(script_dir, "fpr_control_common_volatility_n2000_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lib_dir <- file.path(tempdir(), "mfscreen_common_volatility_sim")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

install_status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "-l", shQuote(lib_dir), shQuote(package_dir)),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(install_status, "status")) &&
    attr(install_status, "status") != 0L) {
  stop(paste(install_status, collapse = "\n"))
}

.libPaths(c(lib_dir, .libPaths()))
library(mfscreen)

pearson_pvalues <- function(X, y) {
  n <- nrow(X)
  yc <- y - mean(y)
  xc <- sweep(X, 2L, colMeans(X), check.margin = FALSE)

  denom <- sqrt(colSums(xc^2) * sum(yc^2))
  r <- colSums(xc * yc) / denom
  r[!is.finite(r)] <- 0
  r <- pmin(pmax(r, -1 + 1e-15), 1 - 1e-15)

  t_stat <- abs(r) * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  2 * pt(-t_stat, df = n - 2)
}

spearman_pvalues <- function(X, y) {
  X_rank <- apply(X, 2L, rank, ties.method = "average")
  y_rank <- rank(y, ties.method = "average")
  pearson_pvalues(X_rank, y_rank)
}

selected_from_scores <- function(scores, method, q) {
  if (method == "mfscreen") {
    gamma <- qnorm(1 - q / 2)
    return(is.finite(scores) & scores <= 1 / gamma^2)
  }

  scores <= q
}

raw_results <- data.frame()
methods <- c("mfscreen", "Pearson SIS", "Spearman SIS")

for (replication in seq_len(n_rep)) {
  if (replication %% 10L == 0L || replication == 1L) {
    message("Common-volatility FPR replication ", replication, "/", n_rep)
  }

  volatility <- exp(scale_sd * rnorm(n))
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) * volatility
  y <- volatility * rnorm(n)

  method_scores <- list(
    mfscreen = screening_test_matrix(X, y, q = 0.10)$scores,
    `Pearson SIS` = pearson_pvalues(X, y),
    `Spearman SIS` = spearman_pvalues(X, y)
  )

  for (method in methods) {
    scores <- method_scores[[method]]
    for (q in q_levels) {
      selected <- selected_from_scores(scores, method, q)
      raw_results <- rbind(
        raw_results,
        data.frame(
          method = method,
          replication = replication,
          q = q,
          fpr = mean(selected)
        )
      )
    }
  }

  if (replication %% 10L == 0L) {
    write.csv(raw_results, file.path(out_dir, "raw_fpr_control.csv"),
              row.names = FALSE)
  }

  rm(X, y, method_scores)
  gc()
}

summary_results <- aggregate(
  fpr ~ method + q,
  data = raw_results,
  FUN = mean
)

summary_sd <- aggregate(
  fpr ~ method + q,
  data = raw_results,
  FUN = sd
)
names(summary_sd)[names(summary_sd) == "fpr"] <- "fpr_sd"

summary_results <- merge(summary_results, summary_sd, by = c("method", "q"))
summary_results$inflation <- summary_results$fpr / summary_results$q

write.csv(raw_results, file.path(out_dir, "raw_fpr_control.csv"),
          row.names = FALSE)
write.csv(summary_results, file.path(out_dir, "summary_fpr_control.csv"),
          row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  p_fpr <- ggplot(
    summary_results,
    aes(x = q, y = fpr, colour = method, shape = method)
  ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = "grey45") +
    geom_line() +
    geom_point(size = 2) +
    scale_x_log10(breaks = q_levels) +
    scale_y_log10() +
    labs(
      x = "Nominal q level",
      y = "Empirical false positive rate",
      colour = "Method",
      shape = "Method",
      title = "FPR control under common volatility",
      subtitle = paste0("Global null, n = ", n, ", p = ", p,
                        ", replications = ", n_rep,
                        ", log-volatility sd = ", scale_sd)
    ) +
    theme_bw(base_size = 11)

  p_inflation <- ggplot(
    summary_results,
    aes(x = q, y = inflation, colour = method, shape = method)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_line() +
    geom_point(size = 2) +
    scale_x_log10(breaks = q_levels) +
    labs(
      x = "Nominal q level",
      y = "FPR / q",
      colour = "Method",
      shape = "Method",
      title = "False-positive inflation under common volatility"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "fpr_control_by_q.png"), p_fpr,
         width = 8, height = 5, dpi = 150)
  ggsave(file.path(out_dir, "fpr_control_by_q.pdf"), p_fpr,
         width = 8, height = 5)
  ggsave(file.path(out_dir, "fpr_inflation_by_q.png"), p_inflation,
         width = 8, height = 5, dpi = 150)
}

print(summary_results[order(summary_results$q, summary_results$method), ])
