set.seed(20260707)

# Competitors:
# 1. mfscreen: reciprocal model-free screening statistic from this package.
# 2. Pearson SIS: classical marginal Pearson correlation screening.
# 3. Spearman SIS: marginal rank-correlation screening.
#
# The simulation fixes n and p, varies the nominal q level, and estimates:
#   FPR = false positives / inactive variables
#   SSP = P(all active variables are selected)
#
# Two data-generating models are used:
#   Parametric: linear Gaussian model.
#   Nonparametric: nonlinear monotone/additive signal with heteroskedastic noise.

q_levels <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30)
n <- 300L
p <- 1000L
n_rep <- 150L
active <- seq_len(5L)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  getwd()
}
script_dir <- dirname(normalizePath(script_path))
package_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(script_dir, "fpr_ssp_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lib_dir <- file.path(tempdir(), "mfscreen_competitor_sim")
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

make_design <- function(n, p) {
  matrix(rnorm(n * p), nrow = n, ncol = p)
}

make_response <- function(X, model) {
  if (model == "Parametric linear") {
    beta <- c(1.2, -1.0, 0.9, -0.8, 0.7)
    return(drop(X[, active] %*% beta) + rnorm(nrow(X), sd = 1.5))
  }

  signal <- 1.6 * sin(X[, 1]) +
    1.2 * tanh(X[, 2]) +
    0.9 * X[, 3]^3 / (1 + X[, 3]^2) -
    0.8 * exp(0.35 * X[, 4]) +
    0.7 * X[, 5]
  sigma <- 0.7 + 0.5 * abs(X[, 1]) + 0.25 * X[, 2]^2

  signal + sigma * rnorm(nrow(X), sd = 0.8)
}

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

mfscreen_scores <- function(X, y) {
  screening_test_matrix(X, y, q = 0.10)$scores
}

selected_from_scores <- function(scores, method, q) {
  if (method == "mfscreen") {
    gamma <- qnorm(1 - q / 2)
    return(is.finite(scores) & scores <= 1 / gamma^2)
  }

  scores <= q
}

evaluate_selection <- function(selected, active, p) {
  inactive <- setdiff(seq_len(p), active)
  c(
    fpr = sum(selected[inactive]) / length(inactive),
    ssp = as.numeric(all(selected[active]))
  )
}

raw_results <- data.frame()

models <- c("Parametric linear", "Nonparametric nonlinear")
methods <- c("mfscreen", "Pearson SIS", "Spearman SIS")

for (model in models) {
  for (replication in seq_len(n_rep)) {
    message("Model: ", model, " | replication ", replication, "/", n_rep)

    X <- make_design(n, p)
    y <- make_response(X, model)

    method_scores <- list(
      mfscreen = mfscreen_scores(X, y),
      `Pearson SIS` = pearson_pvalues(X, y),
      `Spearman SIS` = spearman_pvalues(X, y)
    )

    for (method in methods) {
      scores <- method_scores[[method]]
      for (q in q_levels) {
        selected <- selected_from_scores(scores, method, q)
        metrics <- evaluate_selection(selected, active, p)

        raw_results <- rbind(
          raw_results,
          data.frame(
            model = model,
            method = method,
            replication = replication,
            q = q,
            fpr = metrics[["fpr"]],
            ssp = metrics[["ssp"]]
          )
        )
      }
    }

    if (replication %% 10L == 0L) {
      write.csv(raw_results, file.path(out_dir, "raw_fpr_ssp.csv"),
                row.names = FALSE)
    }

    rm(X, y, method_scores)
    gc()
  }
}

summary_results <- aggregate(
  cbind(fpr, ssp) ~ model + method + q,
  data = raw_results,
  FUN = mean
)

summary_sd <- aggregate(
  cbind(fpr, ssp) ~ model + method + q,
  data = raw_results,
  FUN = sd
)
names(summary_sd)[names(summary_sd) == "fpr"] <- "fpr_sd"
names(summary_sd)[names(summary_sd) == "ssp"] <- "ssp_sd"

summary_results <- merge(
  summary_results,
  summary_sd,
  by = c("model", "method", "q")
)

write.csv(raw_results, file.path(out_dir, "raw_fpr_ssp.csv"),
          row.names = FALSE)
write.csv(summary_results, file.path(out_dir, "summary_fpr_ssp.csv"),
          row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  summary_results$q_label <- ifelse(
    summary_results$q %in% c(0.001, 0.01, 0.05, 0.10, 0.30),
    as.character(summary_results$q),
    ""
  )

  p_curve <- ggplot(
    summary_results,
    aes(x = fpr, y = ssp, colour = method, shape = method)
  ) +
    geom_path(aes(group = method), linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_text(
      aes(label = q_label),
      size = 2.5,
      nudge_y = 0.025,
      show.legend = FALSE
    ) +
    facet_wrap(~ model) +
    coord_cartesian(xlim = c(0, max(summary_results$fpr) * 1.08),
                    ylim = c(0, 1.04)) +
    labs(
      x = "False positive rate",
      y = "Sure screening property",
      colour = "Method",
      shape = "Method",
      title = "FPR versus SSP across nominal q levels",
      subtitle = paste0("n = ", n, ", p = ", p,
                        ", active variables = ", length(active),
                        ", replications = ", n_rep)
    ) +
    theme_bw(base_size = 11)

  p_fpr_q <- ggplot(
    summary_results,
    aes(x = q, y = fpr, colour = method, shape = method)
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ model) +
    scale_x_log10(breaks = q_levels) +
    labs(
      x = "Nominal q level",
      y = "False positive rate",
      colour = "Method",
      shape = "Method",
      title = "False positive rate by q"
    ) +
    theme_bw(base_size = 11)

  p_ssp_q <- ggplot(
    summary_results,
    aes(x = q, y = ssp, colour = method, shape = method)
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ model) +
    scale_x_log10(breaks = q_levels) +
    coord_cartesian(ylim = c(0, 1.04)) +
    labs(
      x = "Nominal q level",
      y = "Sure screening property",
      colour = "Method",
      shape = "Method",
      title = "Sure screening property by q"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "fpr_vs_ssp.png"), p_curve,
         width = 10, height = 5.5, dpi = 150)
  ggsave(file.path(out_dir, "fpr_vs_ssp.pdf"), p_curve,
         width = 10, height = 5.5)
  ggsave(file.path(out_dir, "fpr_by_q.png"), p_fpr_q,
         width = 10, height = 5.5, dpi = 150)
  ggsave(file.path(out_dir, "ssp_by_q.png"), p_ssp_q,
         width = 10, height = 5.5, dpi = 150)
}

print(summary_results[order(summary_results$model, summary_results$method,
                            summary_results$q), ])
