set.seed(20260710)

# Original efficiency study:
# Can a screening method retain true active predictors while avoiding
# false positives created by a shared volatility factor?
#
# Two scenarios are compared:
# 1. Clean homoskedastic noise.
# 2. Shared volatility: inactive predictors and the noise component of Y
#    share the same row-level volatility A_i. They have zero mean association,
#    but ordinary correlation p-values become anti-conservative.

q_levels <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20, 0.30)
n <- 1000L
p <- 1000L
n_rep <- 150L
active <- seq_len(5L)
beta <- 0.30 * c(1.0, -0.9, 0.8, -0.7, 0.6)
volatility_sd <- 0.5

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  getwd()
}
script_dir <- dirname(normalizePath(script_path))
package_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(script_dir, "shared_volatility_active_signal_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lib_dir <- file.path(tempdir(), "mfscreen_shared_volatility_study")
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

make_data <- function(scenario) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X[, active] <- matrix(rnorm(n * length(active)),
                        nrow = n,
                        ncol = length(active))

  signal <- drop(X[, active] %*% beta)

  if (scenario == "Clean homoskedastic") {
    y <- signal + rnorm(n)
    return(list(X = X, y = y))
  }

  volatility <- exp(volatility_sd * rnorm(n))
  X[, -active] <- X[, -active] * volatility
  y <- signal + volatility * rnorm(n)

  list(X = X, y = y)
}

selected_from_scores <- function(scores, method, q) {
  if (method == "mfscreen") {
    gamma <- qnorm(1 - q / 2)
    return(is.finite(scores) & scores <= 1 / gamma^2)
  }

  scores <= q
}

evaluate_selection <- function(selected) {
  inactive <- setdiff(seq_len(p), active)
  c(
    fpr = mean(selected[inactive]),
    ssp = as.numeric(all(selected[active])),
    selected_count = sum(selected)
  )
}

raw_results <- data.frame()
scenarios <- c("Clean homoskedastic", "Shared volatility")
methods <- c("mfscreen", "Pearson SIS", "Spearman SIS")

for (scenario in scenarios) {
  for (replication in seq_len(n_rep)) {
    if (replication %% 10L == 0L || replication == 1L) {
      message("Scenario: ", scenario, " | replication ",
              replication, "/", n_rep)
    }

    data <- make_data(scenario)
    X <- data$X
    y <- data$y

    method_scores <- list(
      mfscreen = screening_test_matrix(X, y, q = 0.10)$scores,
      `Pearson SIS` = pearson_pvalues(X, y),
      `Spearman SIS` = spearman_pvalues(X, y)
    )

    for (method in methods) {
      scores <- method_scores[[method]]
      for (q in q_levels) {
        selected <- selected_from_scores(scores, method, q)
        metrics <- evaluate_selection(selected)

        raw_results <- rbind(
          raw_results,
          data.frame(
            scenario = scenario,
            method = method,
            replication = replication,
            q = q,
            fpr = metrics[["fpr"]],
            ssp = metrics[["ssp"]],
            selected_count = metrics[["selected_count"]]
          )
        )
      }
    }

    if (replication %% 10L == 0L) {
      write.csv(raw_results, file.path(out_dir, "raw_results.csv"),
                row.names = FALSE)
    }

    rm(data, X, y, method_scores)
    gc()
  }
}

summary_results <- aggregate(
  cbind(fpr, ssp, selected_count) ~ scenario + method + q,
  data = raw_results,
  FUN = mean
)

summary_sd <- aggregate(
  cbind(fpr, ssp, selected_count) ~ scenario + method + q,
  data = raw_results,
  FUN = sd
)
names(summary_sd)[names(summary_sd) == "fpr"] <- "fpr_sd"
names(summary_sd)[names(summary_sd) == "ssp"] <- "ssp_sd"
names(summary_sd)[names(summary_sd) == "selected_count"] <- "selected_count_sd"

summary_results <- merge(
  summary_results,
  summary_sd,
  by = c("scenario", "method", "q")
)
summary_results$fpr_inflation <- summary_results$fpr / summary_results$q

write.csv(raw_results, file.path(out_dir, "raw_results.csv"),
          row.names = FALSE)
write.csv(summary_results, file.path(out_dir, "summary_results.csv"),
          row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  summary_results$q_label <- ifelse(
    summary_results$q %in% c(0.001, 0.01, 0.05, 0.10, 0.30),
    as.character(summary_results$q),
    ""
  )

  p_frontier <- ggplot(
    summary_results,
    aes(x = fpr, y = ssp, colour = method, shape = method)
  ) +
    geom_path(aes(group = method), linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_text(aes(label = q_label), size = 2.5, nudge_y = 0.025,
              show.legend = FALSE) +
    facet_wrap(~ scenario) +
    coord_cartesian(xlim = c(0, max(summary_results$fpr) * 1.05),
                    ylim = c(0, 1.04)) +
    labs(
      x = "False positive rate",
      y = "Sure screening property",
      colour = "Method",
      shape = "Method",
      title = "Efficiency frontier under shared volatility",
      subtitle = paste0("n = ", n, ", p = ", p,
                        ", replications = ", n_rep)
    ) +
    theme_bw(base_size = 11)

  p_fpr <- ggplot(
    summary_results,
    aes(x = q, y = fpr, colour = method, shape = method)
  ) +
    geom_line(
      data = data.frame(q = q_levels, fpr = q_levels),
      aes(x = q, y = fpr),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour = "grey45"
    ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ scenario) +
    scale_x_log10(breaks = q_levels) +
    labs(
      x = "Nominal q level",
      y = "False positive rate",
      colour = "Method",
      shape = "Method",
      title = "FPR control by q"
    ) +
    theme_bw(base_size = 11)

  p_selected <- ggplot(
    summary_results,
    aes(x = q, y = selected_count, colour = method, shape = method)
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ scenario) +
    scale_x_log10(breaks = q_levels) +
    labs(
      x = "Nominal q level",
      y = "Average selected variables",
      colour = "Method",
      shape = "Method",
      title = "Model size induced by screening"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "efficiency_frontier.png"), p_frontier,
         width = 10, height = 5.5, dpi = 150)
  ggsave(file.path(out_dir, "efficiency_frontier.pdf"), p_frontier,
         width = 10, height = 5.5)
  ggsave(file.path(out_dir, "fpr_by_q.png"), p_fpr,
         width = 10, height = 5.5, dpi = 150)
  ggsave(file.path(out_dir, "selected_count_by_q.png"), p_selected,
         width = 10, height = 5.5, dpi = 150)
}

print(summary_results[order(summary_results$scenario,
                            summary_results$q,
                            summary_results$method), ])
