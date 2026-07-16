set.seed(123)

q <- 0.10
n_rep <- 12L
min_batch_seconds <- 0.05
max_batch_calls <- 1000L

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

grid <- expand.grid(
  n = c(500L, 2000L, 5000L),
  p = c(200L, 1000L, 5000L),
  KEEP.OUT.ATTRS = FALSE
)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  getwd()
}
script_dir <- dirname(normalizePath(script_path))
package_dir <- normalizePath(file.path(script_dir, "..", ".."))
workspace_dir <- normalizePath(file.path(package_dir, ".."))
mfscreen0_dir <- file.path(workspace_dir, "mfscreen0")
mfscreen_dir <- file.path(workspace_dir, "mfscreen")

out_dir <- file.path(script_dir, "mfscreen0_vs_mfscreen_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lib_mfscreen0 <- file.path(tempdir(), "mfscreen0_benchmark")
lib_mfscreen <- file.path(tempdir(), "mfscreen_benchmark")
dir.create(lib_mfscreen0, recursive = TRUE, showWarnings = FALSE)
dir.create(lib_mfscreen, recursive = TRUE, showWarnings = FALSE)

install_pkg <- function(pkg_dir, lib_dir) {
  cmd <- file.path(R.home("bin"), "R")
  status <- system2(
    cmd,
    c("CMD", "INSTALL", "-l", shQuote(lib_dir), shQuote(pkg_dir)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop(paste(status, collapse = "\n"))
  }
}

unload_screening_packages <- function() {
  for (pkg in c("mfscreen0", "mfscreen")) {
    search_name <- paste0("package:", pkg)
    if (search_name %in% search()) {
      detach(search_name, unload = TRUE, character.only = TRUE)
    }
    if (pkg %in% loadedNamespaces()) {
      unloadNamespace(pkg)
    }
  }
}

time_implementation <- function(package, lib_dir, X, y, q, n_rep) {
  unload_screening_packages()
  library(package, lib.loc = lib_dir, character.only = TRUE)

  first_result <- getExportedValue(package, "screening_test_matrix")(X, y, q = q)
  gc()

  batch_calls <- 1L
  repeat {
    elapsed <- system.time({
      for (i in seq_len(batch_calls)) {
        invisible(getExportedValue(package, "screening_test_matrix")(X, y, q = q))
      }
    })[["elapsed"]]

    if (elapsed >= min_batch_seconds || batch_calls >= max_batch_calls) {
      break
    }
    batch_calls <- min(max_batch_calls, batch_calls * 2L)
  }

  times <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    elapsed <- system.time({
      for (i in seq_len(batch_calls)) {
        invisible(getExportedValue(package, "screening_test_matrix")(X, y, q = q))
      }
    })[["elapsed"]]

    times[r] <- elapsed / batch_calls
  }

  unload_screening_packages()
  list(
    times = times,
    selected = first_result$selected_indices,
    batch_calls = batch_calls
  )
}

install_pkg(mfscreen0_dir, lib_mfscreen0)
install_pkg(mfscreen_dir, lib_mfscreen)

raw_results <- data.frame()
selection_checks <- data.frame()

for (case in seq_len(nrow(grid))) {
  n <- grid$n[case]
  p <- grid$p[case]
  message("Running n = ", n, ", p = ", p)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- 2.0 * X[, 1L] - 1.5 * X[, min(2L, p)] + rnorm(n)

  base <- time_implementation("mfscreen0", lib_mfscreen0, X, y, q, n_rep)
  fast <- time_implementation("mfscreen", lib_mfscreen, X, y, q, n_rep)

  raw_results <- rbind(
    raw_results,
    data.frame(
      n = n,
      p = p,
      replication = seq_len(n_rep),
      implementation = "mfscreen0",
      time_seconds = base$times,
      batch_calls = base$batch_calls
    ),
    data.frame(
      n = n,
      p = p,
      replication = seq_len(n_rep),
      implementation = "mfscreen",
      time_seconds = fast$times,
      batch_calls = fast$batch_calls
    )
  )

  selection_checks <- rbind(
    selection_checks,
    data.frame(
      n = n,
      p = p,
      same_selected_indices = identical(base$selected, fast$selected)
    )
  )

  rm(X, y, base, fast)
  gc()
}

summary_results <- do.call(
  rbind,
  lapply(split(raw_results, list(raw_results$n, raw_results$p,
                                 raw_results$implementation), drop = TRUE),
         function(d) {
           data.frame(
             n = d$n[1L],
             p = d$p[1L],
             implementation = d$implementation[1L],
             mean_seconds = mean(d$time_seconds),
             sd_seconds = sd(d$time_seconds),
             min_seconds = min(d$time_seconds),
             max_seconds = max(d$time_seconds),
             batch_calls = d$batch_calls[1L],
             n_rep = nrow(d)
           )
         })
)
row.names(summary_results) <- NULL

summary_wide <- merge(
  summary_results[summary_results$implementation == "mfscreen0",
                  c("n", "p", "mean_seconds", "sd_seconds")],
  summary_results[summary_results$implementation == "mfscreen",
                  c("n", "p", "mean_seconds", "sd_seconds")],
  by = c("n", "p"),
  suffixes = c("_mfscreen0", "_mfscreen")
)
summary_wide$speedup <- summary_wide$mean_seconds_mfscreen0 /
  summary_wide$mean_seconds_mfscreen

write.csv(raw_results, file.path(out_dir, "raw_timings.csv"), row.names = FALSE)
write.csv(summary_results, file.path(out_dir, "summary_timings.csv"),
          row.names = FALSE)
write.csv(summary_wide, file.path(out_dir, "speedup_summary.csv"),
          row.names = FALSE)
write.csv(selection_checks, file.path(out_dir, "selection_checks.csv"),
          row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  summary_results$case <- paste0("n=", summary_results$n, ", p=",
                                 summary_results$p)
  summary_results$case <- factor(
    summary_results$case,
    levels = paste0("n=", grid$n, ", p=", grid$p)
  )

  p_time <- ggplot(
    summary_results,
    aes(x = case, y = mean_seconds, colour = implementation,
        group = implementation)
  ) +
    geom_point(size = 2) +
    geom_line() +
    geom_errorbar(aes(ymin = pmax(0, mean_seconds - sd_seconds),
                      ymax = mean_seconds + sd_seconds),
                  width = 0.2) +
    scale_y_log10() +
    labs(
      x = "(n, p)",
      y = "Mean elapsed time, seconds (log scale)",
      colour = "Implementation",
      title = "mfscreen0 vs mfscreen elapsed time",
      subtitle = paste("Error bars show +/- 1 standard deviation over",
                       n_rep, "replications")
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_speedup <- ggplot(
    summary_wide,
    aes(x = factor(p), y = factor(n), fill = speedup)
  ) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1fx", speedup)), colour = "white") +
    scale_fill_viridis_c(option = "C", trans = "log10") +
    labs(
      x = "p",
      y = "n",
      fill = "Speedup",
      title = "mfscreen speedup over mfscreen0"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "elapsed_time_by_np.png"), p_time,
         width = 10, height = 6, dpi = 150)
  ggsave(file.path(out_dir, "elapsed_time_by_np.pdf"), p_time,
         width = 10, height = 6)
  ggsave(file.path(out_dir, "speedup_heatmap.png"), p_speedup,
         width = 7, height = 5, dpi = 150)
  ggsave(file.path(out_dir, "speedup_heatmap.pdf"), p_speedup,
         width = 7, height = 5)
} else {
  png(file.path(out_dir, "elapsed_time_by_np.png"), width = 1200, height = 800)
  op <- par(mar = c(8, 5, 4, 2))
  cases <- paste0("n=", summary_results$n, ", p=", summary_results$p)
  plot(seq_along(cases), summary_results$mean_seconds, log = "y",
       xaxt = "n", xlab = "(n, p)", ylab = "Mean elapsed time (seconds)",
       pch = ifelse(summary_results$implementation == "mfscreen0", 16, 17),
       col = ifelse(summary_results$implementation == "mfscreen0", 1, 2))
  axis(1, seq_along(cases), cases, las = 2, cex.axis = 0.7)
  legend("topleft", legend = c("mfscreen0", "mfscreen"),
         pch = c(16, 17), col = c(1, 2), bty = "n")
  par(op)
  dev.off()
}

print(summary_results[order(summary_results$n, summary_results$p,
                            summary_results$implementation), ])
print(summary_wide[order(summary_wide$n, summary_wide$p), ])
print(selection_checks)
