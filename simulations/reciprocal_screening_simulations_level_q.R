# ============================================================================
# Monte Carlo study: level-q marginal screening rules
#
# Methods and level-q decision rules
#   proposed: select D_j^2 <= 1 / qnorm(1 - q/2)^2.
#             Since D_j^2 = 1 / T_HC0,j^2, this is |T_HC0,j| >= z_{1-q/2}.
#   pearson : select |r_j| >= r_q, where
#             r_q = sqrt(t_q^2 / (t_q^2 + n - 2)), t_q = qt(1-q/2,n-2).
#             This is exactly the usual two-sided Pearson t-test at level q.
#   dcor    : select if the distance-covariance independence test has p <= q.
#             Default dcorT.test is fast asymptotic calibration. Set
#             dcor_calibration = "permutation" for finite-sample permutation
#             calibration (recommended when matching null FPR is critical).
#
# Metrics: SSP, TPR, mean FP, FPR, FWER, mean selected size, FDR, precision,
#          active-variable MIP, average active-variable rank, and runtime.
# ============================================================================

# install.packages("energy")
suppressPackageStartupMessages(library(energy))

make_design <- function(n, p, design = c("independent", "ar1"), rho = 0.3,
                        heavy_tails = FALSE, df_t = 5) {
  design <- match.arg(design)
  if (design == "independent") {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  } else {
    X <- matrix(0, n, p)
    X[, 1] <- rnorm(n)
    innovation_sd <- sqrt(1 - rho^2)
    for (j in 2:p) X[, j] <- rho * X[, j - 1] + innovation_sd * rnorm(n)
  }
  if (heavy_tails) X <- scale(qt(pnorm(X), df = df_t))
  colnames(X) <- paste0("X", seq_len(p))
  X
}

make_response <- function(X, scenario = c("N", "L", "C", "Q", "M"),
                          theta = c(1, 1, 1), sigma = 1) {
  scenario <- match.arg(scenario)
  eps <- rnorm(nrow(X), sd = sigma)
  switch(scenario,
    N = list(y = eps, active = integer(0), label = "Global null"),
    L = list(y = theta[1] * X[, 1] + theta[2] * X[, 2] + eps,
             active = c(1L, 2L), label = "Linear"),
    C = list(y = theta[1] * X[, 1] + theta[2] * X[, 2]^3 + eps,
             active = c(1L, 2L), label = "Cubic"),
    Q = list(y = theta[1] * X[, 1] + theta[2] * (X[, 2]^2 - 1) + eps,
             active = c(1L, 2L), label = "Quadratic"),
    M = list(y = theta[1] * X[, 1] + theta[2] * X[, 2]^3 +
               theta[3] * (X[, 3]^2 - 1) + eps,
             active = c(1L, 2L, 3L), label = "Mixed")
  )
}

proposed_scores <- function(X, y) {
  p <- ncol(X)
  xc <- sweep(X, 2, colMeans(X), "-")
  yc <- y - mean(y)
  sxx <- colSums(xc^2)
  sxy <- colSums(xc * yc)
  valid <- sxx > .Machine$double.eps
  beta <- numeric(p); beta[valid] <- sxy[valid] / sxx[valid]
  residuals <- sweep(sweep(xc, 2, beta, "*"), 1, yc, function(f, yy) yy - f)
  sandwich_num <- colSums(xc^2 * residuals^2)
  D2 <- rep(Inf, p); T <- rep(0, p)
  usable <- valid & abs(sxy) > .Machine$double.eps
  finite <- usable & sandwich_num > .Machine$double.eps
  D2[finite] <- sandwich_num[finite] / sxy[finite]^2
  T[finite] <- sxy[finite] / sqrt(sandwich_num[finite])
  exact <- usable & sandwich_num <= .Machine$double.eps
  D2[exact] <- 0; T[exact] <- sign(sxy[exact]) * Inf
  list(D2 = D2, T = T)
}

pearson_scores <- function(X, y) {
  n <- nrow(X)
  xc <- sweep(X, 2, colMeans(X), "-")
  yc <- y - mean(y)
  den <- sqrt(colSums(xc^2) * sum(yc^2))
  r <- colSums(xc * yc) / den
  r[!is.finite(r)] <- 0
  list(r = pmin(1, pmax(-1, r)))
}

dcor_scores <- function(X, y, calibration = c("t", "permutation"), B = 499) {
  calibration <- match.arg(calibration)
  p <- ncol(X); stat <- numeric(p); pval <- numeric(p)
  for (j in seq_len(p)) {
    stat[j] <- tryCatch(energy::dcor(X[, j], y), error = function(e) 0)
    if (calibration == "t") {
      pval[j] <- tryCatch(energy::dcorT.test(X[, j], y)$p.value,
                          error = function(e) 1)
    } else {
      # energy::dcov.test is a permutation/randomization test under independence.
      pval[j] <- tryCatch(energy::dcov.test(X[, j], y, R = B)$p.value,
                          error = function(e) 1)
    }
  }
  stat[!is.finite(stat)] <- 0; pval[!is.finite(pval)] <- 1
  list(dcor = stat, p_value = pval)
}

screen_once <- function(X, y, q = 0.05,
                        dcor_calibration = c("t", "permutation"), B = 499,
                        d = ceiling(nrow(X) / log(nrow(X)))) {
  dcor_calibration <- match.arg(dcor_calibration)
  n <- nrow(X); p <- ncol(X); d <- min(d, p)
  zq <- qnorm(1 - q / 2)
  reciprocal_cutoff <- 1 / zq^2
  tq <- qt(1 - q / 2, df = n - 2)
  pearson_cutoff <- sqrt(tq^2 / (tq^2 + n - 2))

  t0 <- proc.time()[["elapsed"]]; prop <- proposed_scores(X, y)
  prop_time <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]; pear <- pearson_scores(X, y)
  pear_time <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]; dc <- dcor_scores(X, y, dcor_calibration, B)
  dc_time <- proc.time()[["elapsed"]] - t0

  list(
    cutoffs = c(reciprocal_D2 = reciprocal_cutoff, pearson_abs_r = pearson_cutoff,
                dcor_p_value = q),
    ranks = list(proposed = rank(prop$D2, ties.method = "average"),
                 pearson = rank(-abs(pear$r), ties.method = "average"),
                 dcor = rank(-dc$dcor, ties.method = "average")),
    selected = list(
      level_q = list(proposed = which(prop$D2 <= reciprocal_cutoff),
                     pearson = which(abs(pear$r) >= pearson_cutoff),
                     dcor = which(dc$p_value <= q)),
      top_d = list(proposed = order(prop$D2)[seq_len(d)],
                   pearson = order(abs(pear$r), decreasing = TRUE)[seq_len(d)],
                   dcor = order(dc$dcor, decreasing = TRUE)[seq_len(d)])
    ),
    timing = c(proposed = prop_time, pearson = pear_time, dcor = dc_time)
  )
}

metrics_one <- function(selected, ranks, active, p) {
  inactive <- setdiff(seq_len(p), active)
  tp <- sum(selected %in% active); fp <- sum(selected %in% inactive); m <- length(selected)
  list(SSP = if (length(active)) as.numeric(all(active %in% selected)) else NA_real_,
       TPR = if (length(active)) tp / length(active) else NA_real_,
       FP = fp, FPR = fp / length(inactive), FWER = as.numeric(fp > 0),
       selected_size = m, FDP = if (m) fp / m else 0,
       precision = if (m) tp / m else if (length(active)) 0 else 1,
       MIP = if (length(active)) as.numeric(active %in% selected) else numeric(),
       active_rank = if (length(active)) ranks[active] else numeric())
}

run_one_model <- function(n = 200, p = 500, scenario = "M", design = "independent",
                          rho = 0.3, theta = c(1, 1, 1), sigma = 1,
                          R = 500, q = 0.05, dcor_calibration = "t", B = 499,
                          heavy_tails = FALSE, df_t = 5, seed = 2026, progress = TRUE) {
  methods <- c("proposed", "pearson", "dcor"); rules <- c("level_q", "top_d")
  active <- switch(scenario, N=integer(0), L=c(1L,2L), C=c(1L,2L),
                   Q=c(1L,2L), M=c(1L,2L,3L), stop("Unknown scenario"))
  set.seed(seed); d <- min(ceiling(n/log(n)), p)
  init_rule <- function() list(
    SSP=matrix(NA_real_,R,3,dimnames=list(NULL,methods)), TPR=matrix(NA_real_,R,3,dimnames=list(NULL,methods)),
    FP=matrix(NA_real_,R,3,dimnames=list(NULL,methods)), FPR=matrix(NA_real_,R,3,dimnames=list(NULL,methods)),
    FWER=matrix(NA_real_,R,3,dimnames=list(NULL,methods)), size=matrix(NA_real_,R,3,dimnames=list(NULL,methods)),
    FDP=matrix(NA_real_,R,3,dimnames=list(NULL,methods)), precision=matrix(NA_real_,R,3,dimnames=list(NULL,methods)),
    MIP=setNames(lapply(methods, function(.) matrix(NA_real_,R,length(active))),methods),
    rank=setNames(lapply(methods, function(.) matrix(NA_real_,R,length(active))),methods))
  raw <- setNames(lapply(rules, function(.) init_rule()), rules)
  runtime <- matrix(NA_real_,R,3,dimnames=list(NULL,methods))

  for (rep in seq_len(R)) {
    X <- make_design(n,p,design,rho,heavy_tails,df_t); dat <- make_response(X,scenario,theta,sigma)
    out <- screen_once(X,dat$y,q,dcor_calibration,B,d)
    for (rule in rules) for (method in methods) {
      m <- metrics_one(out$selected[[rule]][[method]], out$ranks[[method]], dat$active, p)
      raw[[rule]]$SSP[rep,method] <- m$SSP; raw[[rule]]$TPR[rep,method] <- m$TPR
      raw[[rule]]$FP[rep,method] <- m$FP; raw[[rule]]$FPR[rep,method] <- m$FPR
      raw[[rule]]$FWER[rep,method] <- m$FWER; raw[[rule]]$size[rep,method] <- m$selected_size
      raw[[rule]]$FDP[rep,method] <- m$FDP; raw[[rule]]$precision[rep,method] <- m$precision
      if (length(active)) { raw[[rule]]$MIP[[method]][rep,] <- m$MIP; raw[[rule]]$rank[[method]][rep,] <- m$active_rank }
    }
    runtime[rep,] <- out$timing[methods]
    if (progress && (rep %% max(1L,floor(R/10L)) == 0L || rep == R)) message(sprintf("%s: %d / %d",scenario,rep,R))
  }
  summary <- do.call(rbind,lapply(rules,function(rule) do.call(rbind,lapply(methods,function(method)
    data.frame(scenario=scenario,design=design,rule=rule,q=if(rule=="level_q") q else NA_real_,n=n,p=p,R=R,d=d,
      dcor_calibration=if(rule=="level_q") dcor_calibration else NA_character_,method=method,
      SSP=mean(raw[[rule]]$SSP[,method],na.rm=TRUE),TPR=mean(raw[[rule]]$TPR[,method],na.rm=TRUE),
      mean_FP=mean(raw[[rule]]$FP[,method]),FPR=mean(raw[[rule]]$FPR[,method]),FWER=mean(raw[[rule]]$FWER[,method]),
      mean_selected=mean(raw[[rule]]$size[,method]),FDR=mean(raw[[rule]]$FDP[,method]),
      precision=mean(raw[[rule]]$precision[,method]),mean_runtime_seconds=mean(runtime[,method]))))))
  active_table <- if (!length(active)) data.frame() else do.call(rbind,lapply(rules,function(rule) do.call(rbind,lapply(methods,function(method)
    data.frame(scenario=scenario,design=design,rule=rule,q=if(rule=="level_q")q else NA_real_,method=method,
      active_variable=paste0("X",active),MIP=colMeans(raw[[rule]]$MIP[[method]]),average_rank=colMeans(raw[[rule]]$rank[[method]]))))))
  list(summary=summary,active_results=active_table,raw=raw,runtime=runtime)
}

run_all_models <- function(n=200,p=500,scenarios=c("N","L","C","Q","M"),design="independent",rho=0.3,
                           theta=c(1,1,1),sigma=1,R=500,q=0.05,dcor_calibration="t",B=499,
                           heavy_tails=FALSE,df_t=5,seed=2026,progress=TRUE) {
  ans <- lapply(seq_along(scenarios),function(k) run_one_model(n,p,scenarios[k],design,rho,theta,sigma,R,q,dcor_calibration,B,heavy_tails,df_t,seed+1000L*k,progress))
  names(ans) <- scenarios
  s <- do.call(rbind,lapply(ans,`[[`,"summary")); a <- do.call(rbind,lapply(ans,`[[`,"active_results"))
  list(threshold_report=subset(s,rule=="level_q"),ranking_report=subset(s,rule=="top_d"),active_variable_report=a,all_results=ans)
}

if (sys.nframe()==0L) {
  out <- run_all_models(n=200,p=500,R=100,q=0.05,dcor_calibration="t",B=499)
  cat("\n=== Level-q threshold report ===\n"); print(out$threshold_report,row.names=FALSE)
  cat("\n=== Active-variable report ===\n"); print(out$active_variable_report,row.names=FALSE)
  cat("\n=== Top-d ranking report ===\n"); print(out$ranking_report,row.names=FALSE)
  write.csv(out$threshold_report,"level_q_threshold_report.csv",row.names=FALSE)
  write.csv(out$active_variable_report,"active_variable_report.csv",row.names=FALSE)
  write.csv(out$ranking_report,"top_d_report.csv",row.names=FALSE)
}
