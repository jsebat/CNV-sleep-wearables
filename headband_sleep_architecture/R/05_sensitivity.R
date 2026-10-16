# =====================================================================
# 05_sensitivity.R -- robustness checks for the primary findings
#
# Three things could plausibly drive a spurious NREM-depth result, so each
# is checked explicitly and the results are written alongside the main table:
#
#   (a) DI is an unbounded ratio, N3 / (N1+N2+REM).  A night that is almost
#       all N3 has enormous leverage.  -> refit with increasingly strict caps.
#   (b) The choice of depth-index denominator.       -> compare all three.
#   (c) Movement-time and recording-length QC thresholds. -> vary them.
# =====================================================================

#' Refit the mirror and U-shape models for one outcome under a value cap.
sens_outlier_cap <- function(d, outcome, caps = c(Inf, 2, 1.5, 1),
                             verbose = TRUE) {
  rows <- list()
  for (cap in caps) {
    dd <- d[!is.na(d[[outcome]]) & d[[outcome]] <= cap, , drop = FALSE]
    for (locus in c("16p", "22q")) {
      s <- locus_subset(dd, locus)
      if (nrow(s) < 20 || length(unique(s$CopyNum)) < 3) next
      for (mod in c("mirror", "ushape")) {
        tm  <- if (mod == "mirror") "copy_c" else "copy_q"
        fit <- fit_lmm(build_formula(outcome, tm, s), s)
        if (is.null(fit)) next
        co <- summary(fit)$coefficients
        if (!tm %in% rownames(co)) next
        rows[[length(rows) + 1]] <- data.frame(
          outcome = outcome, check = "value cap",
          setting = if (is.infinite(cap)) "none (primary)" else paste0(outcome, " <= ", cap),
          locus = locus, model = mod, term = tm,
          n_nights = nrow(s), n_dropped = sum(!is.na(d[[outcome]])) - nrow(dd),
          estimate = unname(co[tm, "Estimate"]),
          se = unname(co[tm, "Std. Error"]),
          p_value = unname(co[tm, "Pr(>|t|)"]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  if (verbose && !is.null(out)) {
    msg("  %s: value-cap check over %d settings", outcome, length(caps))
  }
  out
}

#' Compare alternative definitions of the same construct (e.g. the three
#' depth-index denominators), all on the primary QC-passing sample.
sens_alt_definitions <- function(d, outcomes, label, verbose = TRUE) {
  rows <- list()
  for (y in outcomes) {
    if (!y %in% names(d)) next
    for (locus in c("16p", "22q")) {
      s <- locus_subset(d[!is.na(d[[y]]), , drop = FALSE], locus)
      if (nrow(s) < 20 || length(unique(s$CopyNum)) < 3) next
      for (mod in c("mirror", "ushape")) {
        tm  <- if (mod == "mirror") "copy_c" else "copy_q"
        fit <- fit_lmm(build_formula(y, tm, s), s)
        if (is.null(fit)) next
        co <- summary(fit)$coefficients
        if (!tm %in% rownames(co)) next
        rows[[length(rows) + 1]] <- data.frame(
          outcome = y, check = label, setting = y,
          locus = locus, model = mod, term = tm,
          n_nights = nrow(s), n_dropped = NA_integer_,
          estimate = unname(co[tm, "Estimate"]),
          se = unname(co[tm, "Std. Error"]),
          p_value = unname(co[tm, "Pr(>|t|)"]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

#' Re-derive the analysis table under alternative QC thresholds and refit.
sens_qc_thresholds <- function(metrics, link, manifest, outcome,
                               trt_grid = c(30, 60, 120, 180),
                               mt_grid  = c(0.30, 0.50, 0.70),
                               max_pct_n3 = 50,
                               depth_def = "N3_over_others",
                               verbose = TRUE) {
  rows <- list()
  grid <- expand.grid(trt = trt_grid, mt = mt_grid)
  for (i in seq_len(nrow(grid))) {
    mm <- apply_qc(metrics, grid$trt[i], grid$mt[i], TRUE, max_pct_n3,
                   verbose = FALSE)
    d2 <- tryCatch(build_analysis_table(mm, link, manifest, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(d2) || !nrow(d2)) next
    if (identical(depth_def, "N3_over_NREM")) d2$DI <- d2$DI_over_NREM
    if (identical(depth_def, "N3_over_TST"))  d2$DI <- d2$DI_over_TST
    for (locus in c("16p", "22q")) {
      s <- locus_subset(d2[!is.na(d2[[outcome]]), , drop = FALSE], locus)
      if (nrow(s) < 20 || length(unique(s$CopyNum)) < 3) next
      fit <- fit_lmm(build_formula(outcome, "copy_c", s), s)
      if (is.null(fit)) next
      co <- summary(fit)$coefficients
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, check = "QC thresholds",
        setting = sprintf("TRT >= %g min, MT <= %g%%", grid$trt[i], 100 * grid$mt[i]),
        locus = locus, model = "mirror", term = "copy_c",
        n_nights = nrow(s), n_dropped = sum(!mm$qc_pass),
        estimate = unname(co["copy_c", "Estimate"]),
        se = unname(co["copy_c", "Std. Error"]),
        p_value = unname(co["copy_c", "Pr(>|t|)"]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Vary the implausible-architecture QC threshold, including switching it off.
#' This is the direct sensitivity check on the N3 outlier filter now applied
#' as primary QC -- it shows exactly how much the result depends on it.
sens_n3_threshold <- function(metrics, link, manifest, outcome,
                              trt_min = 120, mt_frac = 0.50,
                              n3_grid = c(Inf, 70, 60, 50, 40),
                              verbose = TRUE) {
  rows <- list()
  for (cap in n3_grid) {
    mm <- apply_qc(metrics, trt_min, mt_frac, TRUE, cap, verbose = FALSE)
    d2 <- tryCatch(build_analysis_table(mm, link, manifest, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(d2) || !nrow(d2)) next
    for (locus in c("16p", "22q")) {
      s <- locus_subset(d2[!is.na(d2[[outcome]]), , drop = FALSE], locus)
      if (nrow(s) < 20 || length(unique(s$CopyNum)) < 3) next
      fit <- fit_lmm(build_formula(outcome, "copy_c", s), s)
      if (is.null(fit)) next
      co <- summary(fit)$coefficients
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, check = "N3 outlier threshold",
        setting = if (is.infinite(cap)) "no filter (all nights)"
                  else sprintf("N3 <= %g%% of TST%s", cap,
                               if (cap == 50) "  [primary]" else ""),
        locus = locus, model = "mirror", term = "copy_c",
        n_nights = nrow(s), n_dropped = sum(mm$excl_implausible_n3 & mm$qc_pass_core),
        estimate = unname(co["copy_c", "Estimate"]),
        se = unname(co["copy_c", "Std. Error"]),
        p_value = unname(co["copy_c", "Pr(>|t|)"]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (verbose) msg("  %s: N3 threshold check over %d settings", outcome, length(n3_grid))
  do.call(rbind, rows)
}

#' Run the full sensitivity battery for the primary outcomes.
run_sensitivity <- function(d, metrics, link, manifest, verbose = TRUE,
                            max_pct_n3 = 50) {
  if (verbose) msg("Sensitivity analyses:")
  out <- list(
    sens_n3_threshold(metrics, link, manifest, "DI", verbose = verbose),
    sens_outlier_cap(d, "DI",  c(Inf, 2, 1.5, 1), verbose),
    sens_outlier_cap(d, "SOL", c(Inf, 180, 120, 60), verbose),
    sens_alt_definitions(d, c("DI", "DI_over_NREM", "DI_over_TST"),
                         "depth-index definition", verbose),
    sens_alt_definitions(d, c("SOL", "SOL_log", "SOL_persist"),
                         "onset-latency definition", verbose),
    sens_qc_thresholds(metrics, link, manifest, "DI",
                       max_pct_n3 = max_pct_n3, verbose = verbose)
  )
  res <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (verbose && !is.null(res)) msg("  %d sensitivity fits", nrow(res))
  res
}
