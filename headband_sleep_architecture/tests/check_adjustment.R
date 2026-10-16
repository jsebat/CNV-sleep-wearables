# =====================================================================
# tests/check_adjustment.R
#
# Verifies the properties the covariate adjustment is supposed to have,
# and contrasts them with the two adjustments used in an earlier version of this analysis.
#
#   Rscript tests/check_adjustment.R <results_dir>
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
PROJ <- normalizePath(file.path(dirname(sub("^--file=", "",
          commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), ".."))
RES  <- if (length(args)) args[1] else file.path(PROJ, "results")

source(file.path(PROJ, "config.R"))
for (f in c("00_utils.R", "01_ingest.R", "02_qc_metrics.R", "03_models.R"))
  source(file.path(PROJ, "R", f))

w <- readRDS(file.path(RES, "workspace.rds"))
d <- w$d; fits <- w$fits

ok <- TRUE
check <- function(label, pass, detail = "") {
  ok <<- ok && pass
  msg("  [%s] %s%s", if (pass) "PASS" else "FAIL", label,
      if (nzchar(detail)) paste0("  -- ", detail) else "")
}

for (y in c("DI", "TST", "SOL")) {
  adj <- paste0("adj_", y)
  msg("\n=== %s ===", y)
  s <- d[!is.na(d[[adj]]), ]

  # ---- 1. mean is preserved -------------------------------------
  check("mean(y_adj) == mean(y)",
        isTRUE(all.equal(mean(s[[adj]]), mean(s[[y]]), tolerance = 1e-8)),
        sprintf("%.4f vs %.4f", mean(s[[adj]]), mean(s[[y]])))

  # ---- 2. genotype differences SURVIVE the adjustment ------------
  raw_sd <- stats::sd(tapply(s[[y]],   s$Genotype, mean))
  adj_sd <- stats::sd(tapply(s[[adj]], s$Genotype, mean))
  check("group means still differ after adjustment",
        adj_sd > 0.4 * raw_sd,
        sprintf("sd of group means: raw %.4f -> adjusted %.4f", raw_sd, adj_sd))

  # ---- 3. the adjusted value is a property of the NIGHT ----------
  # The same night must plot at the same height whichever panel it is in.
  s16 <- fits[[y]][["16p_data"]]
  key <- match(s16$recording_id, s$recording_id)
  check("identical value in the 5-group and 16p panels",
        isTRUE(all.equal(s16[[adj]], s[[adj]][key])),
        sprintf("%d shared nights, max diff %.2e", sum(!is.na(key)),
                max(abs(s16[[adj]] - s[[adj]][key]), na.rm = TRUE)))

  # ...and to show what a per-panel refit would have done, redo
  # the adjustment on the 16p subset alone and compare.
  f16  <- build_formula(y, "Genotype", s16)
  m16  <- fit_lmm(f16, s16)
  a16  <- covariate_adjusted(m16, "Genotype", y)
  refit_gap <- max(abs(a16$adj - s16[[adj]][match(a16$.row, rownames(s16))]),
                   na.rm = TRUE)
  msg("       (refitting the adjustment on the 16p subset alone would move")
  msg("        nights by up to %.4f %s -- %.1f%% of an SD)",
      refit_gap, if (y == "TST") "min" else "units",
      100 * refit_gap / stats::sd(s16[[adj]], na.rm = TRUE))

  # ---- 4. adjustment is orthogonal-ish to the covariates ---------
  r_age_raw <- abs(stats::cor(s[[y]],   s$Age, use = "complete.obs"))
  r_age_adj <- abs(stats::cor(s[[adj]], s$Age, use = "complete.obs"))
  check("age association reduced",
        r_age_adj <= r_age_raw + 1e-9,
        sprintf("|r| with Age: %.3f -> %.3f", r_age_raw, r_age_adj))

  # ---- 5. the student's method, for comparison -------------------
  f <- build_formula(y, "Genotype", s)
  m <- fit_lmm(f, s)
  student <- s[[y]] - stats::predict(m, re.form = NA) + mean(s[[y]])
  st_sd <- stats::sd(tapply(student, s$Genotype, mean))
  msg("  ---- earlier version: 'y - predict(re.form=NA) + const' -------------")
  msg("       sd of group means: %.5f  (ours %.5f, raw %.5f)",
      st_sd, adj_sd, raw_sd)
  msg("       -> that method shrinks the between-group signal to %.0f%% of raw;",
      100 * st_sd / raw_sd)
  msg("          ours retains %.0f%%.", 100 * adj_sd / raw_sd)

  # ---- 6. model-implied means line up with the adjusted dots -----
  fit16 <- fits[[y]][["16p_group"]]
  if (!is.null(fit16)) {
    gcols <- intersect(paste0("Genotype", c("16pDel", "16pDup")),
                       colnames(lme4::getME(fit16, "X")))
    tg <- as.data.frame(matrix(0, 3, length(gcols))); names(tg) <- gcols
    tg[1, paste0("Genotype", "16pDel")] <- 1
    tg[3, paste0("Genotype", "16pDup")] <- 1
    mm <- model_means(fit16, "Genotype", tg)
    obs <- tapply(s16[[adj]], factor(as.character(s16$Genotype),
                  levels = c("16pDel", "Control", "16pDup")), mean)
    gap <- max(abs(mm$estimate - as.numeric(obs)))
    check("model line sits within the adjusted dots",
          gap < 0.5 * stats::sd(s16[[adj]], na.rm = TRUE),
          sprintf("max |EMM - mean(adjusted)| = %.4f (sd = %.4f)",
                  gap, stats::sd(s16[[adj]], na.rm = TRUE)))
  }
}

msg("\n%s", if (ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***")
if (!ok) quit(status = 1)
