# =====================================================================
# run_all.R -- end-to-end sleep architecture x CNV analysis
#
#   Rscript run_all.R
#
# Everything configurable lives in config.R.
# =====================================================================

PROJ <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(".")
}, error = function(e) normalizePath("."))

source(file.path(PROJ, "config.R"))
for (f in c("00_utils.R", "01_ingest.R", "02_qc_metrics.R",
            "03_models.R", "04_figures.R", "05_sensitivity.R",
            "06_tradeoff.R")) {
  source(file.path(PROJ, "R", f))
}

# resolve config paths: absolute used as-is, relative taken from the project dir
HYPNO_DIR <- resolve_path(HYPNO_DIR, PROJ)
MANIFEST  <- resolve_path(MANIFEST,  PROJ)

OUT <- ensure_dir(resolve_path(OUT_DIR, PROJ))
FIG <- ensure_dir(file.path(OUT, "figures"))
# clear previously generated figures so a rename never leaves a stale file
invisible(file.remove(list.files(FIG, pattern = "[.]png$", full.names = TRUE)))

set.seed(1)

msg("=====================================================================")
msg("Sleep architecture x 16p11.2 / 22q11.2 copy number")
msg("=====================================================================")

# ---------------------------------------------------------------------
# 1. Ingest
# ---------------------------------------------------------------------
msg("\n[1] INGEST")
manifest <- read_manifest(MANIFEST)

epochs <- read_all_hypnograms(
  HYPNO_DIR, epoch_sec = EPOCH_SEC, glob = HYPNO_GLOB,
  stage_col = HYPNO_STAGE_COL, coding = NUMERIC_STAGE_CODING,
  recursive = HYPNO_RECURSIVE
)

link <- link_recordings_to_subjects(
  epochs$recording_id, manifest,
  mode = SUBJECT_ID_MODE, regex = SUBJECT_ID_REGEX,
  device_regex = DEVICE_ID_REGEX
)

# ---------------------------------------------------------------------
# 2. Metrics + QC
# ---------------------------------------------------------------------
msg("\n[2] METRICS AND QC")
metrics <- compute_recording_metrics(epochs, EPOCH_SEC, PERSIST_SLEEP_MIN)
metrics <- apply_qc(metrics, QC_MIN_TRT_MIN, QC_MAX_MT_FRAC, QC_DROP_TST_ZERO,
                    QC_MAX_PCT_N3)

write.csv(metrics, file.path(OUT, "per_recording_metrics_with_qc.csv"), row.names = FALSE)
write.csv(qc_venn_counts(metrics), file.path(OUT, "qc_venn_counts.csv"), row.names = FALSE)
qc_sum <- qc_summary_table(metrics, QC_MIN_TRT_MIN, QC_MAX_MT_FRAC, QC_MAX_PCT_N3)
write.csv(qc_sum, file.path(OUT, "qc_summary.csv"), row.names = FALSE)
plot_qc_venn(metrics, file.path(FIG, "00_qc_venn.png"),
             min_trt_min = QC_MIN_TRT_MIN, max_mt_frac = QC_MAX_MT_FRAC)

d <- build_analysis_table(metrics, link, manifest)
if (!nrow(d)) stop("Analysis table is empty -- check the linkage step.", call. = FALSE)

# `DI` is computed as N3 / (N1+N2+R).  If a different definition is chosen as
# primary, swap it into `DI` so the rest of the pipeline is unchanged.
di_alt <- c(N3_over_NREM = "DI_over_NREM", N3_over_TST = "DI_over_TST")
if (DEPTH_INDEX_DEF %in% names(di_alt)) {
  other <- di_alt[[DEPTH_INDEX_DEF]]
  tmp <- d$DI; d$DI <- d[[other]]; d[[other]] <- tmp
  msg("NOTE: DEPTH_INDEX_DEF = %s -- primary DI swapped with %s.",
      DEPTH_INDEX_DEF, other)
} else if (!identical(DEPTH_INDEX_DEF, "N3_over_others")) {
  stop("Unrecognised DEPTH_INDEX_DEF: ", DEPTH_INDEX_DEF, call. = FALSE)
}
msg("Depth index (primary) = %s;  observed range %.3f - %.3f, median %.3f",
    DEPTH_INDEX_DEF, min(d$DI, na.rm = TRUE), max(d$DI, na.rm = TRUE),
    stats::median(d$DI, na.rm = TRUE))

# ---------------------------------------------------------------------
# 3. Covariate adjustment (ONE reference model per outcome)
# ---------------------------------------------------------------------
msg("\n[3] COVARIATE ADJUSTMENT")
d <- add_adjusted_outcomes(d, OUTCOMES)

write.csv(d, file.path(OUT, "analysis_dataset_per_night.csv"), row.names = FALSE)
desc <- describe_groups(d, OUTCOMES)
write.csv(desc, file.path(OUT, "group_descriptives.csv"), row.names = FALSE)
msg("")
print(desc[, c("Genotype", "n_nights", "n_subjects", "age_mean", "pct_female")])

# ---------------------------------------------------------------------
# 4. Models
# ---------------------------------------------------------------------
msg("\n[4] MODELS")
all_tests <- list()
all_fits  <- list()
for (o in OUTCOMES) {
  if (!o$name %in% names(d)) next
  r <- run_outcome(d, o$name)
  all_tests[[o$name]] <- r$tests
  all_fits[[o$name]]  <- r$fits
}
tests <- do.call(rbind, all_tests[!vapply(all_tests, is.null, logical(1))])
tests <- finalise_results(tests, OUTCOMES, FDR_METHOD)

write.csv(tests, file.path(OUT, "results_all_tests.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4b. Sensitivity analyses
# ---------------------------------------------------------------------
msg("")
sens <- run_sensitivity(d, metrics, link, manifest, max_pct_n3 = QC_MAX_PCT_N3)
write.csv(sens, file.path(OUT, "sensitivity_analyses.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4c. Depth / duration trade-off
# ---------------------------------------------------------------------
msg("")
tr <- run_tradeoff(d, "DI", epochs = epochs)
write.csv(tr$correlations,  file.path(OUT, "tradeoff_correlations.csv"), row.names = FALSE)
write.csv(tr$by_group,      file.path(OUT, "tradeoff_slopes_by_group.csv"), row.names = FALSE)
write.csv(tr$interactions,  file.path(OUT, "tradeoff_group_comparisons.csv"), row.names = FALSE)
write.csv(tr$loglog,        file.path(OUT, "tradeoff_loglog_check.csv"), row.names = FALSE)
write.csv(tr$null_calibration, file.path(OUT, "tradeoff_null_calibration.csv"), row.names = FALSE)
write.csv(tr$adjusted_main, file.path(OUT, "tradeoff_main_effect_adjusted_for_TST.csv"), row.names = FALSE)
if (!is.null(tr$fixed_correlations))
  write.csv(tr$fixed_correlations, file.path(OUT, "tradeoff_fixedwindow_correlations.csv"),
            row.names = FALSE)
if (!is.null(tr$fixed_by_group))
  write.csv(tr$fixed_by_group, file.path(OUT, "tradeoff_fixedwindow_slopes.csv"),
            row.names = FALSE)

if (requireNamespace("writexl", quietly = TRUE)) {
  sheets <- list(all_tests = tests, descriptives = desc, sensitivity = sens,
                 qc_summary = qc_sum, qc_venn = qc_venn_counts(metrics),
                 tradeoff_corr = tr$correlations,
                 tradeoff_slopes = tr$by_group,
                 tradeoff_groupdiff = tr$interactions,
                 tradeoff_loglog = tr$loglog,
                 tradeoff_null = tr$null_calibration,
                 tradeoff_adjmain = tr$adjusted_main)
  if (!is.null(tr$fixed_correlations))
    sheets$tradeoff_fixedwin <- tr$fixed_correlations
  writexl::write_xlsx(sheets, file.path(OUT, "results_all_tests.xlsx"))
}

msg("\nPrimary results (mirror and U-shape, single-term models):")
key <- tests[tests$primary_outcome & tests$model %in% c("mirror_single", "ushape_single"), ]
print(data.frame(
  outcome = key$outcome, locus = key$cohort, model = key$model,
  estimate = signif(key$estimate, 3), p = fmt_p(key$p_value), q = fmt_p(key$q_value)
), row.names = FALSE)

# ---------------------------------------------------------------------
# 5. Figures
# ---------------------------------------------------------------------
msg("\n[5] FIGURES")
i <- 0
for (o in OUTCOMES) {
  fits <- all_fits[[o$name]]
  if (is.null(fits)) next
  i <- i + 1
  tag <- sprintf("%02d_%s", i, o$name)
  save_fig(plot_five_group(d, fits, o$name, o$label, tests),
           file.path(FIG, paste0(tag, "_A_five_group.png")))
  for (locus in c("16p", "22q")) {
    save_fig(plot_locus(d, fits, o$name, o$label, tests, locus, "joint"),
             file.path(FIG, paste0(tag, "_B_", locus, "_joint.png")))
    save_fig(plot_locus(d, fits, o$name, o$label, tests, locus, "mirror"),
             file.path(FIG, paste0(tag, "_C_", locus, "_mirror.png")))
  }
  msg("  figures for %s", o$name)
}
save_fig(plot_copy_number_lines(all_fits, OUTCOMES, tests),
         file.path(FIG, "99_copy_number_lines.png"), w = 12, h = 4.6)

save_fig(plot_tradeoff_scatter(d, tr$correlations, "DI", "NREM depth index"),
         file.path(FIG, "90_tradeoff_raw_scatter.png"), w = 16, h = 5.0)
save_fig(plot_tradeoff_within(d, "DI", "NREM depth index"),
         file.path(FIG, "91_tradeoff_within_participant.png"), w = 14, h = 5.2)
save_fig(plot_tradeoff_slopes(tr$by_group, "NREM depth index"),
         file.path(FIG, "92_tradeoff_slopes_by_group.png"), w = 11, h = 5.2)

# copy-number effect with sleep duration removed
d <- add_tst_adjusted_outcome(d, "DI")
for (locus in c("16p", "22q")) {
  save_fig(plot_locus_tst_adjusted(d, "DI", "NREM depth index", locus),
           file.path(FIG, paste0("93_DI_", locus, "_TSTadjusted.png")))
}
msg("  trade-off figures")

saveRDS(list(d = d, tests = tests, sens = sens, tradeoff = tr, fits = all_fits,
             metrics = metrics, link = link, manifest = manifest),
        file.path(OUT, "workspace.rds"))

msg("\nDONE. Outputs in: %s", OUT)
msg("  results_all_tests.csv / .xlsx   -- every statistical test")
msg("  analysis_dataset_per_night.csv  -- the cleaned per-night table")
msg("  group_descriptives.csv          -- N, age, sex by group")
msg("  figures/                        -- QC Venn, box/scatter panels, line plot")
