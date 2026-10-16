# =====================================================================
# 02_qc_metrics.R -- per-recording sleep metrics and QC
#
# Metric definitions (all in minutes unless noted):
#
#   TRT   total recording time  = all epochs x epoch length
#   TST   total sleep time      = (N1 + N2 + N3 + R) epochs
#   SE    sleep efficiency      = TST / TRT
#   WASO  wake after sleep onset
#   DI    NREM depth index      = N3 / (N1 + N2 + R)           [primary]
#                                 i.e. N3 over the sum of the OTHER sleep
#                                 stages -- a ratio, not a proportion, so it
#                                 is unbounded above.
#   DI_over_NREM                = N3 / (N1 + N2 + N3)          [sensitivity]
#   DI_over_TST                 = N3 / (N1 + N2 + N3 + R)      [sensitivity]
#   SOL   sleep onset latency   = recording start -> first sleep epoch
#         (no lights-out diary is available, so the recording start is used)
#   SOL_persist                 = start -> onset of PERSIST_SLEEP_MIN of
#                                 near-continuous sleep                [sensitivity]
#
# QC exclusions, matching the original study's QC:
#   (a) TST == 0                 "No Sleep"
#   (b) TRT < QC_MIN_TRT_MIN     "Too Short"
#   (c) MT fraction > 50%        "High Movement"
# =====================================================================

SLEEP_STAGES <- c("N1", "N2", "N3", "R")
NREM_STAGES  <- c("N1", "N2", "N3")

#' Compute per-recording sleep metrics from the tidy epoch table.
compute_recording_metrics <- function(epochs, epoch_sec = 30,
                                      persist_min = 10) {
  epoch_min <- epoch_sec / 60

  # Per-epoch durations, when the source provides them (Dreem does), are used
  # for all time totals rather than assuming a fixed epoch length.
  if (!"dur_sec" %in% names(epochs)) epochs$dur_sec <- epoch_sec

  epochs <- epochs[order(epochs$recording_id, epochs$epoch_idx), , drop = FALSE]
  split_ep  <- split(epochs$stage,   epochs$recording_id)
  split_dur <- split(epochs$dur_sec, epochs$recording_id)

  n_persist <- max(1L, ceiling(persist_min / epoch_min))

  res <- lapply(names(split_ep), function(rid) {
    s <- split_ep[[rid]]
    dmin <- split_dur[[rid]] / 60
    n <- length(s)
    stage_chr <- as.character(s)
    cnt <- table(factor(stage_chr, levels = CANON_STAGES))
    # minutes per stage, from the actual epoch durations
    mins <- vapply(CANON_STAGES, function(st) sum(dmin[stage_chr == st]), numeric(1))

    n_sleep <- sum(cnt[SLEEP_STAGES])
    n_nrem  <- sum(cnt[NREM_STAGES])
    min_sleep <- sum(mins[SLEEP_STAGES])

    is_sleep <- stage_chr %in% SLEEP_STAGES
    first_sleep <- which(is_sleep)[1]

    # latency to persistent sleep: first index starting a window of
    # n_persist epochs that is >=90% sleep
    sol_persist <- NA_real_
    if (n >= n_persist && any(is_sleep)) {
      run <- stats::filter(as.numeric(is_sleep), rep(1, n_persist), sides = 1)
      ok <- which(run >= 0.9 * n_persist)
      if (length(ok)) sol_persist <- (ok[1] - n_persist) * epoch_min
    }

    data.frame(
      recording_id = rid,
      n_epochs     = n,
      TRT          = sum(dmin),
      TST          = min_sleep,
      WASO         = if (!is.na(first_sleep))
                       sum(dmin[first_sleep:n][stage_chr[first_sleep:n] == "W"])
                     else NA_real_,
      SOL          = if (!is.na(first_sleep) && first_sleep > 1)
                       sum(dmin[seq_len(first_sleep - 1)])
                     else if (!is.na(first_sleep)) 0 else NA_real_,
      SOL_persist  = sol_persist,
      min_N1  = mins[["N1"]],
      min_N2  = mins[["N2"]],
      min_N3  = mins[["N3"]],
      min_R   = mins[["R"]],
      min_W   = mins[["W"]],
      min_MT  = mins[["MT"]],
      min_UNSCORED = mins[["UNSCORED"]],
      frac_MT = as.numeric(cnt[["MT"]]) / n,
      SE      = if (sum(dmin) > 0) min_sleep / sum(dmin) else NA_real_,
      # primary: N3 relative to the sum of the other sleep stages
      DI           = if (min_sleep - mins[["N3"]] > 0)
                       mins[["N3"]] / (min_sleep - mins[["N3"]]) else NA_real_,
      DI_over_NREM = if (sum(mins[NREM_STAGES]) > 0)
                       mins[["N3"]] / sum(mins[NREM_STAGES]) else NA_real_,
      DI_over_TST  = if (min_sleep > 0) mins[["N3"]] / min_sleep else NA_real_,
      pct_N1 = if (min_sleep > 0) 100 * mins[["N1"]] / min_sleep else NA_real_,
      pct_N2 = if (min_sleep > 0) 100 * mins[["N2"]] / min_sleep else NA_real_,
      pct_N3 = if (min_sleep > 0) 100 * mins[["N3"]] / min_sleep else NA_real_,
      pct_R  = if (min_sleep > 0) 100 * mins[["R"]]  / min_sleep else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  m <- do.call(rbind, res)
  m$SOL_log <- log1p(m$SOL)
  rownames(m) <- NULL
  m
}

#' Flag QC exclusions.  Returns `metrics` with three logical columns plus
#' `qc_pass`, in the same order as the Venn diagram.
apply_qc <- function(metrics, min_trt_min = 120, max_mt_frac = 0.50,
                     drop_tst_zero = TRUE, max_pct_n3 = Inf, verbose = TRUE) {
  # --- the three original criteria (these are what the Venn diagram shows) ---
  metrics$excl_no_sleep     <- drop_tst_zero & (metrics$TST == 0 | is.na(metrics$TST))
  metrics$excl_too_short    <- metrics$TRT < min_trt_min
  metrics$excl_high_movement <- !is.na(metrics$frac_MT) & metrics$frac_MT > max_mt_frac

  metrics$qc_pass_core <- !(metrics$excl_no_sleep |
                            metrics$excl_too_short |
                            metrics$excl_high_movement)

  # --- implausible sleep architecture, applied after the core three ---------
  metrics$excl_implausible_n3 <- is.finite(max_pct_n3) &
    !is.na(metrics$pct_N3) & metrics$pct_N3 > max_pct_n3

  metrics$qc_pass <- metrics$qc_pass_core & !metrics$excl_implausible_n3

  if (verbose) {
    n <- nrow(metrics)
    msg("QC on %d recordings", n)
    msg("  No Sleep (TST = 0)           : %d", sum(metrics$excl_no_sleep))
    msg("  Too Short (TRT < %g min)      : %d", min_trt_min, sum(metrics$excl_too_short))
    msg("  High Movement (MT > %g%%)      : %d", 100 * max_mt_frac,
        sum(metrics$excl_high_movement))
    msg("  ---------------------------------")
    msg("  excluded by the three (unique): %d", sum(!metrics$qc_pass_core))
    msg("  pass the three               : %d", sum(metrics$qc_pass_core))
    if (is.finite(max_pct_n3)) {
      msg("  Implausible N3 (> %g%% of TST) : %d additional",
          max_pct_n3, sum(metrics$excl_implausible_n3 & metrics$qc_pass_core))
    }
    msg("  ---------------------------------")
    msg("  PASS                         : %d", sum(metrics$qc_pass))
  }
  metrics
}

#' Sequential QC summary: how many recordings each filter removes, in order.
qc_summary_table <- function(metrics, min_trt_min, max_mt_frac, max_pct_n3) {
  core   <- metrics$qc_pass_core
  n3_on  <- is.finite(max_pct_n3)
  n3_extra <- if (n3_on) sum(metrics$excl_implausible_n3 & core) else 0L

  add <- function(step, excluded = NA, remaining = NA, note = "") {
    data.frame(step = step, excluded = excluded, remaining = remaining,
               note = note, stringsAsFactors = FALSE)
  }

  rows <- list(
    add("Hypnogram files parsed", NA, nrow(metrics)),
    add("No sleep (TST = 0)", sum(metrics$excl_no_sleep), NA,
        "the three criteria overlap, so they do not sum to the unique total"),
    add(sprintf("Too short (TRT < %g min)", min_trt_min),
        sum(metrics$excl_too_short)),
    add(sprintf("High movement (MT > %g%% of epochs)", 100 * max_mt_frac),
        sum(metrics$excl_high_movement))
  )

  if (n3_on) {
    rows <- c(rows, list(
      add("= passing the three criteria above", sum(!core), sum(core),
          "matches the original analysis"),
      add(sprintf("Implausible architecture (N3 > %g%% of TST)", max_pct_n3),
          n3_extra, NA, "physiological bound; see the sensitivity table"),
      add("= passing all QC", sum(!metrics$qc_pass), sum(metrics$qc_pass))))
  } else {
    rows <- c(rows, list(
      add("= passing all QC", sum(!metrics$qc_pass), sum(metrics$qc_pass),
          "matches the original analysis; no outlier filter applied")))
  }
  do.call(rbind, rows)
}

#' Counts for each of the 7 Venn regions, in the order plotted.
qc_venn_counts <- function(metrics) {
  a <- metrics$excl_no_sleep
  b <- metrics$excl_too_short
  c <- metrics$excl_high_movement
  data.frame(
    region = c("NoSleep only", "TooShort only", "HighMove only",
               "NoSleep & TooShort", "NoSleep & HighMove",
               "TooShort & HighMove", "all three"),
    n = c(sum(a & !b & !c), sum(!a & b & !c), sum(!a & !b & c),
          sum(a & b & !c), sum(a & !b & c), sum(!a & b & c), sum(a & b & c)),
    stringsAsFactors = FALSE
  )
}

#' Build the final per-night analysis table:
#' QC-passing recordings joined to subject-level covariates.
build_analysis_table <- function(metrics, link, manifest, verbose = TRUE) {
  d <- merge(metrics, link, by = "recording_id", all.x = TRUE)
  n0 <- nrow(d)

  d <- d[d$qc_pass, , drop = FALSE]
  n_qc <- nrow(d)

  d <- d[!is.na(d$GUID), , drop = FALSE]
  n_linked <- nrow(d)

  d <- merge(d, manifest, by = "GUID", all.x = FALSE)
  n_mani <- nrow(d)

  # complete covariates -- required by the mixed models
  ok <- !is.na(d$Age) & !is.na(d$Sex) & nzchar(d$Sex) &
        !is.na(d$Site) & nzchar(d$Site)
  d <- d[ok, , drop = FALSE]

  if (verbose) {
    msg("Analysis table assembly:")
    msg("  recordings parsed                : %d", n0)
    msg("  after QC                         : %d", n_qc)
    msg("  linked to a manifest subject     : %d", n_linked)
    msg("  in one of the 5 analysis groups  : %d", n_mani)
    msg("  with complete age/sex/site       : %d", nrow(d))
    msg("  subjects                         : %d", length(unique(d$GUID)))
  }

  d$GUID     <- factor(d$GUID)
  d$Sex      <- factor(d$Sex)
  d$Site     <- factor(d$Site)
  d$Genotype <- factor(as.character(d$Genotype), levels = GENO_LEVELS)
  d$Genotype <- stats::relevel(d$Genotype, ref = "Control")

  # Age is mean-centred before squaring: this leaves the model span (and so
  # every fitted value and genotype coefficient) unchanged, but removes the
  # near-collinearity between Age and Age^2 and stops lme4 warning that the
  # predictors are on very different scales.
  d$Age_c <- d$Age - mean(d$Age, na.rm = TRUE)
  attr(d, "age_mean") <- mean(d$Age, na.rm = TRUE)

  # copy-number contrasts
  #   copy_c : mirror / linear   Del = -1, Control = 0, Dup = +1
  #   copy_q : U-shape           Del = +1, Control = 0, Dup = +1
  d$copy_c <- d$CopyNum - 2L
  d$copy_q <- d$copy_c^2

  d[order(d$GUID, d$recording_id), , drop = FALSE]
}

#' Per-group descriptive summary (nights, subjects, age, sex, outcome means).
describe_groups <- function(d, outcomes) {
  vars <- vapply(outcomes, function(o) o$name, character(1))
  vars <- intersect(vars, names(d))
  out <- lapply(levels(d$Genotype), function(g) {
    s <- d[d$Genotype == g, , drop = FALSE]
    base <- data.frame(
      Genotype   = g,
      n_nights   = nrow(s),
      n_subjects = length(unique(as.character(s$GUID))),
      nights_per_subject = if (nrow(s)) round(nrow(s) / length(unique(as.character(s$GUID))), 2) else NA,
      age_mean   = round(mean(s$Age, na.rm = TRUE), 2),
      age_sd     = round(stats::sd(s$Age, na.rm = TRUE), 2),
      pct_female = round(100 * mean(grepl("^F", s$Sex, ignore.case = TRUE)), 1),
      stringsAsFactors = FALSE
    )
    for (v in vars) {
      base[[paste0(v, "_mean")]] <- round(mean(s[[v]], na.rm = TRUE), 3)
      base[[paste0(v, "_sd")]]   <- round(stats::sd(s[[v]], na.rm = TRUE), 3)
    }
    base
  })
  do.call(rbind, out)
}
