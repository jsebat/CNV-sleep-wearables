# =====================================================================
# tests/make_synthetic_data.R
#
# Generates a fake hypnogram directory + sample manifest with the same
# SHAPE as the real inputs, so the pipeline can be exercised end to end
# without the real data.  Used to verify the code, not to produce results.
#
#   Rscript tests/make_synthetic_data.R <outdir>
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else file.path(tempdir(), "synth")
hyp <- file.path(outdir, "sfari_hypnograms")
dir.create(hyp, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

GROUPS <- c("16p11.2 deletion", "16p11.2 duplication",
            "22q11.2 deletion", "22q11.2 duplication", "Control")
N_SUBJ <- c(38, 34, 40, 30, 95)
SITES  <- c("UCSD", "Vanderbilt", "UCLA", "Penn")

# true effects, on the NREM depth index scale
mirror_di <- c("16p11.2 deletion" = -0.045, "16p11.2 duplication" =  0.045,
               "22q11.2 deletion" =  0.040, "22q11.2 duplication" = -0.040,
               "Control" = 0)
# Total sleep time.  These numbers exist ONLY to give the pipeline a known
# signal to recover -- they are not estimates and carry no scientific claim.
# Directions are set to the expected pattern (16p11.2 deletion sleeps longer)
# so the demo figures are not misread as results.
ushape_tst <- c("16p11.2 deletion" =  35, "16p11.2 duplication" = -15,
                "22q11.2 deletion" = -20, "22q11.2 duplication" = -10,
                "Control" = 0)

subj <- do.call(rbind, lapply(seq_along(GROUPS), function(g) {
  n <- N_SUBJ[g]
  data.frame(
    GUID = sprintf("NDAR%s%03d", LETTERS[g], seq_len(n)),
    Diagnosis = GROUPS[g],
    Age = round(runif(n, 6, 42), 1),
    Sex = sample(c("M", "F"), n, TRUE, c(0.55, 0.45)),
    Site = sample(SITES, n, TRUE),
    stringsAsFactors = FALSE
  )
}))
# a few manifest rows the analysis should ignore
subj <- rbind(subj, data.frame(
  GUID = sprintf("NDARZ%03d", 1:12), Diagnosis = "1q21.1 deletion",
  Age = round(runif(12, 8, 30), 1), Sex = "F", Site = "UCSD",
  stringsAsFactors = FALSE))

subj$cohort_wave <- sample(c("w1", "w2"), nrow(subj), TRUE)   # decoy column
write.csv(subj, file.path(outdir, "full_sample_manifest.csv"), row.names = FALSE)

# subject random intercepts
u_di  <- setNames(rnorm(nrow(subj), 0, 0.055), subj$GUID)
u_tst <- setNames(rnorm(nrow(subj), 0, 42),    subj$GUID)

site_di  <- setNames(c(0, 0.02, -0.015, 0.008), SITES)
site_tst <- setNames(c(0, 15, -10, 6), SITES)

n_rec <- 0
for (i in seq_len(nrow(subj))) {
  s <- subj[i, ]
  if (!s$Diagnosis %in% GROUPS) next
  n_nights <- sample(1:6, 1, prob = c(.10, .18, .24, .22, .16, .10))

  for (night in seq_len(n_nights)) {
    n_rec <- n_rec + 1

    # ---- occasional QC failures, to populate the Venn ----
    fail <- sample(c("none", "short", "movement", "nosleep"), 1,
                   prob = c(0.885, 0.020, 0.085, 0.010))

    trt_min <- if (fail == "short") runif(1, 20, 110) else rnorm(1, 505, 55)
    trt_min <- max(trt_min, 12)
    n_ep <- max(round(trt_min * 2), 24)          # 30 s epochs

    # target metrics for this night
    di <- 0.32 + mirror_di[[s$Diagnosis]] + u_di[[s$GUID]] +
      site_di[[s$Site]] - 0.0025 * (s$Age - 20) + rnorm(1, 0, 0.06)
    di <- min(max(di, 0.02), 0.85)

    tst_min <- 430 + ushape_tst[[s$Diagnosis]] + u_tst[[s$GUID]] +
      site_tst[[s$Site]] - 1.6 * (s$Age - 20) + rnorm(1, 0, 45)
    tst_min <- min(max(tst_min, 30), trt_min * 0.97)
    if (fail == "nosleep") tst_min <- 0

    mt_frac <- if (fail == "movement") runif(1, 0.52, 0.85) else rbeta(1, 1.2, 40)
    sol_min <- min(rgamma(1, shape = 1.6, scale = 14), trt_min * 0.4)

    n_mt    <- round(n_ep * mt_frac)
    n_sol   <- round(sol_min * 2)
    n_sleep <- round(tst_min * 2)
    n_rem   <- round(n_sleep * runif(1, 0.15, 0.26))
    n_nrem  <- n_sleep - n_rem
    n_n3    <- round(n_nrem * di)
    n_n1    <- round((n_nrem - n_n3) * runif(1, 0.12, 0.22))
    n_n2    <- max(n_nrem - n_n3 - n_n1, 0)
    n_wake  <- max(n_ep - n_mt - n_sol - n_sleep, 0)

    body <- sample(c(rep("N1", n_n1), rep("N2", n_n2), rep("N3", n_n3),
                     rep("R", n_rem), rep("W", n_wake), rep("MT", n_mt)))
    stages <- c(rep("W", n_sol), body)
    stages <- stages[seq_len(min(length(stages), n_ep))]
    if (length(stages) < n_ep) stages <- c(stages, rep("W", n_ep - length(stages)))

    fn <- file.path(hyp, sprintf("%s_night%d_hypnogram.txt", s$GUID, night))
    write.table(
      data.frame(epoch = seq_along(stages),
                 clock_time = format(as.POSIXct("2024-01-01 22:00:00", tz = "UTC") +
                                       (seq_along(stages) - 1) * 30, "%H:%M:%S"),
                 stage = stages),
      fn, sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}

cat(sprintf("Synthetic data written to %s\n  %d recordings, %d manifest rows\n",
            outdir, n_rec, nrow(subj)))
