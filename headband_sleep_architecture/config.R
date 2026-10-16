# =====================================================================
# config.R -- the ONLY file you should normally need to edit
# Sleep architecture x 16p11.2 / 22q11.2 copy number
# =====================================================================

# ---------------------------------------------------------------------
# 1. INPUT PATHS
# ---------------------------------------------------------------------
# DATA_ROOT holds the raw inputs and receives the outputs.  It is deliberately
# OUTSIDE the Dropbox-synced tree: the extracted hypnograms are ~900 small
# files (~75 MB) that are simply a re-extraction of sfari_hypnograms.zip and
# demographics_distributions.zip, and syncing them serves no purpose.
#
# This is the one machine-specific line in the repo.  Override it without
# editing the file by setting the SLEEP_HYPNO_DIR / SLEEP_MANIFEST /
# SLEEP_OUT_DIR environment variables (see the bottom of this section).
#
# Expected layout under DATA_ROOT:
#   data/sfari_hypnograms/                     <- from sfari_hypnograms.zip
#   data/demographics_distributions/           <- from demographics_distributions.zip
#   results/                                   <- created by run_all.R
DATA_ROOT <- "C:/Users/Jonathan Sebat/sleep_arch_data"

# Paths are resolved against the project directory if relative, used as-is if
# absolute.
HYPNO_DIR <- file.path(DATA_ROOT, "data/sfari_hypnograms")
MANIFEST  <- file.path(DATA_ROOT, "data/demographics_distributions/full_sample_manifest.csv")

# Output directory.  Also outside Dropbox by default -- run_all.R writes ~25 MB
# here on every run, including a 9 MB workspace.rds.  The curated deliverables
# live in G2MH/SleepArchitectureResults instead.
OUT_DIR   <- file.path(DATA_ROOT, "results")

# Environment-variable overrides (used by the test harness; ignore otherwise).
if (nzchar(Sys.getenv("SLEEP_HYPNO_DIR"))) HYPNO_DIR <- Sys.getenv("SLEEP_HYPNO_DIR")
if (nzchar(Sys.getenv("SLEEP_MANIFEST")))  MANIFEST  <- Sys.getenv("SLEEP_MANIFEST")
if (nzchar(Sys.getenv("SLEEP_OUT_DIR")))   OUT_DIR   <- Sys.getenv("SLEEP_OUT_DIR")

# ---------------------------------------------------------------------
# 2. HYPNOGRAM FORMAT
#    Leave as NULL to auto-detect; set explicitly to override.
# ---------------------------------------------------------------------
EPOCH_SEC        <- 30      # seconds per scored epoch
HYPNO_GLOB       <- NULL    # e.g. "*.txt"; NULL = any of txt/csv/tsv/dat/stg
HYPNO_STAGE_COL  <- NULL    # e.g. "stage"; NULL = auto-detect
HYPNO_RECURSIVE  <- TRUE

# Numeric stage coding convention.
#   "RK"   : 0=W 1=N1 2=N2 3=N3 4=N3(S4) 5=REM 6=MT   (Rechtschaffen & Kales)
#   "AASM" : 0=W 1=N1 2=N2 3=N3 4=REM
NUMERIC_STAGE_CODING <- "RK"

# How to get a subject ID from a hypnogram filename.
#   "device_id"      : filenames carry a Dreem device ID, looked up in the
#                      manifest's `dreem_id` column.  <- this cohort
#   "match_manifest" : find the manifest GUID inside the filename.
#   "regex"          : use SUBJECT_ID_REGEX below, capture group 1.
SUBJECT_ID_MODE  <- "device_id"
SUBJECT_ID_REGEX <- "^([^_.]+)"
# Dreem: SFARI_Dreem_Prod_Hypno_SITE0000_2025-01-01T23_00_00-05_00__1.txt
DEVICE_ID_REGEX  <- "^SFARI_Dreem_Prod_Hypno_([A-Za-z]+[0-9]+)_"

# ---------------------------------------------------------------------
# 3. QC THRESHOLDS
#    Defaults reproduce the original study's QC (881 -> 782 recordings).
# ---------------------------------------------------------------------
QC_MIN_TRT_MIN <- 120    # exclude recordings with total recording time < this
                         #   NOTE: the Venn diagram used 2 h (=120 min).
                         #   Set to 30 to use the 30-minute rule instead.
QC_MAX_MT_FRAC <- 0.50   # exclude recordings with > 50% movement-time epochs
QC_DROP_TST_ZERO <- TRUE # exclude recordings with no sleep at all

# Implausible sleep architecture.  N3 does not exceed half of total sleep time
# in real human sleep at any age (typical adult 10-25%, children up to ~35%),
# so a night scored above this is an artefact rather than deep sleep.  It also
# matters mechanically: the depth index N3/(TST-N3) is an unbounded ratio, so
# one such night carries enormous leverage.
#
# NOT APPLIED BY DEFAULT.  Extreme nights are kept in the analysis -- that is
# the more conservative choice, and the figures handle them by clipping the
# axis (see FIG_CLIP_IQR below) rather than by dropping data.
# R/05_sensitivity.R reports what applying this filter would do, at every
# threshold from 40% to 70%, in `sensitivity_analyses.csv`.
# Set to e.g. 50 to apply it.
QC_MAX_PCT_N3 <- Inf

# ---------------------------------------------------------------------
# 4. OUTCOME DEFINITIONS
# ---------------------------------------------------------------------
# NREM depth index.  Primary = N3 over the sum of the OTHER sleep stages.
# All three are always computed; this only chooses which one is called "DI"
# and gets the primary (FDR-corrected) tests.
#   "N3_over_others" : N3 / (N1 + N2 + R)       <- primary; a ratio, unbounded
#   "N3_over_NREM"   : N3 / (N1 + N2 + N3)      <- sensitivity
#   "N3_over_TST"    : N3 / (N1 + N2 + N3 + R)  <- sensitivity
DEPTH_INDEX_DEF <- "N3_over_others"

# Sleep onset latency: no lights-out diary, so measured from recording start.
PERSIST_SLEEP_MIN <- 10   # minutes of ~continuous sleep defining "persistent
                          # sleep" for the secondary SOL_persist outcome

# ---------------------------------------------------------------------
# 5. MODEL SPECIFICATION
# ---------------------------------------------------------------------
# Age_c is age centred on the analysis-sample mean (built in 02_qc_metrics.R).
# Centring before squaring is numerically equivalent to Age + I(Age^2) but
# avoids the Age/Age^2 collinearity warning from lme4.
COVARIATES   <- c("Age_c", "I(Age_c^2)", "Sex", "Site")
RANDOM_TERM  <- "(1 | GUID)"

# Outcomes to run.  name = column in the analysis table, label = axis label.
OUTCOMES <- list(
  list(name = "DI",           label = "NREM Depth Index (N3 / [N1+N2+REM])", primary = TRUE),
  list(name = "TST",          label = "Total Sleep Time (min)",        primary = TRUE),
  list(name = "SOL",          label = "Sleep Onset Latency (min)",     primary = TRUE),
  list(name = "DI_over_NREM", label = "N3 / NREM",                     primary = FALSE),
  list(name = "DI_over_TST",  label = "N3 / TST",                      primary = FALSE),
  list(name = "SOL_log",      label = "log(1 + Sleep Onset Latency)",  primary = FALSE),
  list(name = "SOL_persist",  label = "Latency to Persistent Sleep (min)", primary = FALSE)
)

# Multiple-testing correction is applied within the set of primary tests.
FDR_METHOD <- "BH"

# ---------------------------------------------------------------------
# 6. PLOTTING
# ---------------------------------------------------------------------
FIG_WIDTH  <- 9
FIG_HEIGHT <- 7
FIG_DPI    <- 400

# Y-axis clipping.  A handful of nights sit far outside the bulk of the data
# (the depth index is an unbounded ratio, and sleep onset latency is strongly
# right-skewed), which compresses every panel into a thin band.  Panels are
# therefore drawn with the axis limited to
#     [Q1 - FIG_CLIP_IQR*IQR , Q3 + FIG_CLIP_IQR*IQR]
# via coord_cartesian(), which hides those points WITHOUT removing them: the
# boxplots, the model fit and every p-value still use all the data.  Each
# caption states how many points fall outside the axis.
# Set to Inf to show the full range.
FIG_CLIP_IQR <- 3
