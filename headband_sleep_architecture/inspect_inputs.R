# =====================================================================
# inspect_inputs.R -- run this FIRST, before run_all.R
#
#   Rscript inspect_inputs.R
#
# Prints what the pipeline sees in your two input sources, so any format
# mismatch shows up immediately instead of halfway through a model fit.
# Nothing is written and nothing is modelled.
# =====================================================================

PROJ <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(".")
}, error = function(e) normalizePath("."))

source(file.path(PROJ, "config.R"))
source(file.path(PROJ, "R", "00_utils.R"))
source(file.path(PROJ, "R", "01_ingest.R"))

HYPNO_DIR <- resolve_path(HYPNO_DIR, PROJ)
MANIFEST  <- resolve_path(MANIFEST,  PROJ)

hr <- function(t) msg("\n---- %s %s", t, strrep("-", max(0, 62 - nchar(t))))

# ---------------------------------------------------------------------
hr("HYPNOGRAM DIRECTORY")
msg("HYPNO_DIR = %s", HYPNO_DIR)
if (!dir.exists(HYPNO_DIR)) {
  msg("  *** DOES NOT EXIST ***")
} else {
  all_files <- list.files(HYPNO_DIR, recursive = TRUE, full.names = TRUE)
  msg("  %d files total", length(all_files))
  ext <- tolower(tools::file_ext(all_files))
  print(sort(table(ext), decreasing = TRUE)[1:min(10, length(unique(ext)))])
  msg("  first 10 names:")
  for (f in utils::head(basename(all_files), 10)) msg("    %s", f)

  files <- tryCatch(list_hypnogram_files(HYPNO_DIR, HYPNO_GLOB, HYPNO_RECURSIVE),
                    error = function(e) { msg("  !! %s", conditionMessage(e)); character(0) })

  if (length(files)) {
    hr("FIRST FILE, RAW")
    msg("  %s", files[1])
    writeLines(paste0("    | ", utils::head(readLines(files[1], warn = FALSE), 12)))

    hr("FIRST FILE, AS PARSED")
    p <- tryCatch(read_one_hypnogram(files[1], EPOCH_SEC, HYPNO_STAGE_COL,
                                     NUMERIC_STAGE_CODING),
                  error = function(e) { msg("  !! %s", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      msg("  %d epochs -> %.1f min at %d s/epoch", nrow(p),
          nrow(p) * EPOCH_SEC / 60, EPOCH_SEC)
      msg("  raw stage values seen: %s",
          paste(utils::head(sort(unique(p$stage_raw)), 15), collapse = " | "))
      msg("  normalised stage counts:")
      print(table(p$stage))
      msg("  --> If 'UNSCORED' dominates, the wrong column was picked or the")
      msg("      numeric coding is wrong. Set HYPNO_STAGE_COL /")
      msg("      NUMERIC_STAGE_CODING in config.R.")
    }

    hr("PARSE CHECK ACROSS 25 FILES")
    samp <- files[unique(round(seq(1, length(files), length.out = min(25, length(files)))))]
    chk <- lapply(samp, function(f) {
      p <- tryCatch(read_one_hypnogram(f, EPOCH_SEC, HYPNO_STAGE_COL,
                                       NUMERIC_STAGE_CODING),
                    error = function(e) NULL)
      if (is.null(p)) return(data.frame(file = basename(f), epochs = NA,
                                        min = NA, pct_unscored = NA))
      data.frame(file = basename(f), epochs = nrow(p),
                 min = round(nrow(p) * EPOCH_SEC / 60, 1),
                 pct_unscored = round(100 * mean(p$stage == "UNSCORED"), 1))
    })
    print(do.call(rbind, chk), row.names = FALSE)
  }
}

# ---------------------------------------------------------------------
hr("SAMPLE MANIFEST")
msg("MANIFEST = %s", MANIFEST)
if (!file.exists(MANIFEST)) {
  msg("  *** DOES NOT EXIST ***")
} else {
  raw <- utils::read.csv(MANIFEST, stringsAsFactors = FALSE, check.names = FALSE)
  msg("  %d rows x %d columns", nrow(raw), ncol(raw))
  msg("  columns: %s", paste(names(raw), collapse = " | "))
  msg("\n  head:")
  print(utils::head(raw, 4))

  hr("MANIFEST, AS MAPPED")
  m <- tryCatch(read_manifest(MANIFEST),
                error = function(e) { msg("  !! %s", conditionMessage(e)); NULL })
  if (!is.null(m)) {
    print(utils::head(m, 6))
    msg("\n  subjects per group:")
    print(table(m$Genotype))
    msg("  sites: %s", paste(sort(unique(m$Site)), collapse = " | "))
    msg("  sex codes: %s", paste(sort(unique(m$Sex)), collapse = " | "))
    msg("  age range: %.1f - %.1f", min(m$Age, na.rm = TRUE), max(m$Age, na.rm = TRUE))
  }
}

# ---------------------------------------------------------------------
hr("RECORDING -> SUBJECT LINKAGE")
if (dir.exists(HYPNO_DIR) && file.exists(MANIFEST)) {
  files <- tryCatch(list_hypnogram_files(HYPNO_DIR, HYPNO_GLOB, HYPNO_RECURSIVE),
                    error = function(e) character(0))
  m <- tryCatch(read_manifest(MANIFEST, verbose = FALSE), error = function(e) NULL)
  if (length(files) && !is.null(m)) {
    rid <- tools::file_path_sans_ext(basename(files))
    link <- link_recordings_to_subjects(rid, m, SUBJECT_ID_MODE,
                                        SUBJECT_ID_REGEX, DEVICE_ID_REGEX)
    msg("  example links:")
    print(utils::head(link[!is.na(link$GUID), ], 8), row.names = FALSE)
    if (any(is.na(link$GUID))) {
      msg("  unlinked recording ids (first 10):")
      print(utils::head(link$recording_id[is.na(link$GUID)], 10))
      msg("  --> adjust SUBJECT_ID_MODE / SUBJECT_ID_REGEX in config.R")
    }
    n <- table(link$GUID[!is.na(link$GUID)])
    msg("  nights per subject: median %.0f, range %d-%d, %d subjects with >=1 night",
        stats::median(n), min(n), max(n), length(n))
  }
}

msg("\nInspection complete.")
