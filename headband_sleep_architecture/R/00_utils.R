# =====================================================================
# 00_utils.R -- shared helpers: stage vocabulary, file sniffing, logging
# =====================================================================

msg <- function(...) cat(sprintf(...), "\n", sep = "")

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------
# Canonical sleep stages used everywhere downstream:
#   W, N1, N2, N3, R, MT, UNSCORED
# ---------------------------------------------------------------------
CANON_STAGES <- c("W", "N1", "N2", "N3", "R", "MT", "UNSCORED")

# Dreem exports (this cohort) use SLEEP-S0 / S1 / S2 / S3 / REM / MT, where
# S0 is wake.  Other conventions are included so the reader stays general.
.stage_alias <- list(
  W        = c("W", "WAKE", "AWAKE", "WK", "SLEEP STAGE W", "STAGE W",
               "SLEEP-S0", "SLEEP-W", "SLEEP-WAKE", "S0", "0"),
  N1       = c("N1", "S1", "STAGE 1", "SLEEP STAGE N1", "SLEEP STAGE 1",
               "SLEEP-S1", "SLEEP-N1", "NREM1", "1"),
  N2       = c("N2", "S2", "STAGE 2", "SLEEP STAGE N2", "SLEEP STAGE 2",
               "SLEEP-S2", "SLEEP-N2", "NREM2", "2"),
  N3       = c("N3", "S3", "S4", "N4", "STAGE 3", "STAGE 4", "SLEEP STAGE N3",
               "SLEEP STAGE 3", "SLEEP STAGE 4", "SLEEP-S3", "SLEEP-S4",
               "SLEEP-N3", "NREM3", "SWS", "DEEP", "3"),
  R        = c("R", "REM", "STAGE R", "SLEEP STAGE R", "SLEEP-REM", "SLEEP-R", "5"),
  MT       = c("MT", "MOVEMENT", "MOVEMENT TIME", "MOVE", "BODY MOVEMENT",
               "SLEEP-MT", "ARTIFACT", "ART", "6"),
  UNSCORED = c("UNSCORED", "UNKNOWN", "NS", "NOT SCORED", "?", "-1", "9", "999",
               "NA", "", "L", "LIGHTS OFF", "LIGHTS ON")
)

#' Normalise a vector of raw stage labels to CANON_STAGES.
#'
#' Numeric codes are ambiguous, so they are resolved via `coding`:
#'   "RK"   0=W 1=N1 2=N2 3=N3 4=N3 5=R 6=MT
#'   "AASM" 0=W 1=N1 2=N2 3=N3 4=R
normalise_stage <- function(x, coding = "RK") {
  raw <- toupper(trimws(as.character(x)))
  raw[is.na(raw)] <- ""

  out <- rep(NA_character_, length(raw))

  # Numeric codes first -- they take precedence over the alias table so the
  # `coding` convention is respected (e.g. "4" is N3 under RK, R under AASM).
  num_like <- grepl("^-?[0-9]+$", raw)
  if (any(num_like)) {
    n <- suppressWarnings(as.integer(raw[num_like]))
    map <- if (identical(coding, "AASM")) {
      c("0" = "W", "1" = "N1", "2" = "N2", "3" = "N3", "4" = "R",
        "5" = "R", "6" = "MT", "-1" = "UNSCORED", "9" = "UNSCORED")
    } else {
      c("0" = "W", "1" = "N1", "2" = "N2", "3" = "N3", "4" = "N3",
        "5" = "R", "6" = "MT", "7" = "MT", "-1" = "UNSCORED", "9" = "UNSCORED")
    }
    out[num_like] <- unname(map[as.character(n)])
  }

  # Then string aliases for anything still unresolved.
  todo <- is.na(out)
  if (any(todo)) {
    for (st in names(.stage_alias)) {
      hit <- todo & raw %in% .stage_alias[[st]]
      out[hit] <- st
      todo <- todo & !hit
    }
  }

  out[is.na(out)] <- "UNSCORED"
  factor(out, levels = CANON_STAGES)
}

#' Fraction of a character vector that looks like a scorable sleep stage.
#' Used to pick the stage column when the file has no informative header.
stage_likeness <- function(x, coding = "RK") {
  s <- normalise_stage(x, coding)
  mean(s != "UNSCORED", na.rm = TRUE)
}

# ---------------------------------------------------------------------
# Delimited-file sniffing
# ---------------------------------------------------------------------

#' Guess the field delimiter of a text file from its first data lines.
sniff_delim <- function(lines) {
  cand <- c("\t" = "\t", "," = ",", ";" = ";", "|" = "|", " " = "[ ]+")
  probe <- utils::head(lines[nzchar(trimws(lines))], 20)
  if (!length(probe)) return(NULL)
  counts <- vapply(cand, function(d) {
    n <- vapply(strsplit(probe, d), length, integer(1))
    # a good delimiter yields >1 field, consistently
    if (stats::median(n) > 1 && length(unique(n)) <= 2) stats::median(n) else 0
  }, numeric(1))
  if (all(counts == 0)) return(NULL)          # single-column file
  names(cand)[which.max(counts)]
}

#' Does the first row look like a header rather than data?
looks_like_header <- function(first_fields, coding = "RK") {
  f <- trimws(first_fields)
  if (!length(f)) return(FALSE)
  numeric_frac <- mean(grepl("^-?[0-9.]+$", f))
  stagey <- stage_likeness(f, coding)
  # A header row is mostly non-numeric AND mostly not stage tokens.
  numeric_frac < 0.5 && stagey < 0.5
}

# ---------------------------------------------------------------------
# Column-name matching for the sample manifest
# ---------------------------------------------------------------------

#' Find the first column in `df` whose (normalised) name matches any candidate.
pick_col <- function(df, candidates, what = "column", required = TRUE) {
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  nm <- norm(names(df))
  # exact, case-sensitive match wins -- this manifest has BOTH `guid` and
  # `GUID` (different ID formats), so case must not be collapsed first
  for (cand in candidates) {
    i <- match(cand, names(df))
    if (!is.na(i)) return(names(df)[i])
  }
  for (cand in candidates) {
    i <- match(norm(cand), nm)
    if (!is.na(i)) return(names(df)[i])
  }
  # fall back to substring match
  for (cand in candidates) {
    i <- grep(norm(cand), nm, fixed = TRUE)
    if (length(i)) return(names(df)[i[1]])
  }
  if (required) {
    stop(sprintf(
      "Could not find a %s column. Looked for: %s.\nColumns present: %s",
      what, paste(candidates, collapse = ", "), paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  NA_character_
}

# ---------------------------------------------------------------------
# Genotype-label normalisation
# ---------------------------------------------------------------------

#' Map free-text genotype/diagnosis labels onto the five analysis groups.
#'
#' Returns NA for anything that is not exactly one of the five; those rows are
#' dropped.  In particular:
#'   * carriers of TWO CNVs ("22q11.2_Deletion & 16p11.2_Deletion") are NA --
#'     they belong to neither locus's dose-response series;
#'   * triplications are NA -- copy number 4 is outside the 1/2/3 contrast;
#'   * "no_cnv" is Control.
normalise_genotype <- function(x) {
  raw <- as.character(x)
  s   <- toupper(gsub("[^A-Za-z0-9]", "", raw))

  is16  <- grepl("16P", s)
  is22  <- grepl("22Q", s)
  isdel <- grepl("DEL", s)
  isdup <- grepl("DUP", s)
  isctl <- grepl("^NOCNV$|^NONE$|CONTROL|CTRL|NONCARRIER|NEGATIVE|^FS$|FAMILYCONTROL|SIBLING|TYPICAL|TD$|HC$", s)

  # Reject anything naming two loci, two dosage directions, or a triplication.
  ambiguous <- (is16 & is22) | (isdel & isdup) | grepl("&|\\+|AND", toupper(raw)) |
               grepl("TRIPLIC|TETRASOM", s)

  out <- rep(NA_character_, length(s))
  out[!ambiguous & is16 & isdel] <- "16pDel"
  out[!ambiguous & is16 & isdup] <- "16pDup"
  out[!ambiguous & is22 & isdel] <- "22qDel"
  out[!ambiguous & is22 & isdup] <- "22qDup"
  out[is.na(out) & !ambiguous & isctl] <- "Control"
  out
}

GENO_LEVELS <- c("16pDel", "16pDup", "22qDel", "22qDup", "Control")

#' Locus ("16p" / "22q") and copy number (1/2/3) implied by a genotype label.
geno_locus <- function(g) {
  ifelse(grepl("^16p", g), "16p", ifelse(grepl("^22q", g), "22q", "Control"))
}
geno_copy <- function(g) {
  ifelse(grepl("Del$", g), 1L, ifelse(grepl("Dup$", g), 3L, 2L))
}

# ---------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
           formatC(p, format = "f", digits = 4)))
}

#' Resolve a path from config.R: used as-is if absolute (drive letter, leading
#' slash, or UNC), otherwise taken relative to the project directory.
resolve_path <- function(p, base) {
  if (grepl("^([A-Za-z]:|/|\\\\)", p)) p else file.path(base, p)
}

ensure_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
