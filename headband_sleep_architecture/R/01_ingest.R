# =====================================================================
# 01_ingest.R -- read hypnograms and the sample manifest
#
# This is the ONLY format-dependent file in the pipeline.  Everything
# downstream consumes two tidy tables:
#
#   epochs   : recording_id | file | epoch_idx | stage      (one row / epoch)
#   manifest : GUID | Genotype | Age | Sex | Site           (one row / subject)
#
# If auto-detection picks the wrong column, set HYPNO_STAGE_COL in
# config.R, or edit read_one_hypnogram() below.
# =====================================================================

HYPNO_EXTS <- c("txt", "csv", "tsv", "dat", "stg", "hyp", "annot", "sco")

#' List candidate hypnogram files.
list_hypnogram_files <- function(dir, glob = NULL, recursive = TRUE) {
  if (!dir.exists(dir)) {
    stop("Hypnogram directory does not exist:\n  ", dir,
         "\nEdit HYPNO_DIR in config.R.", call. = FALSE)
  }
  if (!is.null(glob)) {
    f <- list.files(dir, pattern = utils::glob2rx(glob),
                    recursive = recursive, full.names = TRUE)
  } else {
    pat <- paste0("\\.(", paste(HYPNO_EXTS, collapse = "|"), ")$")
    f <- list.files(dir, pattern = pat, recursive = recursive,
                    full.names = TRUE, ignore.case = TRUE)
  }
  f <- f[!grepl("(^|/)(readme|README|\\.)", basename(f))]
  if (!length(f)) {
    stop("No hypnogram files found under:\n  ", dir,
         "\nSet HYPNO_GLOB in config.R to match your file naming.", call. = FALSE)
  }
  f
}

#' Read a single hypnogram file into a data.frame of epochs.
#'
#' Handles:
#'   * one-stage-per-line files (no delimiter)
#'   * delimited files with or without a header
#'   * an explicit epoch/time column, if present
read_one_hypnogram <- function(path, epoch_sec = 30,
                               stage_col = NULL, coding = "RK") {

  raw_lines <- readLines(path, warn = FALSE)

  # ---- Dreem export? ------------------------------------------------
  # A metadata preamble, then a tab-delimited table introduced by
  #   Sleep Stage <TAB> Time [hh:mm:ss] <TAB> Event <TAB> Duration[s]
  hdr <- grep("^Sleep Stage\t", raw_lines)[1]
  if (!is.na(hdr)) {
    body <- raw_lines[seq.int(hdr + 1L, length(raw_lines))]
    body <- body[nzchar(trimws(body))]
    if (!length(body)) return(NULL)
    parts <- strsplit(body, "\t", fixed = TRUE)

    fld <- function(i) vapply(parts, function(x)
      if (length(x) >= i) trimws(x[i]) else NA_character_, character(1))

    dur <- suppressWarnings(as.numeric(fld(4)))
    dur[is.na(dur) | dur <= 0] <- epoch_sec

    df <- data.frame(
      recording_id = tools::file_path_sans_ext(basename(path)),
      file         = path,
      epoch_idx    = seq_along(body),
      stage        = normalise_stage(fld(1), coding),
      stage_raw    = fld(1),
      dur_sec      = dur,
      stringsAsFactors = FALSE
    )
    return(df[, c("recording_id", "file", "epoch_idx", "stage", "stage_raw", "dur_sec")])
  }

  # ---- generic delimited / one-stage-per-line ------------------------
  lines <- raw_lines[nzchar(trimws(raw_lines))]
  # drop common comment/preamble lines
  lines <- lines[!grepl("^\\s*[#%;]", lines)]
  if (!length(lines)) return(NULL)

  delim <- sniff_delim(lines)

  if (is.null(delim)) {
    # ---- single column: every line is a stage ----
    hdr <- looks_like_header(lines[1], coding)
    vals <- if (hdr) lines[-1] else lines
    df <- data.frame(stage_raw = vals, stringsAsFactors = FALSE)

  } else {
    sep <- if (delim == " ") "" else delim   # "" = any whitespace in read.table
    first <- strsplit(lines[1], if (delim == " ") "[ ]+" else delim)[[1]]
    hdr <- looks_like_header(first, coding)

    df <- utils::read.table(
      text = paste(lines, collapse = "\n"),
      sep = sep, header = hdr, stringsAsFactors = FALSE,
      quote = "\"", comment.char = "", fill = TRUE,
      check.names = FALSE, colClasses = "character"
    )
    if (!hdr) names(df) <- paste0("V", seq_along(df))

    # ---- choose the stage column ----
    sc <- stage_col
    if (is.null(sc) || !(sc %in% names(df))) {
      byname <- pick_col(
        df,
        c("stage", "sleep_stage", "sleepstage", "hypnogram", "epoch_stage",
          "score", "scoring", "annotation", "event", "value"),
        what = "stage", required = FALSE
      )
      if (!is.na(byname) && stage_likeness(df[[byname]], coding) > 0.5) {
        sc <- byname
      } else {
        # pick whichever column looks most like sleep stages
        lik <- vapply(df, stage_likeness, numeric(1), coding = coding)
        sc <- names(df)[which.max(lik)]
        if (max(lik) < 0.5) {
          warning("No convincing stage column in ", basename(path),
                  " (best = '", sc, "', ", round(100 * max(lik)), "% stage-like)",
                  call. = FALSE)
        }
      }
    }
    df <- data.frame(stage_raw = df[[sc]], stringsAsFactors = FALSE)
  }

  df$stage     <- normalise_stage(df$stage_raw, coding)
  df$epoch_idx <- seq_len(nrow(df))
  df$file      <- path
  df$recording_id <- tools::file_path_sans_ext(basename(path))
  df$dur_sec   <- epoch_sec
  df[, c("recording_id", "file", "epoch_idx", "stage", "stage_raw", "dur_sec")]
}

#' Read every hypnogram under `dir`.
read_all_hypnograms <- function(dir, epoch_sec = 30, glob = NULL,
                                stage_col = NULL, coding = "RK",
                                recursive = TRUE, verbose = TRUE) {
  files <- list_hypnogram_files(dir, glob, recursive)
  if (verbose) msg("Reading %d hypnogram files from %s", length(files), dir)

  out <- vector("list", length(files))
  bad <- character(0)
  for (i in seq_along(files)) {
    out[[i]] <- tryCatch(
      read_one_hypnogram(files[i], epoch_sec, stage_col, coding),
      error = function(e) { bad <<- c(bad, files[i]); NULL }
    )
    if (verbose && i %% 100 == 0) msg("  ... %d / %d", i, length(files))
  }
  if (length(bad)) {
    warning(length(bad), " file(s) could not be parsed, e.g. ",
            basename(bad[1]), call. = FALSE)
  }
  ep <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(ep) || !nrow(ep)) stop("No epochs parsed.", call. = FALSE)

  # duplicate recording_ids would silently merge two nights -- make unique
  dup <- duplicated(unique(ep[, c("recording_id", "file")])$recording_id)
  if (any(dup)) {
    warning("Duplicate recording_id values across files; disambiguating by path.",
            call. = FALSE)
    ep$recording_id <- paste0(ep$recording_id, "__",
                              as.integer(factor(ep$file)))
  }
  if (verbose) {
    msg("  parsed %d recordings, %d epochs", length(unique(ep$recording_id)), nrow(ep))
    msg("  stage distribution: %s",
        paste(sprintf("%s=%.1f%%", names(table(ep$stage)),
                      100 * prop.table(table(ep$stage))), collapse = "  "))
  }
  ep
}

# ---------------------------------------------------------------------
# Sample manifest
# ---------------------------------------------------------------------

#' Read the sample manifest and standardise its columns.
read_manifest <- function(path, verbose = TRUE) {
  if (!file.exists(path)) {
    stop("Manifest file does not exist:\n  ", path,
         "\nEdit MANIFEST in config.R.", call. = FALSE)
  }
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (verbose) msg("Manifest: %d rows, %d columns", nrow(m), ncol(m))

  col_id   <- pick_col(m, c("GUID", "guid", "subject_id", "subjectid", "subject",
                            "participant_id", "sample_id", "IID", "id"), "subject ID")
  col_geno <- pick_col(m, c("final_group_class", "group_class", "Diagnosis",
                            "diagnosis", "genotype", "group", "cnv",
                            "carrier_status", "cnv_status", "locus_group",
                            "cohort_group"), "genotype/diagnosis")
  col_age  <- pick_col(m, c("age_final", "Age", "age", "age_years",
                            "age_at_visit", "age_at_scan", "visit_age"), "age")
  col_sex  <- pick_col(m, c("sex_final", "Sex", "sex", "gender",
                            "sex_at_birth"), "sex")
  col_site <- pick_col(m, c("redcap_data_access_group", "Site", "site",
                            "study_site", "recruitment_site", "cohort",
                            "center", "centre"), "site")
  col_dreem <- pick_col(m, c("dreem_id", "dev_beacon_dreem_id", "device_id"),
                        "dreem/device ID", required = FALSE)

  if (verbose) {
    msg("  mapped columns: id=%s  genotype=%s  age=%s  sex=%s  site=%s  device=%s",
        col_id, col_geno, col_age, col_sex, col_site, col_dreem %||% "<none>")
  }

  out <- data.frame(
    GUID     = trimws(as.character(m[[col_id]])),
    Genotype = normalise_genotype(m[[col_geno]]),
    GenotypeRaw = as.character(m[[col_geno]]),
    Age      = suppressWarnings(as.numeric(m[[col_age]])),
    Sex      = trimws(as.character(m[[col_sex]])),
    Site     = trimws(as.character(m[[col_site]])),
    stringsAsFactors = FALSE
  )
  out$DeviceID <- if (!is.na(col_dreem)) {
    # Dreem IDs may be stored as "SITE0000@beacon.study", sometimes padded
    sub("@.*$", "", trimws(as.character(m[[col_dreem]])))
  } else NA_character_
  out$DeviceID[is.na(out$DeviceID) | !nzchar(out$DeviceID)] <- NA_character_

  if (verbose) {
    unmapped <- sort(unique(out$GenotypeRaw[is.na(out$Genotype)]))
    if (length(unmapped)) {
      msg("  genotype labels NOT mapped to an analysis group (rows dropped): %s",
          paste(utils::head(unmapped, 20), collapse = " | "))
    }
    msg("  mapped group counts: %s",
        paste(sprintf("%s=%d", names(table(out$Genotype)), table(out$Genotype)),
              collapse = "  "))
  }

  # Record every device ID *before* genotype filtering, so the linkage step can
  # distinguish "device unknown to the manifest" from "device known but the
  # subject's genotype is not one of the five analysis groups".
  all_dev <- out[!is.na(out$DeviceID), c("DeviceID", "GUID", "GenotypeRaw", "Genotype")]

  out <- out[!is.na(out$Genotype), , drop = FALSE]

  # keep the row carrying a device ID when a subject appears more than once
  out <- out[order(is.na(out$DeviceID)), , drop = FALSE]

  # one row per subject; if the manifest is per-visit, take the first
  # non-missing value of each field
  out <- out[!duplicated(out$GUID), , drop = FALSE]
  out$Genotype <- factor(out$Genotype, levels = GENO_LEVELS)
  out$Locus    <- geno_locus(as.character(out$Genotype))
  out$CopyNum  <- geno_copy(as.character(out$Genotype))
  attr(out, "all_devices") <- all_dev
  out
}

# ---------------------------------------------------------------------
# Linking recordings to subjects
# ---------------------------------------------------------------------

#' Attach a subject ID (GUID) to each recording_id.
#'
#' mode = "device_id": the recording filename carries a device ID
#'   (Dreem: SFARI_Dreem_Prod_Hypno_<DEVICE>_<ISO datetime>__1), which is
#'   looked up in the manifest's device-ID column.  This is the mode for the
#'   SFARI hypnograms.
#' mode = "match_manifest": find the manifest GUID that occurs as a substring
#'   of the recording_id (longest match wins).  Robust when the file naming
#'   scheme is unknown.
#' mode = "regex": capture group 1 of SUBJECT_ID_REGEX.
#'
#' `manifest` must be the data.frame from read_manifest() for "device_id";
#' a plain vector of IDs suffices for the other two modes.
link_recordings_to_subjects <- function(recording_ids, manifest,
                                        mode = "match_manifest",
                                        regex = "^([^_.]+)",
                                        device_regex = "^SFARI_Dreem_Prod_Hypno_([A-Za-z]+[0-9]+)_",
                                        verbose = TRUE) {
  rid <- unique(recording_ids)
  manifest_ids <- if (is.data.frame(manifest)) manifest$GUID else manifest

  if (identical(mode, "device_id")) {
    if (!is.data.frame(manifest) || !"DeviceID" %in% names(manifest)) {
      stop("mode='device_id' needs the manifest data.frame with a DeviceID column.",
           call. = FALSE)
    }
    dev <- rep(NA_character_, length(rid))
    hit <- grepl(device_regex, rid)
    dev[hit] <- sub(paste0(device_regex, ".*$"), "\\1", rid[hit])
    dev <- toupper(trimws(dev))

    key <- toupper(trimws(manifest$DeviceID))
    ok  <- !is.na(key) & nzchar(key)
    guid <- manifest$GUID[ok][match(dev, key[ok])]

    if (verbose) {
      nd <- length(unique(dev[!is.na(dev)]))
      unmatched <- sort(unique(dev[!is.na(dev) & is.na(guid)]))
      msg("  device IDs parsed from filenames: %d unique", nd)
      if (length(unmatched)) {
        # split "absent from manifest" from "present but excluded genotype"
        allv <- attr(manifest, "all_devices")
        known <- if (!is.null(allv)) toupper(trimws(allv$DeviceID)) else character(0)
        excl  <- unmatched[unmatched %in% known]
        absent <- setdiff(unmatched, known)
        if (length(absent)) {
          msg("  device IDs absent from the manifest (%d): %s",
              length(absent), paste(absent, collapse = ", "))
        }
        if (length(excl)) {
          lab <- vapply(excl, function(d) {
            g <- allv$GenotypeRaw[known == d][1]
            sprintf("%s [%s]", d, if (is.na(g) || !nzchar(g)) "no genotype" else g)
          }, character(1))
          msg("  device IDs in the manifest but NOT in an analysis group (%d): %s",
              length(excl), paste(lab, collapse = ", "))
        }
      }
      if (any(is.na(dev))) {
        msg("  filenames not matching the device pattern: %d", sum(is.na(dev)))
      }
    }
    link <- data.frame(recording_id = rid, DeviceID = dev, GUID = guid,
                       stringsAsFactors = FALSE)
    if (verbose) {
      n_ok <- sum(!is.na(link$GUID))
      msg("Linked %d / %d recordings to a manifest subject (%.1f%%)",
          n_ok, nrow(link), 100 * n_ok / nrow(link))
    }
    return(link)
  }

  if (identical(mode, "regex")) {
    guid <- sub(regex, "\\1", rid)
    guid[guid == rid & !grepl(regex, rid)] <- NA_character_
  } else {
    # longest manifest ID that appears inside the recording id
    ord <- manifest_ids[order(nchar(manifest_ids), decreasing = TRUE)]
    ord <- ord[nchar(ord) >= 3]
    guid <- vapply(rid, function(r) {
      hit <- ord[vapply(ord, function(k) grepl(k, r, fixed = TRUE), logical(1))]
      if (length(hit)) hit[1] else NA_character_
    }, character(1), USE.NAMES = FALSE)
  }

  link <- data.frame(recording_id = rid, GUID = guid, stringsAsFactors = FALSE)
  if (verbose) {
    n_ok <- sum(!is.na(link$GUID))
    msg("Linked %d / %d recordings to a manifest subject (%.1f%%)",
        n_ok, nrow(link), 100 * n_ok / nrow(link))
    if (n_ok < nrow(link)) {
      msg("  unlinked examples: %s",
          paste(utils::head(link$recording_id[is.na(link$GUID)], 5), collapse = ", "))
    }
  }
  link
}
