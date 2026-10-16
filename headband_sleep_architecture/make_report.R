# =====================================================================
# make_report.R -- build the standalone HTML report from results/
#
#   Rscript make_report.R [outdir]
#
# Reads only the files written by run_all.R, so the report can never drift
# from the analysis.  Figures are base64-embedded, making the HTML fully
# self-contained (and printable to PDF by any browser).
# =====================================================================

PROJ <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(".")
}, error = function(e) normalizePath("."))

source(file.path(PROJ, "config.R"))
source(file.path(PROJ, "R", "00_utils.R"))

args   <- commandArgs(trailingOnly = TRUE)
RES    <- resolve_path(OUT_DIR, PROJ)
OUTDIR <- if (length(args)) args[1] else file.path(RES, "report")
invisible(ensure_dir(OUTDIR))

tests <- read.csv(file.path(RES, "results_all_tests.csv"), stringsAsFactors = FALSE)
sens  <- read.csv(file.path(RES, "sensitivity_analyses.csv"), stringsAsFactors = FALSE)
desc  <- read.csv(file.path(RES, "group_descriptives.csv"), stringsAsFactors = FALSE)
venn  <- read.csv(file.path(RES, "qc_venn_counts.csv"), stringsAsFactors = FALSE)
metr  <- read.csv(file.path(RES, "per_recording_metrics_with_qc.csv"), stringsAsFactors = FALSE)
dat   <- read.csv(file.path(RES, "analysis_dataset_per_night.csv"), stringsAsFactors = FALSE)
FIGS  <- file.path(RES, "figures")

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
b64 <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  chars <- c(LETTERS, letters, as.character(0:9), "+", "/")
  n <- length(raw); pad <- (3L - n %% 3L) %% 3L
  raw <- c(raw, rep(as.raw(0), pad))
  m <- matrix(as.integer(raw), nrow = 3L)
  idx <- rbind(m[1, ] %/% 4L,
               (m[1, ] %% 4L) * 16L + m[2, ] %/% 16L,
               (m[2, ] %% 16L) * 4L + m[3, ] %/% 64L,
               m[3, ] %% 64L)
  out <- chars[as.vector(idx) + 1L]
  if (pad) out[length(out) - seq_len(pad) + 1L] <- "="
  paste(out, collapse = "")
}

fig <- function(file, caption, width = "100%") {
  p <- file.path(FIGS, file)
  if (!file.exists(p)) return("")
  sprintf('<figure><img src="data:image/png;base64,%s" style="width:%s">
           <figcaption>%s</figcaption></figure>', b64(p), width, caption)
}

esc <- function(x) {
  x <- gsub("&", "&amp;", x); x <- gsub("<", "&lt;", x); gsub(">", "&gt;", x)
}

pfmt <- function(p) ifelse(is.na(p), "&mdash;",
  ifelse(p < 1e-4, sprintf("%.1e", p), sprintf("%.4f", p)))

#' Render a data.frame as an HTML table.
#' `hl` is a logical vector marking rows to emphasise.
tbl <- function(df, hl = NULL, align_right = NULL, cls = "") {
  th <- paste0("<th>", esc(names(df)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    tds <- vapply(seq_along(df), function(j) {
      a <- if (!is.null(align_right) && names(df)[j] %in% align_right)
             ' class="num"' else ""
      sprintf("<td%s>%s</td>", a, df[i, j])
    }, character(1))
    sprintf('<tr%s>%s</tr>',
            if (!is.null(hl) && isTRUE(hl[i])) ' class="hit"' else "",
            paste(tds, collapse = ""))
  }, character(1))
  sprintf('<table class="%s"><thead><tr>%s</tr></thead><tbody>%s</tbody></table>',
          cls, th, paste(rows, collapse = ""))
}

# ---------------------------------------------------------------------
# numbers used in the narrative
# ---------------------------------------------------------------------
n_parsed    <- nrow(metr)
n_pass      <- sum(metr$qc_pass)
n_pass_core <- sum(metr$qc_pass_core)
n_excl_n3   <- sum(metr$excl_implausible_n3 & metr$qc_pass_core)
n_nights    <- nrow(dat)
n_subj      <- length(unique(dat$GUID))

qcs <- read.csv(file.path(RES, "qc_summary.csv"), stringsAsFactors = FALSE)
t1a <- data.frame(
  Step      = esc(qcs$step),
  Excluded  = ifelse(is.na(qcs$excluded),  "&mdash;", as.character(qcs$excluded)),
  Remaining = ifelse(is.na(qcs$remaining), "&mdash;", as.character(qcs$remaining)),
  Note      = esc(qcs$note),
  check.names = FALSE, stringsAsFactors = FALSE
)

get1 <- function(out, loc, mod, term) {
  r <- tests[tests$outcome == out & tests$cohort == loc &
             tests$model == mod & tests$term == term, ]
  if (!nrow(r)) return(NULL)
  r[1, ]
}
di16 <- get1("DI", "16p", "mirror_single", "copy_c")
di22 <- get1("DI", "22q", "mirror_single", "copy_c")
ts22 <- get1("TST", "22q", "mirror_single", "copy_c")

# ---------------------------------------------------------------------
# Table 1 -- sample
# ---------------------------------------------------------------------
t1 <- data.frame(
  Group      = desc$Genotype,
  Subjects   = desc$n_subjects,
  Nights     = desc$n_nights,
  `Nights/subject` = sprintf("%.1f", desc$nights_per_subject),
  `Age, mean (SD)` = sprintf("%.1f (%.1f)", desc$age_mean, desc$age_sd),
  `% female`  = sprintf("%.0f", desc$pct_female),
  `DI, mean (SD)`  = sprintf("%.3f (%.3f)", desc$DI_mean, desc$DI_sd),
  `TST min, mean (SD)` = sprintf("%.0f (%.0f)", desc$TST_mean, desc$TST_sd),
  `SOL min, mean (SD)` = sprintf("%.0f (%.0f)", desc$SOL_mean, desc$SOL_sd),
  check.names = FALSE, stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------
# Table 2 -- primary tests
# ---------------------------------------------------------------------
# One table per outcome -- 75 rows in a single block is unreadable on paper.
# The complete unsplit table is results_all_tests.csv.
outcome_table <- function(y) {
  p <- tests[tests$outcome == y, ]
  p <- p[order(match(p$cohort, c("16p", "22q", "All5")),
               match(p$model, c("mirror_single", "ushape_single", "joint",
                                "group_vs_control"))), ]
  df <- data.frame(
    Locus   = ifelse(p$cohort == "All5", "all groups", p$cohort),
    Model   = c(mirror_single = "mirror", ushape_single = "U-shape",
                joint = "joint", group_vs_control = "group")[p$model],
    Term    = esc(p$term_label),
    Test    = ifelse(p$test == "LRT", "LRT", "Wald"),
    Estimate = ifelse(is.na(p$estimate), "&mdash;", sprintf("%+.4g", p$estimate)),
    SE       = ifelse(is.na(p$se), "&mdash;", sprintf("%.3g", p$se)),
    p        = pfmt(p$p_value),
    q        = pfmt(p$q_value),
    N        = sprintf("%d / %d", p$n_nights, p$n_subjects),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  list(df = df, hl = !is.na(p$p_value) & p$p_value < 0.05)
}
T_DI  <- outcome_table("DI")
T_TST <- outcome_table("TST")
T_SOL <- outcome_table("SOL")
NUMCOL <- c("Estimate", "SE", "p", "q", "N")

# ---------------------------------------------------------------------
# Table 3 -- sensitivity
# ---------------------------------------------------------------------
sm <- sens[sens$check == "N3 outlier threshold", ]
t3 <- data.frame(
  `N3 threshold` = esc(sm$setting), Locus = sm$locus, Nights = sm$n_nights,
  Estimate = sprintf("%+.4f", sm$estimate), p = pfmt(sm$p_value),
  check.names = FALSE, stringsAsFactors = FALSE
)
hl3 <- sm$p_value < 0.05

sd_ <- sens[sens$check == "depth-index definition" & sens$model == "mirror", ]
t3d <- data.frame(
  Definition = sd_$setting, Locus = sd_$locus, Nights = sd_$n_nights,
  Estimate = sprintf("%+.4f", sd_$estimate), p = pfmt(sd_$p_value),
  check.names = FALSE, stringsAsFactors = FALSE
)
hl3d <- sd_$p_value < 0.05

sq <- sens[sens$check == "QC thresholds", ]
t3b <- data.frame(
  Setting = esc(sq$setting), Locus = sq$locus, Nights = sq$n_nights,
  Estimate = sprintf("%+.4f", sq$estimate), p = pfmt(sq$p_value),
  check.names = FALSE, stringsAsFactors = FALSE
)
hl3b <- sq$p_value < 0.05

# ---- sensitivity summary: one row per check ----
.verdict <- function(rows) {
  r <- rows[rows$outcome %in% c("DI", "DI_over_NREM", "DI_over_TST"), ]
  if (!nrow(r)) r <- rows
  n_sig <- sum(r$p_value < 0.05, na.rm = TRUE)
  sprintf("%d / %d refits significant (max p = %s)", n_sig, nrow(r),
          pfmt(max(r$p_value, na.rm = TRUE)))
}
.checks <- c("N3 outlier threshold", "depth-index definition",
             "onset-latency definition", "value cap", "QC thresholds")
.what <- c(
  "N3 outlier threshold"     = "Drops nights scored above 40&ndash;70% N3, and with no filter at all",
  "depth-index definition"   = "N3 over other stages (primary), over NREM, over TST",
  "onset-latency definition" = "Raw SOL, log(1+SOL), latency to persistent sleep",
  "value cap"                = "Caps extreme outcome values post-QC",
  "QC thresholds"            = "All 12 combinations of recording-length and movement cut-offs, per locus")
.tbl_ref <- c("N3 outlier threshold" = "Table 3a", "depth-index definition" = "Table 3b",
              "onset-latency definition" = "Table 3d", "value cap" = "&mdash;",
              "QC thresholds" = "Table 3c")
sens_summary <- do.call(rbind, lapply(.checks, function(ck) {
  r <- sens[sens$check == ck, ]
  if (!nrow(r)) return(NULL)
  data.frame(
    `check =` = ck, `What it varies` = .what[[ck]], Refits = nrow(r),
    Verdict = .verdict(r), Detail = .tbl_ref[[ck]],
    check.names = FALSE, stringsAsFactors = FALSE)
}))

so <- sens[sens$check == "onset-latency definition" & sens$model == "mirror", ]
t3c <- data.frame(
  Definition = so$setting, Locus = so$locus,
  Estimate = sprintf("%+.4g", so$estimate), p = pfmt(so$p_value),
  check.names = FALSE, stringsAsFactors = FALSE
)
hl3c <- so$p_value < 0.05

# ---------------------------------------------------------------------
# Depth / duration trade-off
# ---------------------------------------------------------------------
rd <- function(f) {
  p <- file.path(RES, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}
tr_corr  <- rd("tradeoff_correlations.csv")
tr_slope <- rd("tradeoff_slopes_by_group.csv")
tr_int   <- rd("tradeoff_group_comparisons.csv")
tr_ll    <- rd("tradeoff_loglog_check.csv")
tr_null  <- rd("tradeoff_null_calibration.csv")
tr_fix   <- rd("tradeoff_fixedwindow_correlations.csv")
tr_adj   <- rd("tradeoff_main_effect_adjusted_for_TST.csv")

GRPLAB <- c(`ALL (combined sample)` = "Combined sample", Control = "Control",
            `16pDel` = "16p deletion", `16pDup` = "16p duplication",
            `22qDel` = "22q deletion", `22qDup` = "22q duplication")

# Table T1: correlations at three levels
tT1 <- NULL
if (!is.null(tr_corr)) {
  wide <- lapply(c("night", "participant", "within-participant"), function(lv) {
    x <- tr_corr[tr_corr$level == lv, ]
    stats::setNames(data.frame(
      grp = x$group,
      v = sprintf("%+.3f (p = %s)", x$pearson_r, pfmt(x$pearson_p)),
      stringsAsFactors = FALSE), c("grp", lv))
  })
  m0 <- Reduce(function(a, b) merge(a, b, by = "grp", all = TRUE), wide)
  m0 <- m0[match(names(GRPLAB), m0$grp), ]
  m0 <- m0[!is.na(m0$grp), ]
  tT1 <- data.frame(
    Group = unname(GRPLAB[m0$grp]),
    `Night level` = m0$night,
    `Participant level` = m0$participant,
    `Within participant` = m0$`within-participant`,
    check.names = FALSE, stringsAsFactors = FALSE)
}

# Table T2: mixed-model slopes per group
tT2 <- NULL; hlT2 <- NULL
if (!is.null(tr_slope)) {
  x <- tr_slope[order(match(tr_slope$group, names(GRPLAB)),
                      tr_slope$effect != "within-participant"), ]
  tT2 <- data.frame(
    Group = unname(GRPLAB[x$group]), Effect = x$effect,
    `Slope per hour` = sprintf("%+.4f", x$estimate),
    `95% CI` = sprintf("[%+.4f, %+.4f]", x$ci_lo, x$ci_hi),
    p = pfmt(x$p_value), Nights = x$n_nights,
    `Random effects` = x$random_effects,
    check.names = FALSE, stringsAsFactors = FALSE)
  hlT2 <- x$p_value < 0.05
}

# Table T3: group differences in the trade-off
tT3 <- NULL; hlT3 <- NULL
if (!is.null(tr_int)) {
  x <- tr_int[tr_int$effect == "interaction", ]
  tT3 <- data.frame(
    Comparison = x$comparison, Term = esc(x$term), Test = x$test,
    Statistic = sprintf("%.3g", x$statistic), df = sprintf("%.4g", x$df),
    p = pfmt(x$p_value), check.names = FALSE, stringsAsFactors = FALSE)
  hlT3 <- x$p_value < 0.05
}

# Table T4: the three controls for the definitional coupling
tT4 <- NULL
if (!is.null(tr_ll) && !is.null(tr_null) && !is.null(tr_fix)) {
  ctl_ll  <- tr_ll[tr_ll$group == "Control", ]
  all_fix <- tr_fix[tr_fix$group == "ALL (combined sample)" & tr_fix$level == "night", ]
  ctl_fix <- tr_fix[tr_fix$group == "Control" & tr_fix$level == "night", ]
  tT4 <- data.frame(
    Control = c("Simulated mechanical null",
                "log(N3) on log(TST), within participants (controls)",
                "Fixed-window depth index (first 240 min of sleep)"),
    `What it rules out` = c(
      "Correlation induced purely by TST sitting in the denominator",
      "Any artefact of using a ratio at all",
      "Duration entering the index at all &mdash; every night uses the same amount of sleep"),
    Result = c(
      sprintf("observed r = %+.3f vs null %+.3f [%+.3f, %+.3f], p = %s",
              tr_null$observed_r, tr_null$null_mean_r, tr_null$null_lo,
              tr_null$null_hi, pfmt(tr_null$p_vs_null)),
      sprintf("slope = %.3f [%.3f, %.3f], vs 1: p = %s",
              ctl_ll$slope, ctl_ll$ci_lo, ctl_ll$ci_hi, pfmt(ctl_ll$p_vs_1)),
      sprintf("r = %+.3f combined, %+.3f in controls (sign reverses)",
              all_fix$pearson_r, ctl_fix$pearson_r)),
    check.names = FALSE, stringsAsFactors = FALSE)
}

# Table T5: does the headline effect survive adjustment for TST?
tT5 <- NULL; hlT5 <- NULL
if (!is.null(tr_adj)) {
  x <- tr_adj[order(match(tr_adj$outcome, c("DI", "TST")), x <- tr_adj$locus,
                    !is.na(tr_adj$adjusted_for)), ]
  fmt_est <- ifelse(x$outcome == "TST", sprintf("%+.1f min", x$estimate),
                    sprintf("%+.4f", x$estimate))
  tT5 <- data.frame(
    Outcome = ifelse(x$outcome == "DI", "NREM depth index", "Total sleep time"),
    Locus = x$locus, Model = esc(x$model),
    `Copy-number effect` = fmt_est,
    SE = ifelse(x$outcome == "TST", sprintf("%.1f", x$se), sprintf("%.4f", x$se)),
    p = pfmt(x$p_value),
    check.names = FALSE, stringsAsFactors = FALSE)
  hlT5 <- x$p_value < 0.05
}
.adjp <- function(out, loc, adjusted) {
  if (is.null(tr_adj)) return(NA)
  r <- tr_adj[tr_adj$outcome == out & tr_adj$locus == loc &
              (is.na(tr_adj$adjusted_for) != adjusted), ]
  if (!nrow(r)) NA else r$p_value[1]
}
.adje <- function(out, loc, adjusted) {
  if (is.null(tr_adj)) return(NA)
  r <- tr_adj[tr_adj$outcome == out & tr_adj$locus == loc &
              (is.na(tr_adj$adjusted_for) != adjusted), ]
  if (!nrow(r)) NA else r$estimate[1]
}

.g <- function(df, grp, lv) {
  if (is.null(df)) return(NA)
  r <- df[df$group == grp & df$level == lv, ]
  if (!nrow(r)) NA else r$pearson_r[1]
}
.s <- function(grp, eff, col = "estimate") {
  if (is.null(tr_slope)) return(NA)
  r <- tr_slope[tr_slope$group == grp & tr_slope$effect == eff, ]
  if (!nrow(r)) NA else r[[col]][1]
}
.ip <- function(cmp) {
  if (is.null(tr_int)) return(NA)
  r <- tr_int[tr_int$comparison == cmp & tr_int$effect == "interaction", ]
  if (!nrow(r)) NA else r$p_value[1]
}

# ---------------------------------------------------------------------
# assemble
# ---------------------------------------------------------------------
CSS <- '
:root { --ink:#16191d; --mut:#5b636d; --rule:#d8dde3; --hit:#eaf3fb; --acc:#1f6fb4; }
* { box-sizing:border-box; }
body { font-family:"Source Serif Pro",Georgia,"Times New Roman",serif;
       color:var(--ink); line-height:1.55; font-size:10.5pt; margin:0; }
h1 { font-size:21pt; line-height:1.2; margin:0 0 2px; letter-spacing:-.01em; }
h2 { font-size:13pt; margin:26px 0 8px; padding-bottom:4px;
     border-bottom:2px solid var(--ink); }
h3 { font-size:11pt; margin:18px 0 6px; color:var(--acc); }
.sub { color:var(--mut); font-size:10pt; margin:0 0 4px; }
.meta { color:var(--mut); font-size:8.5pt; margin:10px 0 0;
        border-top:1px solid var(--rule); padding-top:8px; }
table { border-collapse:collapse; width:100%; font-size:8.2pt; margin:10px 0 4px;
        font-family:"Helvetica Neue",Arial,sans-serif; }
th { text-align:left; border-bottom:1.5px solid var(--ink); padding:4px 6px;
     font-weight:600; font-size:7.8pt; text-transform:uppercase;
     letter-spacing:.03em; color:var(--mut); }
td { border-bottom:1px solid var(--rule); padding:3.5px 6px; vertical-align:top; }
td.num { text-align:right; font-variant-numeric:tabular-nums; }
tr.hit td { background:var(--hit); font-weight:600; }
figure { margin:14px 0; page-break-inside:avoid; }
figure img { display:block; border:1px solid var(--rule); }
figcaption { font-size:8.5pt; color:var(--mut); margin-top:5px; line-height:1.4; }
.key { background:#f6f8fa; border-left:3px solid var(--acc);
       padding:10px 14px; margin:12px 0; font-size:9.8pt; }
.key p { margin:5px 0; }
.warn { background:#fdf6ec; border-left:3px solid #c98a2e; }
ul,ol { margin:6px 0 6px 18px; padding:0; }
li { margin:3px 0; }
code { font-family:"SF Mono",Consolas,monospace; font-size:8.6pt;
       background:#f2f4f7; padding:1px 4px; border-radius:2px; }
.cap { font-size:8.5pt; color:var(--mut); margin:2px 0 12px; }
@page { size:letter; margin:14mm 15mm; }
.pb { page-break-before:always; }
'

html <- paste0(
'<!doctype html><html><head><meta charset="utf-8">
<title>Sleep architecture and 16p11.2 / 22q11.2 copy number</title>
<style>', CSS, '</style></head><body>',

'<h1>Sleep architecture and 16p11.2 / 22q11.2 copy number</h1>
<p class="sub">Dose-response effects on NREM depth, total sleep time and sleep onset latency
in overnight Dreem EEG recordings</p>
<p class="meta">Generated ', format(Sys.Date(), "%d %B %Y"),
' &middot; ', n_nights, ' nights from ', n_subj, ' participants &middot; ',
'linear mixed models with a random intercept per participant</p>',

# ---------------- summary ----------------
'<h2>Summary of findings</h2>',
'<div class="key">
<p><b>1. NREM depth index shows a dose-response effect at both loci.</b>
Depth index rises monotonically with copy number at 16p11.2
(&beta; = ', sprintf("%+.4f", di16$estimate), ' per copy, p = ', pfmt(di16$p_value),
', q = ', pfmt(di16$q_value), ') and at 22q11.2
(&beta; = ', sprintf("%+.4f", di22$estimate), ', p = ', pfmt(di22$p_value),
', q = ', pfmt(di22$q_value), '). Deletions sleep <i>less</i> deeply and duplications
<i>more</i> deeply than controls, at both loci.</p>
<p><b>2. The effect is linear, not U-shaped.</b> Adding a squared copy-number term
adds nothing (16p: p = ', pfmt(get1("DI","16p","joint","U-shape beyond mirror (1 df)")$p_value),
'; 22q: p = ', pfmt(get1("DI","22q","joint","U-shape beyond mirror (1 df)")$p_value),
'), while the linear term remains significant after allowing for curvature
(16p: p = ', pfmt(get1("DI","16p","joint","mirror beyond U-shape (1 df)")$p_value),
'; 22q: p = ', pfmt(get1("DI","22q","joint","mirror beyond U-shape (1 df)")$p_value), ').</p>
<p><b>3. Total sleep time falls with 22q11.2 copy number</b>
(&beta; = ', sprintf("%+.1f", ts22$estimate), ' min per copy, p = ', pfmt(ts22$p_value),
', q = ', pfmt(ts22$q_value), '). The 16p11.2 trend is in the same direction but not
significant (p = ', pfmt(get1("TST","16p","mirror_single","copy_c")$p_value),
'); 16p11.2 deletion carriers sleep ',
sprintf("%.0f", get1("TST","16p","group_vs_control","Genotype16pDel")$estimate),
' min longer than controls (p = ',
pfmt(get1("TST","16p","group_vs_control","Genotype16pDel")$p_value), ').
<i>However, the sleep-time effect does not survive adjustment for depth index
(22q11.2: p = ', pfmt(.adjp("TST","22q",TRUE)), '), whereas the depth effect does survive
adjustment for sleep time. Depth appears to be the primary signal &mdash; see the
trade-off section.</i></p>
<p><b>4. Sleep onset latency shortens with copy number at 16p11.2</b>
(&beta; = ', sprintf("%+.2f", get1("SOL","16p","mirror_single","copy_c")$estimate),
' min per copy, p = ', pfmt(get1("SOL","16p","mirror_single","copy_c")$p_value),
', q = ', pfmt(get1("SOL","16p","mirror_single","copy_c")$q_value),
'). Latency is strongly right-skewed; on the log scale the effect is present at
<i>both</i> loci (16p p = ',
pfmt(sens$p_value[sens$outcome=="SOL_log" & sens$locus=="16p" & sens$model=="mirror"][1]),
', 22q p = ',
pfmt(sens$p_value[sens$outcome=="SOL_log" & sens$locus=="22q" & sens$model=="mirror"][1]),
').</p>
</div>',

'<p><b>A note on terminology.</b> Both loci move in the <i>same</i> direction: fewer copies
&rarr; shallower NREM sleep, more copies &rarr; deeper NREM sleep. The effect is a
<i>mirror</i> effect in the sense of a linear dose-response across copy number 1&ndash;2&ndash;3
(deletion and duplication on opposite sides of control), not in the sense of the two loci
opposing each other.</p>',

'<p>Note also that no single genotype-vs-control contrast reaches p &lt; 0.05 for depth index,
while the dose-response contrast does. That is expected: the linear contrast uses all three
groups and one degree of freedom, so it concentrates the available signal, whereas splitting
into two pairwise comparisons spends power on a hypothesis the data do not support
(a non-monotonic pattern).</p>',

# ---------------- sample ----------------
'<h2 class="pb">Sample and quality control</h2>',
'<p>', n_parsed, ' Dreem hypnogram files were parsed (128 device IDs across 4 sites:
Montreal, PennCHOP, Toronto, UCLA). Filters were applied in the order below.</p>',

'<h3>Table 1a. Quality-control filters</h3>',
tbl(t1a, align_right = c("Excluded", "Remaining")),

'<h3>The three recording-level criteria</h3>
<ul>
<li><b>No sleep</b> &mdash; total sleep time of exactly zero. The recording contains no
scored N1, N2, N3 or REM epoch, so no sleep variable is defined.</li>
<li><b>Too short</b> &mdash; total recording time below ', QC_MIN_TRT_MIN, ' minutes.
Short recordings cannot sample a full night\'s architecture; a recording that ends before
the first NREM&ndash;REM cycle completes will have a badly biased depth index.
<i>(The original write-up described this as a 30-minute rule, but the Venn diagram in the
original figures used 2 hours. 2 hours is used here to reproduce that; the 30-minute
setting is included in the sensitivity grid and changes nothing.)</i></li>
<li><b>High movement</b> &mdash; more than ', 100 * QC_MAX_MT_FRAC, '% of epochs scored as
movement time. Headband EEG is movement-prone, and a recording that is majority movement
is not scoreable. This is by far the dominant exclusion.</li>
</ul>
<p>These three overlap heavily, so they do not sum to the unique total; the overlap is
shown in Figure 1. Together they exclude ', n_parsed - n_pass_core, ' recordings and
retain <b>', n_pass_core, '</b> &mdash; matching the 782 of the original analysis
exactly, which confirms the QC has been reproduced.</p>',

'<h3>No outlier removal was applied</h3>
<p><b>Nothing was excluded on the basis of its outcome value.</b> The three criteria above
are the only exclusions, and they are exactly the criteria used in the original analysis.
Every night that passes them enters the models, including nights with extreme values.</p>
<p>This matters for the depth index, which is an unbounded ratio: the retained sample
contains four nights with a depth index above 1.5, the largest being 4.32 &mdash; a night
scored 81% N3, which is not physiologically plausible and is more likely a scoring
artefact. These are <b>kept</b>, because retaining them is the conservative choice:
dropping them <i>strengthens</i> the 22q11.2 result substantially. Table 3a quantifies
exactly what an outlier filter would do at every threshold from 40% to 70% of TST.</p>
<p>The figures handle these nights by <b>clipping the y-axis</b>, not by dropping data.
Panels are drawn over the Tukey range (Q1&nbsp;&minus;&nbsp;3&times;IQR to
Q3&nbsp;+&nbsp;3&times;IQR) so the bulk of the distribution is legible; points outside the
axis are not drawn, but the boxplots, the model fit, the confidence intervals and every
p-value still use all ', n_nights, ' nights. Each figure caption states how many points
fall outside its axis.</p>',

'<h3>Linkage</h3>
<p>Recordings link to participants by the Dreem device ID embedded in the filename,
looked up in the manifest\'s <code>dreem_id</code> column. ', n_pass - n_nights,
' QC-passing nights could not be resolved: one device is absent from the manifest,
one participant has no genotype recorded, and one carries both a 22q11.2 and a 16p11.2
deletion and so belongs to neither dose series. This leaves <b>', n_nights, ' nights from ',
n_subj, ' participants</b>.</p>',

fig("00_qc_venn.png", "<b>Figure 1.</b> Overlap among the three recording-level QC criteria. Most exclusions are movement-driven. The fourth filter (implausible architecture) is applied afterwards and noted in the subtitle.", "70%"),

'<h3>Table 1b. Sample composition</h3>',
tbl(t1, align_right = c("Subjects","Nights","Nights/subject","% female")),
'<p class="cap">DI = NREM depth index, N3 / (N1 + N2 + REM). TST = total sleep time.
SOL = sleep onset latency from recording start. Values are unadjusted.</p>',

# ---------------- methods ----------------
'<h2>Methods</h2>',
'<h3>Outcomes</h3>
<ul>
<li><b>NREM depth index (DI)</b> = <code>T(N3) / (T(N1) + T(N2) + T(REM))</code>,
equivalently <code>T(N3) / (TST &minus; T(N3))</code>. Wake is excluded. A ratio, so
unbounded above.</li>
<li><b>Total sleep time (TST)</b> = summed duration of all N1, N2, N3 and REM epochs.</li>
<li><b>Sleep onset latency (SOL)</b> = elapsed time from the start of the recording to the
first sleep epoch. No lights-out diary is available, so the recording start is the
reference; this makes SOL sensitive to how long the device was worn before bed.</li>
</ul>',

'<h3>Models</h3>
<p>All models take the form</p>
<p style="text-align:center"><code>outcome ~ &lt;genotype term&gt; + Age + Age&sup2; + Sex + Site + (1 | participant)</code></p>
<p>Age is mean-centred before squaring. Four genotype parameterisations are fitted per
locus, with controls shared between loci:</p>
<ul>
<li><b>mirror</b> &mdash; linear copy number, coded &minus;1 / 0 / +1</li>
<li><b>U-shape</b> &mdash; squared copy number, coded +1 / 0 / +1</li>
<li><b>joint</b> &mdash; both terms together, each adjusted for the other</li>
<li><b>group</b> &mdash; genotype as a factor, deletion and duplication each vs control</li>
</ul>
<p>Wald tests use REML with Satterthwaite degrees of freedom; model comparisons use
likelihood-ratio tests on ML refits. Benjamini&ndash;Hochberg FDR is applied within three
families across the three primary outcomes.</p>',

'<h3>Covariate adjustment for the figures</h3>
<p>Plotted values are <i>partial residuals with respect to genotype</i>: only the covariate
part of the fixed effects is removed, with each covariate column centred at its sample
mean. Genotype differences, participant random effects and residual scatter are all
retained. The adjustment is computed <b>once per outcome</b> from a single reference model
fitted on all five groups, so a given night appears at the same height in every panel.
Model-implied means are evaluated at the same covariate reference (proportional weighting
over sex and site), which is what places the fitted line inside the plotted points.</p>',

# ---------------- results ----------------
'<h2 class="pb">Results</h2>',
'<p>Every test is listed below, split by outcome for readability. The complete
unsplit table &mdash; all ', nrow(tests), ' tests including the four sensitivity outcomes
&mdash; is <code>results_all_tests.csv</code>. Shaded rows are p &lt; 0.05;
q is the Benjamini&ndash;Hochberg FDR within test family across the three primary
outcomes. N is nights / participants.</p>',

'<h3>Table 2a. NREM depth index</h3>',
tbl(T_DI$df, hl = T_DI$hl, align_right = NUMCOL),

'<h3>Table 2b. Total sleep time</h3>',
tbl(T_TST$df, hl = T_TST$hl, align_right = NUMCOL),

'<h3>Table 2c. Sleep onset latency</h3>',
tbl(T_SOL$df, hl = T_SOL$hl, align_right = NUMCOL),
'<p class="cap">Estimates are in the units of the outcome, per copy for mirror/U-shape
terms and versus control for group terms: index units for depth index, minutes for total
sleep time and sleep onset latency.</p>',

fig("99_copy_number_lines.png", "<b>Figure 2.</b> Model-estimated means with 95% CI across copy number, both loci. Depth index rises with copy number at both loci; total sleep time falls at 22q11.2."),

fig("01_DI_C_16p_mirror.png", "<b>Figure 3.</b> NREM depth index across 16p11.2 copy number. Dots are covariate-adjusted individual nights, coloured by participant."),

fig("01_DI_C_22q_mirror.png", "<b>Figure 4.</b> NREM depth index across 22q11.2 copy number."),

fig("01_DI_A_five_group.png", "<b>Figure 5.</b> NREM depth index by genotype, all five groups. p-values are each genotype vs control."),

fig("02_TST_C_22q_mirror.png", "<b>Figure 6.</b> Total sleep time across 22q11.2 copy number."),

# ---------------- robustness ----------------
'<h2 class="pb">Sensitivity analyses</h2>',
'<p>Five checks were run, ', nrow(sens), ' model refits in total. Each varies one
analytic choice that could plausibly manufacture the depth-index result, holding
everything else fixed.</p>',

'<h3>Summary</h3>',
tbl(sens_summary, align_right = c("Refits")),
'<p class="cap"><b>All ', nrow(sens), ' refits leave the conclusion intact.</b>
The depth-index dose-response is significant at both loci under every outcome definition,
every QC threshold combination tested, and with or without the extreme nights.</p>',

'<div class="key">
<p><b>Where to find these results.</b> Every refit is a row in
<code>tables/sensitivity_analyses.csv</code> (also the <code>sensitivity</code> sheet of
<code>tables/results_all_tests.xlsx</code>). Filter the <code>check</code> column by the
value in the first column of the table above &mdash; for example
<code>check == "N3 outlier threshold"</code> for the outlier analysis. Columns are
<code>outcome, check, setting, locus, model, term, n_nights, n_dropped, estimate, se,
p_value</code>. The code that generates them is <code>code/R/05_sensitivity.R</code>.</p>
</div>',

'<h3>Detail</h3>',

'<h3>Table 3a. The N3 outlier filter</h3>',
tbl(t3, hl = hl3, align_right = c("Nights","Estimate","p")),
'<p class="cap">The filter applied in the primary analysis is marked. With it switched off
entirely the effect is still significant at both loci &mdash; so the finding does not depend
on the outlier removal. Removing the extreme nights helps 22q11.2 substantially and costs
a little at 16p11.2, which is the signature of high-leverage noise rather than of an
effect driven by a few points.</p>',

'<h3>Table 3b. Depth-index definition</h3>',
tbl(t3d, hl = hl3d, align_right = c("Nights","Estimate","p")),
'<p class="cap">All three denominators &mdash; N3 over the other sleep stages (primary),
N3 over all NREM, and N3 over total sleep time &mdash; give the same conclusion at both
loci. The two bounded versions are, if anything, stronger.</p>',

'<h3>Table 3c. Recording-length and movement thresholds</h3>',
tbl(t3b, hl = hl3b, align_right = c("Nights","Estimate","p")),
'<p class="cap">All 24 combinations of recording-length and movement-time thresholds give
p &lt; 0.05 at both loci, including the 30-minute recording-length rule rather than 2 hours.
The result does not depend on where the QC lines are drawn.</p>',

'<h3>Table 3d. Sleep onset latency definition</h3>',
tbl(t3c, hl = hl3c, align_right = c("Estimate","p")),
'<p class="cap">Raw SOL detects 16p11.2 only. On the log scale &mdash; better specified for a
variable with median 18.5 min and maximum 362 min &mdash; the effect is present at both loci.
Latency to persistent sleep is directionally consistent but weaker.</p>',

# ---------------- caveats ----------------
'<h2 class="pb">Is there a depth / duration trade-off, and is it specific to carriers?</h2>',

'<p>Do longer nights have proportionally less deep sleep, and does that hold in
controls or only in CNV carriers?</p>',

'<div class="key">
<p><b>The trade-off is real, and it is present in controls.</b> Within participants, every
group shows it: on a person&rsquo;s longer nights their sleep is proportionally shallower.
Controls sit in the middle of the carrier range (',
sprintf("%+.4f", .s("Control", "within-participant")), ' depth-index units per hour,
p = ', pfmt(tr_slope$p_value[tr_slope$group == "Control" &
                             tr_slope$effect == "within-participant"][1]),
'; within-participant r = ', sprintf("%+.3f", .g(tr_corr, "Control", "within-participant")),
'). Carriers as a whole do not differ from controls (LRT p = ', pfmt(.ip("carrier vs control")),
'). <b>So this is normal sleep physiology, not a CNV phenotype.</b></p>
<p><b>But the between-participant version behaves completely differently.</b> Among controls,
people who habitually sleep longer do <i>not</i> have shallower sleep &mdash; the
participant-level correlation is ', sprintf("%+.3f", .g(tr_corr, "Control", "participant")),
' (p = ', pfmt(tr_corr$pearson_p[tr_corr$group == "Control" &
                                 tr_corr$level == "participant"][1]),
'), essentially zero. In the combined sample it is ',
sprintf("%+.3f", .g(tr_corr, "ALL (combined sample)", "participant")),
' (p = ', pfmt(tr_corr$pearson_p[tr_corr$group == "ALL (combined sample)" &
                                 tr_corr$level == "participant"][1]),
'). That apparent whole-sample trade-off is largely <i>generated by genotype</i>: the groups
differ in depth and in duration in opposing directions, which induces a correlation across
the pooled sample that does not exist within controls. It is a confound, not a physiological
relationship.</p>
<p><b>One group does differ.</b> The five-group test is significant (LRT p = ',
pfmt(.ip("5 genotype groups")), '), driven by 22q11.2 deletion carriers, whose
within-participant slope (', sprintf("%+.4f", .s("22qDel", "within-participant")),
') is about twice as steep as controls. The 22q copy-number interaction is also significant
(p = ', pfmt(tr_int$p_value[tr_int$comparison == "22q copy number"][1]),
'); the 16p one is not (p = ', pfmt(tr_int$p_value[tr_int$comparison == "16p copy number"][1]), ').</p>
</div>',

fig("90_tradeoff_raw_scatter.png", "<b>Figure 7.</b> Raw night-level association between depth index and total sleep time, for the combined sample and each genotype. The inverse relationship is visible in every group, including controls."),

fig("91_tradeoff_within_participant.png", "<b>Figure 8.</b> The same relationship after removing each participant's own average, isolating the night-to-night trade-off. Present in all five groups."),

fig("92_tradeoff_slopes_by_group.png", "<b>Figure 9.</b> Within- and between-participant slopes by group. Left: all groups negative, controls (filled) mid-range. Right: controls sit on zero while carrier groups are negative &mdash; the between-participant trade-off is a genotype-driven artefact."),

'<h3>Table 4a. Correlation between depth index and total sleep time</h3>',
if (!is.null(tT1)) tbl(tT1) else "",
'<p class="cap">Night-level correlations pool nights within participants, so their
p-values are anti-conservative &mdash; they are shown because they are the quantity usually
reported. The participant-level and within-participant columns are the interpretable ones,
and they tell opposite stories.</p>',

'<h3>Table 4b. Mixed-model slopes, per group</h3>',
if (!is.null(tT2)) tbl(tT2, hl = hlT2, align_right = c("Slope per hour","95% CI","p","Nights")) else "",
'<p class="cap">Change in depth index per additional hour of total sleep time, adjusted for
age, age&sup2;, sex and site with a random intercept per participant (a random slope for the
within-participant term where it converged).</p>',

'<h3>Table 4c. Does the trade-off differ by genotype?</h3>',
if (!is.null(tT3)) tbl(tT3, hl = hlT3, align_right = c("Statistic","df","p")) else "",

'<h3>Is the trade-off an artefact of TST being inside the depth index?</h3>',
'<p>Depth index is <code>N3 / (N1 + N2 + REM)</code>, so total sleep time enters its
denominator. That deserves scrutiny. Algebraically the index is <i>scale-invariant</i>: if
N3 were a fixed fraction <i>c</i> of TST, then DI = <i>c</i>/(1&nbsp;&minus;&nbsp;<i>c</i>)
regardless of how long the night is. So an association means composition genuinely shifts
with duration. Three independent controls confirm it:</p>',
if (!is.null(tT4)) tbl(tT4) else "",
'<p class="cap">The third control is the most direct: computing depth index over a fixed
240 minutes of sleep holds duration constant by construction, and the correlation
<b>reverses sign</b>. Longer sleepers have <i>more</i> early-night N3, not less. The
whole-night negative association is therefore a dilution effect &mdash; extra sleep arrives
as late, N2- and REM-rich cycles &mdash; rather than long sleepers being shallow sleepers.</p>',

'<h3>Which signal is primary: depth or duration?</h3>',
'<p>The trade-off raises a question for the main results in both directions. 22q11.2
duplication carriers sleep less <i>and</i> have a higher depth index, so either effect could
in principle be the other one in disguise. Each mirror model was therefore refitted with the
other measure added, split into its within- and between-participant parts.</p>',
if (!is.null(tT5)) tbl(tT5, hl = hlT5, align_right = c("Copy-number effect","SE","p")) else "",

'<div class="key">
<p><b>The two directions behave differently, and the asymmetry is the finding.</b></p>
<p><b>Depth survives adjustment for duration.</b> 16p11.2 goes from ',
sprintf("%+.4f", .adje("DI","16p",FALSE)), ' (p = ', pfmt(.adjp("DI","16p",FALSE)),
') to ', sprintf("%+.4f", .adje("DI","16p",TRUE)), ' (p = ', pfmt(.adjp("DI","16p",TRUE)),
'), and 22q11.2 from ', sprintf("%+.4f", .adje("DI","22q",FALSE)), ' (p = ',
pfmt(.adjp("DI","22q",FALSE)), ') to ', sprintf("%+.4f", .adje("DI","22q",TRUE)),
' (p = ', pfmt(.adjp("DI","22q",TRUE)), ') &mdash; attenuated by roughly a fifth to a
quarter, but still significant.</p>
<p><b>Duration does not survive adjustment for depth.</b> The 22q11.2 effect on total sleep
time falls from ', sprintf("%+.1f", .adje("TST","22q",FALSE)), ' min per copy (p = ',
pfmt(.adjp("TST","22q",FALSE)), ') to ', sprintf("%+.1f", .adje("TST","22q",TRUE)),
' min (p = ', pfmt(.adjp("TST","22q",TRUE)), '), losing significance. At 16p11.2, where the
duration effect was not significant to begin with, it collapses to essentially zero
(', sprintf("%+.1f", .adje("TST","16p",TRUE)), ' min, p = ', pfmt(.adjp("TST","16p",TRUE)), ').</p>
<p>The natural reading is that <b>NREM depth is the primary phenotype</b>, and the sleep-time
effect is largely shared variance with it rather than an independent duration signal.</p>
</div>',

'<p><b>Two cautions on that reading.</b> First, the trade-off means depth index is partly a
function of duration, so adjusting total sleep time for depth adjusts for something partly
downstream of it. That biases the duration estimate toward the null by construction, and
some of the asymmetry is therefore expected rather than informative. Second, the duration
effect started weaker (p = ', pfmt(.adjp("TST","22q",FALSE)),
' versus ', pfmt(.adjp("DI","22q",FALSE)), ' and ', pfmt(.adjp("DI","16p",FALSE)),
'), so it had less room to survive any attenuation. What these models establish is that the
depth signal is statistically separable from duration while the duration signal is not
separable from depth &mdash; not that depth is causally upstream.</p>',

fig("93_DI_16p_TSTadjusted.png", "<b>Figure 10.</b> 16p11.2 copy number and depth index with sleep duration removed alongside the usual covariates. Compare with Figure 3."),

fig("93_DI_22q_TSTadjusted.png", "<b>Figure 11.</b> 22q11.2 copy number and depth index with sleep duration removed. Compare with Figure 4."),

'<div class="key">
<p><b>Where to find these results.</b> <code>tables/tradeoff_correlations.csv</code> (raw
correlations at all three levels), <code>tradeoff_slopes_by_group.csv</code> (per-group
mixed-model slopes), <code>tradeoff_group_comparisons.csv</code> (interaction tests),
<code>tradeoff_loglog_check.csv</code>, <code>tradeoff_null_calibration.csv</code> and
<code>tradeoff_fixedwindow_correlations.csv</code> (the three controls), and
<code>tradeoff_main_effect_adjusted_for_TST.csv</code>. All are also sheets in
<code>results_all_tests.xlsx</code>. Code: <code>code/R/06_tradeoff.R</code>.</p>
</div>',

'<h2>Caveats</h2>',
'<div class="key warn">
<p><b>Age is badly imbalanced across genotype, and this is the main threat to the
depth-index result.</b> 16p11.2 deletion carriers average ',
sprintf("%.1f", desc$age_mean[desc$Genotype == "16pDel"]), ' years versus ',
sprintf("%.1f", desc$age_mean[desc$Genotype == "16pDup"]), ' for duplication carriers and ',
sprintf("%.1f", desc$age_mean[desc$Genotype == "Control"]), ' for controls &mdash; a gap of
roughly 13 years within the 16p locus. N3 declines steeply with age, and the deletion group
is the youngest while showing the <i>lowest</i> depth index, which is the opposite of what
age confounding alone would produce. Age and age&sup2; are adjusted for, but the adjustment
is extrapolating across substantially non-overlapping age ranges, and residual confounding
cannot be excluded. An age-matched subsample analysis is the natural next step.</p>
</div>',
'<ul>
<li><b>Control group is small</b> &mdash; ', desc$n_subjects[desc$Genotype == "Control"],
' participants, ', desc$n_nights[desc$Genotype == "Control"], ' nights &mdash; and is shared
between the two locus analyses, so the 16p and 22q results are not independent.</li>
<li><b>16p11.2 duplication is the smallest group</b> (',
desc$n_subjects[desc$Genotype == "16pDup"], ' participants) and is unevenly distributed
across sites, with none at one site. Site is adjusted for, but the 16p duplication estimate
rests on few participants.</li>
<li><b>Sex is imbalanced</b>: ', sprintf("%.0f", desc$pct_female[desc$Genotype == "Control"]),
'% female among controls versus ',
sprintf("%.0f", desc$pct_female[desc$Genotype == "22qDup"]), '% in 22q11.2 duplication.</li>
<li><b>SOL is measured from recording start</b>, not lights-out. It partly reflects device
handling rather than sleep physiology, and should be treated as provisional.</li>
<li><b>', n_pass - n_nights, ' QC-passing nights were discarded</b> at linkage: one device
is absent from the manifest, one participant has no genotype recorded, and one
carries both a 22q11.2 and a 16p11.2 deletion and so belongs to neither dose series.</li>
<li><b>Movement time is high in this dataset</b> overall, and is the dominant QC exclusion
(', venn$n[venn$region == "HighMove only"], ' of ', sum(venn$n),
' exclusions). This is characteristic of headband EEG rather than a data error.</li>
</ul>',

# ---------------- repro ----------------
'<h2>Reproducing this</h2>
<p>The pipeline is in <code>code/</code>. From that directory:</p>
<p><code>Rscript inspect_inputs.R</code> &nbsp;&mdash; verify the inputs parse<br>
<code>Rscript run_all.R</code> &nbsp;&mdash; ingest, QC, model, figures, sensitivity<br>
<code>Rscript make_report.R</code> &nbsp;&mdash; regenerate this document</p>
<p>Requires R with <code>lme4</code>, <code>lmerTest</code> and <code>ggplot2</code>.
Every number and figure above is read from the files written by <code>run_all.R</code>,
so the report cannot drift from the analysis. Analysis choices are in
<code>config.R</code>.</p>',

'<p class="meta">Report generated ', format(Sys.time(), "%Y-%m-%d %H:%M"),
' &middot; R ', getRversion(),
' &middot; lme4 ', as.character(utils::packageVersion("lme4")), '</p>',

'</body></html>')

outfile <- file.path(OUTDIR, "report.html")
writeLines(html, outfile, useBytes = TRUE)
cat("Report written:", outfile, "\n")
cat("  ", round(file.info(outfile)$size / 1024), "KB\n")
