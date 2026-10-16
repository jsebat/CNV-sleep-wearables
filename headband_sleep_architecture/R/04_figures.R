# =====================================================================
# 04_figures.R -- QC Venn, box/scatter panels, copy-number line plots
#
# Every panel plots the SAME covariate-adjusted nightly value
# (`adj_<outcome>`, computed once in 03_models.R), so a given night sits at
# the same height in the single-genotype figure, the mirror figure and the
# joint figure.  Only the overlaid model line changes between panels.
# =====================================================================

suppressPackageStartupMessages(library(ggplot2))

GENO_PLOT_ORDER <- c("16pDel", "22qDel", "Control", "22qDup", "16pDup")
GENO_PLOT_LABEL <- c("16p Deletion", "22q Deletion", "Control",
                     "22q Duplication", "16p Duplication")

.base_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.title     = element_text(face = "bold", size = base_size + 3),
      plot.subtitle  = element_text(size = base_size - 1, margin = margin(b = 10)),
      axis.title.x   = element_text(face = "bold"),
      axis.title.y   = element_text(face = "bold"),
      axis.text.x    = element_text(face = "bold"),
      plot.caption   = element_text(size = base_size - 4, hjust = 0,
                                    lineheight = 1.2, margin = margin(t = 14)),
      plot.margin    = margin(15, 25, 15, 15)
    )
}

#' Pastel colour per participant.
#'
#' The palette is permuted with a fixed seed: subject IDs are usually sorted
#' by genotype, so assigning the ramp in ID order would give each group a
#' single hue and make the figure look group-coloured.
.participant_colors <- function(ids) {
  ids <- unique(as.character(ids))
  n <- max(length(ids), 2)
  pal <- grDevices::hcl.colors(n, "Pastel 1")
  ord <- withr_seed(1, sample.int(n))
  stats::setNames(pal[ord][seq_along(ids)], ids)
}

#' Evaluate `expr` under a fixed RNG seed, restoring the caller's RNG state.
withr_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  expr
}

#' Robust y-axis limits: Tukey fences at `mult` IQRs, never tighter than the
#' range of `must_show` (the model estimates and their CIs, which must stay
#' visible), and never wider than the data itself.
#'
#' Applied with coord_cartesian(), so points outside are hidden but retained:
#' boxplots, model fits and p-values all still use every observation.
robust_ylim <- function(values, must_show = numeric(0), mult = 3, pad = 0.04) {
  v <- values[is.finite(values)]
  if (!length(v)) return(NULL)
  if (!is.finite(mult)) return(NULL)
  q <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
  iqr <- q[2] - q[1]
  if (iqr <= 0) return(NULL)
  lo <- max(min(v), q[1] - mult * iqr)
  hi <- min(max(v), q[2] + mult * iqr)
  ms <- must_show[is.finite(must_show)]
  if (length(ms)) { lo <- min(lo, min(ms)); hi <- max(hi, max(ms)) }
  if (!(hi > lo)) return(NULL)
  span <- hi - lo
  c(lo - pad * span, hi + pad * span)
}

#' How many values fall outside the drawn axis, as a caption line.
clip_note <- function(values, ylim) {
  if (is.null(ylim)) return("")
  n <- sum(is.finite(values) & (values < ylim[1] | values > ylim[2]))
  if (!n) return("")
  sprintf(paste("\n%d night%s fall outside the plotted axis and %s not drawn;",
                "all data are retained in the boxplots, the model fit and every p-value."),
          n, if (n == 1) "" else "s", if (n == 1) "is" else "are")
}

CAPTION <- paste(
  "Coloured dots are covariate-adjusted individual nights; one colour per participant.",
  "Adjustment removes age, age², sex and site only — genotype differences and",
  "between-participant variation are retained, and are identical across all panels.",
  "Open circles / line: model-estimated means with 95% CI at the mean covariate profile.",
  "Mixed model: age + age² + sex + site, random intercept per participant.",
  sep = "\n")

# ---------------------------------------------------------------------
# QC Venn diagram (base graphics -- no extra package needed)
# ---------------------------------------------------------------------
plot_qc_venn <- function(metrics, file, n_entered = NULL,
                         min_trt_min = NA, max_mt_frac = NA) {
  a <- metrics$excl_no_sleep
  b <- metrics$excl_too_short
  cc <- metrics$excl_high_movement
  n_entered <- n_entered %||% nrow(metrics)
  n_excl <- sum(a | b | cc)
  # a 4th criterion (implausible architecture) is applied after these three;
  # it is reported in the subtitle rather than drawn, to keep the Venn legible
  n_extra <- if ("excl_implausible_n3" %in% names(metrics))
    sum(metrics$excl_implausible_n3 & metrics$qc_pass_core) else 0L

  grDevices::png(file, width = 2000, height = 2000, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(4, 2, 6, 2)); on.exit(graphics::par(op), add = TRUE)

  graphics::plot(NA, xlim = c(-1.9, 1.9), ylim = c(-2.2, 1.9),
                 asp = 1, axes = FALSE, xlab = "", ylab = "")
  r <- 1.15
  cx <- c(-0.62, 0.62, 0); cy <- c(0.55, 0.55, -0.55)
  fills <- c("#F2F9FB", "#EAF4F8", "#AFD9E4")
  for (i in 1:3) {
    graphics::symbols(cx[i], cy[i], circles = r, inches = FALSE, add = TRUE,
                      fg = "black", bg = grDevices::adjustcolor(fills[i], alpha.f = 0.45),
                      lwd = 2)
  }
  tx <- function(x, y, s) graphics::text(x, y, s, cex = 1.15)
  tx(-0.95,  0.80, sum(a & !b & !cc))
  tx( 0.95,  0.80, sum(!a & b & !cc))
  tx( 0.00, -1.25, sum(!a & !b & cc))
  tx( 0.00,  0.95, sum(a & b & !cc))
  tx(-0.72, -0.35, sum(a & !b & cc))
  tx( 0.72, -0.35, sum(!a & b & cc))
  tx( 0.00, -0.05, sum(a & b & cc))

  lab_trt <- if (is.na(min_trt_min)) "TRT below threshold" else
    if (min_trt_min >= 60 && min_trt_min %% 60 == 0)
      sprintf("TRT < %g hours", min_trt_min / 60) else
      sprintf("TRT < %g min", min_trt_min)
  lab_mt <- if (is.na(max_mt_frac)) "MT above threshold" else
    sprintf("MT > %g%%", 100 * max_mt_frac)

  graphics::text(-1.35, 1.90, "No Sleep\n(TST = 0)", cex = 1.05)
  graphics::text( 1.35, 1.90, paste0("Too Short\n(", lab_trt, ")"), cex = 1.05)
  graphics::text( 0.00, -2.05, paste0("High Movement\n(", lab_mt, ")"), cex = 1.05)
  graphics::mtext("Sleep EEG QC Exclusion Overlap", side = 3, line = 3.4,
                  font = 2, cex = 1.5)
  graphics::mtext(sprintf("%d recordings entered QC | %d excluded by these three | %d retained",
                          n_entered, n_excl, n_entered - n_excl),
                  side = 3, line = 1.9, cex = 1.0)
  if (n_extra > 0) {
    graphics::mtext(sprintf("a further %d excluded for implausible architecture (N3 > 50%% of TST); %d retained overall",
                            n_extra, sum(metrics$qc_pass)),
                    side = 3, line = 0.7, cex = 0.85, col = "grey35")
  }
  invisible(file)
}

# ---------------------------------------------------------------------
# Panel A: all five genotypes vs Control
# ---------------------------------------------------------------------
plot_five_group <- function(d, fits, outcome, label, tests) {
  s <- fits[["all_data"]]
  fit <- fits[["all_group"]]
  adj <- paste0("adj_", outcome)
  s <- s[!is.na(s[[adj]]), , drop = FALSE]
  if (!nrow(s) || is.null(fit)) return(NULL)

  s$grp <- factor(as.character(s$Genotype), levels = GENO_PLOT_ORDER,
                  labels = GENO_PLOT_LABEL)

  gcols <- paste0("Genotype", setdiff(GENO_PLOT_ORDER, "Control"))
  gcols <- intersect(gcols, colnames(lme4::getME(fit, "X")))
  targets <- as.data.frame(matrix(0, nrow = length(GENO_PLOT_ORDER),
                                  ncol = length(gcols)))
  names(targets) <- gcols
  for (i in seq_along(GENO_PLOT_ORDER)) {
    cn <- paste0("Genotype", GENO_PLOT_ORDER[i])
    if (cn %in% gcols) targets[i, cn] <- 1
  }
  mm <- model_means(fit, "Genotype", targets)
  mm$grp <- factor(GENO_PLOT_LABEL, levels = GENO_PLOT_LABEL)

  tt <- tests[tests$outcome == outcome & tests$cohort == "All5" &
              tests$model == "group_vs_control" & tests$test != "LRT", ]
  tt$grp <- factor(GENO_PLOT_LABEL[match(sub(" vs Control$", "", tt$term_label),
                                         GENO_PLOT_ORDER)],
                   levels = GENO_PLOT_LABEL)
  tt <- tt[!is.na(tt$grp), , drop = FALSE]

  # axis clipped to the bulk of the data; model means and CIs must stay visible
  yl <- robust_ylim(s[[adj]], c(mm$lower, mm$upper), mult = FIG_CLIP_IQR)
  if (is.null(yl)) yl <- range(s[[adj]], na.rm = TRUE)
  # leave headroom for the p-value labels, then place them inside it
  yl[2] <- yl[2] + 0.10 * diff(yl)
  tt$y <- yl[2] - 0.03 * diff(yl)

  omni <- tests[tests$outcome == outcome & tests$cohort == "All5" &
                tests$test == "LRT", ]

  ggplot(s, aes(x = grp, y = .data[[adj]])) +
    geom_jitter(aes(colour = GUID), width = 0.13, height = 0, size = 1.7, alpha = 0.6) +
    geom_boxplot(width = 0.46, outlier.shape = NA, fill = NA,
                 colour = "black", linewidth = 0.7) +
    geom_errorbar(data = mm, aes(x = grp, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.8) +
    geom_point(data = mm, aes(x = grp, y = estimate), inherit.aes = FALSE,
               shape = 21, size = 4, fill = "white", colour = "black", stroke = 1) +
    { if (nrow(tt)) geom_text(data = tt, aes(x = grp, y = y, label = paste0("p = ", fmt_p(p_value))),
                              inherit.aes = FALSE, size = 3.4, fontface = "bold") } +
    scale_colour_manual(values = .participant_colors(s$GUID)) +
    coord_cartesian(ylim = yl) +
    labs(
      title = paste0("Single-genotype effects on ", label),
      subtitle = if (nrow(omni))
        sprintf("Each CNV genotype vs Control  |  omnibus LRT p = %s", fmt_p(omni$p_value[1]))
        else "Each CNV genotype vs Control",
      x = "Genotype", y = paste0("Covariate-adjusted ", label),
      caption = paste0(CAPTION, clip_note(s[[adj]], yl))
    ) +
    .base_theme()
}

# ---------------------------------------------------------------------
# Panel B: one locus, Del -> Control -> Dup, with the joint model line
# ---------------------------------------------------------------------
plot_locus <- function(d, fits, outcome, label, tests, locus,
                       which_model = c("joint", "mirror", "ushape", "group")) {
  which_model <- match.arg(which_model)
  s   <- fits[[paste0(locus, "_data")]]
  fit <- fits[[paste0(locus, "_",
                      switch(which_model, joint = "joint", mirror = "mirror",
                             ushape = "ushape", group = "group"))]]
  if (is.null(s) || is.null(fit)) return(NULL)

  adj <- paste0("adj_", outcome)
  s <- s[!is.na(s[[adj]]), , drop = FALSE]
  if (!nrow(s)) return(NULL)

  lv <- c(paste0(locus, "Del"), "Control", paste0(locus, "Dup"))
  lb <- c(paste0(substr(locus, 1, 3), " Deletion"), "Control",
          paste0(substr(locus, 1, 3), " Duplication"))
  s$grp <- factor(as.character(s$Genotype), levels = lv, labels = lb)

  # model-implied means at copy number 1 / 2 / 3
  if (which_model == "group") {
    gcols <- intersect(paste0("Genotype", lv[c(1, 3)]),
                       colnames(lme4::getME(fit, "X")))
    targets <- as.data.frame(matrix(0, nrow = 3, ncol = length(gcols)))
    names(targets) <- gcols
    for (i in c(1, 3)) {
      cn <- paste0("Genotype", lv[i])
      if (cn %in% gcols) targets[i, cn] <- 1
    }
    mm <- model_means(fit, "Genotype", targets)
    geno_terms <- "Genotype"
  } else {
    targets <- switch(
      which_model,
      joint  = data.frame(copy_c = c(-1, 0, 1), copy_q = c(1, 0, 1)),
      mirror = data.frame(copy_c = c(-1, 0, 1)),
      ushape = data.frame(copy_q = c(1, 0, 1))
    )
    geno_terms <- names(targets)
    mm <- model_means(fit, geno_terms, targets)
  }
  mm$grp <- factor(lb, levels = lb)

  yl <- robust_ylim(s[[adj]], c(mm$lower, mm$upper), mult = FIG_CLIP_IQR)

  tt <- tests[tests$outcome == outcome & tests$cohort == locus, ]
  sub_bits <- character(0)
  pm <- tt$p_value[tt$model == "mirror_single" & tt$term == "copy_c"]
  pu <- tt$p_value[tt$model == "ushape_single" & tt$term == "copy_q"]
  pjc <- tt$p_value[tt$model == "joint" & tt$term == "copy_c"]
  pjq <- tt$p_value[tt$model == "joint" & tt$term == "copy_q"]
  if (length(pm))  sub_bits <- c(sub_bits, paste0("mirror p = ", fmt_p(pm[1])))
  if (length(pu))  sub_bits <- c(sub_bits, paste0("U-shape p = ", fmt_p(pu[1])))
  if (which_model == "joint" && length(pjc) && length(pjq)) {
    sub_bits <- c(sub_bits,
                  paste0("joint: mirror p = ", fmt_p(pjc[1]),
                         ", U p = ", fmt_p(pjq[1])))
  }

  ggplot(s, aes(x = grp, y = .data[[adj]])) +
    geom_jitter(aes(colour = GUID), width = 0.12, height = 0, size = 1.9, alpha = 0.6) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA,
                 colour = "black", linewidth = 0.7) +
    geom_line(data = mm, aes(x = grp, y = estimate, group = 1),
              inherit.aes = FALSE, linewidth = 0.9) +
    geom_errorbar(data = mm, aes(x = grp, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.8) +
    geom_point(data = mm, aes(x = grp, y = estimate), inherit.aes = FALSE,
               shape = 21, size = 4, fill = "white", colour = "black", stroke = 1) +
    scale_colour_manual(values = .participant_colors(s$GUID)) +
    coord_cartesian(ylim = yl) +
    labs(
      title = sprintf("%s copy number and %s", substr(locus, 1, 3), label),
      subtitle = paste(sub_bits, collapse = "   |   "),
      x = "Copy number group", y = paste0("Covariate-adjusted ", label),
      caption = paste0(CAPTION, clip_note(s[[adj]], yl))
    ) +
    .base_theme()
}

# ---------------------------------------------------------------------
# Summary line plot: model-implied means across copy number, both loci,
# all primary outcomes.
# ---------------------------------------------------------------------
plot_copy_number_lines <- function(all_fits, outcomes, tests) {
  rows <- list()
  for (o in outcomes) {
    if (!isTRUE(o$primary)) next
    fits <- all_fits[[o$name]]
    if (is.null(fits)) next
    for (locus in c("16p", "22q")) {
      fit <- fits[[paste0(locus, "_group")]]
      if (is.null(fit)) next
      lv <- c(paste0(locus, "Del"), "Control", paste0(locus, "Dup"))
      gcols <- intersect(paste0("Genotype", lv[c(1, 3)]),
                         colnames(lme4::getME(fit, "X")))
      targets <- as.data.frame(matrix(0, nrow = 3, ncol = length(gcols)))
      names(targets) <- gcols
      for (i in c(1, 3)) {
        cn <- paste0("Genotype", lv[i])
        if (cn %in% gcols) targets[i, cn] <- 1
      }
      mm <- model_means(fit, "Genotype", targets)
      mm$CopyNumber <- c(1, 2, 3)
      mm$Locus   <- paste0(locus, "11.2")
      mm$Outcome <- o$label
      rows[[length(rows) + 1]] <- mm[, c("CopyNumber", "Locus", "Outcome",
                                         "estimate", "lower", "upper")]
    }
  }
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  df$Outcome <- factor(df$Outcome,
                       levels = vapply(outcomes[vapply(outcomes, function(o) isTRUE(o$primary), logical(1))],
                                       function(o) o$label, character(1)))

  ggplot(df, aes(x = CopyNumber, y = estimate, colour = Locus, group = Locus)) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08, linewidth = 0.7) +
    geom_point(size = 3.2, shape = 21, fill = "white", stroke = 1.1) +
    facet_wrap(~ Outcome, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:3,
                       labels = c("1\n(del)", "2\n(control)", "3\n(dup)"),
                       expand = expansion(mult = 0.18)) +
    scale_colour_manual(values = c("16p11.2" = "#1F6FB4", "22q11.2" = "#C1442E")) +
    labs(
      title = "Model-estimated sleep measures across copy number",
      subtitle = "Estimated means with 95% CI at the mean covariate profile",
      x = "Copy number", y = "Model-estimated value", colour = "Locus",
      caption = "Mixed model: age + age² + sex + site, random intercept per participant."
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.caption = element_text(size = 9, hjust = 0)
    )
}

# ---------------------------------------------------------------------
# Depth / duration trade-off
# ---------------------------------------------------------------------

#' Within-participant trade-off: person-mean-centred DI against
#' person-mean-centred TST, one facet per genotype group.
#' Centring removes stable between-person differences, so each panel shows
#' the night-to-night relationship only.
plot_tradeoff_within <- function(d, outcome = "DI",
                                 label = "NREM depth index",
                                 groups = TRADEOFF_GROUPS) {
  s <- d[!is.na(d[[outcome]]) & !is.na(d$TST), , drop = FALSE]
  s <- s[as.character(s$Genotype) %in% groups, , drop = FALSE]
  if (!nrow(s)) return(NULL)

  ymean <- tapply(s[[outcome]], as.character(s$GUID), mean)
  tmean <- tapply(s$TST,        as.character(s$GUID), mean)
  s$y_w <- s[[outcome]] - as.numeric(ymean[as.character(s$GUID)])
  s$t_w <- (s$TST       - as.numeric(tmean[as.character(s$GUID)])) / 60
  s <- s[as.character(s$GUID) %in% names(which(table(as.character(s$GUID)) > 1)), ]
  if (!nrow(s)) return(NULL)

  s$grp <- factor(as.character(s$Genotype), levels = groups,
                  labels = c("Control", "16p deletion", "16p duplication",
                             "22q deletion", "22q duplication")[
                               match(groups, TRADEOFF_GROUPS)])

  yl <- robust_ylim(s$y_w, mult = FIG_CLIP_IQR)
  xl <- robust_ylim(s$t_w, mult = FIG_CLIP_IQR)

  ggplot(s, aes(x = t_w, y = y_w)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_point(aes(colour = GUID), size = 1.5, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                fill = "grey70", linewidth = 0.9) +
    facet_wrap(~ grp, nrow = 1) +
    scale_colour_manual(values = .participant_colors(s$GUID)) +
    coord_cartesian(ylim = yl, xlim = xl) +
    labs(
      title = sprintf("Within-participant trade-off between %s and sleep duration", label),
      subtitle = "Each point is one night, expressed as a deviation from that participant's own average",
      x = "Total sleep time, deviation from participant mean (hours)",
      y = sprintf("%s,\ndeviation from participant mean", label),
      caption = paste(
        "Person-mean centring removes stable between-participant differences, so each panel shows only the night-to-night relationship.",
        "A downward slope means that on a participant's longer nights, their sleep is proportionally less deep.",
        "Line is an ordinary least-squares fit with 95% CI; the tests in the results table use mixed models with a random intercept per participant.",
        sep = "\n")
    ) +
    .base_theme(12) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

#' Forest plot of within- and between-participant slopes by group.
plot_tradeoff_slopes <- function(by_group, label = "NREM depth index") {
  if (is.null(by_group) || !nrow(by_group)) return(NULL)
  b <- by_group
  b$grp <- factor(b$group, levels = rev(TRADEOFF_GROUPS))
  b$eff <- factor(b$effect,
                  levels = c("within-participant", "between-participant"),
                  labels = c("Within participant\n(night to night)",
                             "Between participants\n(who sleeps longer overall)"))

  ggplot(b, aes(x = estimate, y = grp)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y",
                  width = 0.18, linewidth = 0.8) +
    geom_point(aes(fill = grp == "Control"), shape = 21, size = 3.6,
               colour = "black", stroke = 0.9) +
    scale_fill_manual(values = c(`TRUE` = "#1F6FB4", `FALSE` = "white"),
                      guide = "none") +
    facet_wrap(~ eff, scales = "free_x") +
    labs(
      title = sprintf("%s per additional hour of sleep", label),
      subtitle = "Mixed model per group: age + age² + sex + site, random intercept per participant",
      x = sprintf("Change in %s per hour of total sleep time (95%% CI)", label),
      y = NULL,
      caption = paste(
        "Controls are filled. Intervals crossing zero indicate no detectable relationship in that group.",
        "The two panels answer different questions and can differ in sign.",
        sep = "\n")
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 15),
          plot.caption = element_text(size = 8.5, hjust = 0, lineheight = 1.3),
          axis.text.y = element_text(face = "bold"))
}

#' Raw DI-against-TST scatter: does the inverse correlation exist at all?
#' Left facet is the combined sample, then one facet per genotype group.
#' Annotated with the night-level Pearson r so the reader can see immediately
#' whether it is present in controls as well as carriers.
plot_tradeoff_scatter <- function(d, corr, outcome = "DI",
                                  label = "NREM depth index",
                                  groups = TRADEOFF_GROUPS) {
  s <- d[!is.na(d[[outcome]]) & !is.na(d$TST), , drop = FALSE]
  if (!nrow(s)) return(NULL)
  labs_ <- c(Control = "Control", `16pDel` = "16p deletion",
             `16pDup` = "16p duplication", `22qDel` = "22q deletion",
             `22qDup` = "22q duplication")
  panels <- rbind(
    transform(s, panel = "Combined sample"),
    transform(s[as.character(s$Genotype) %in% groups, ],
              panel = unname(labs_[as.character(
                s$Genotype[as.character(s$Genotype) %in% groups])])))
  panels$panel <- factor(panels$panel,
                         levels = c("Combined sample", unname(labs_[groups])))

  key <- c(`Combined sample` = "ALL (combined sample)",
           stats::setNames(groups, unname(labs_[groups])))
  ann <- do.call(rbind, lapply(levels(panels$panel), function(p) {
    r <- corr[corr$group == key[[p]] & corr$level == "night", ]
    if (!nrow(r)) return(NULL)
    data.frame(panel = factor(p, levels = levels(panels$panel)),
               lab = sprintf("r = %+.3f\np = %s", r$pearson_r, fmt_p(r$pearson_p)),
               stringsAsFactors = FALSE)
  }))

  yl <- robust_ylim(panels[[outcome]], mult = FIG_CLIP_IQR)
  xl <- range(panels$TST, na.rm = TRUE)

  ggplot(panels, aes(x = TST, y = .data[[outcome]])) +
    geom_point(aes(colour = GUID), size = 1.3, alpha = 0.55) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                fill = "grey70", linewidth = 0.9) +
    geom_text(data = ann, aes(x = xl[1] + 0.03 * diff(xl),
                              y = yl[2] - 0.02 * diff(yl), label = lab),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2,
              fontface = "bold", lineheight = 1.05) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = .participant_colors(panels$GUID)) +
    coord_cartesian(ylim = yl) +
    labs(
      title = sprintf("Raw association between %s and total sleep time", label),
      subtitle = "Every night plotted; the leftmost panel repeats the whole sample, then one panel per genotype",
      x = "Total sleep time (min)", y = label,
      caption = paste(
        "Unadjusted night-level values and ordinary least-squares fits, shown to make the raw association visible.",
        "These r values pool nights within participants, so their p-values are anti-conservative; the mixed models in the",
        "results table separate the within- and between-participant components properly.",
        sep = "\n")
    ) +
    .base_theme(12) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

#' Locus panel using the TST-ADJUSTED depth index: the copy-number effect with
#' sleep duration removed alongside the usual covariates.
plot_locus_tst_adjusted <- function(d, outcome, label, locus,
                                    adj_col = paste0("adjTST_", outcome)) {
  s <- locus_subset(d, locus)
  s <- add_person_centered(s, "TST")
  s <- s[!is.na(s[[adj_col]]) & !is.na(s$TST_w), , drop = FALSE]
  if (!nrow(s)) return(NULL)

  fit <- fit_lmm(build_formula(outcome, c("copy_c", "TST_w", "TST_b"), s), s)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  mm <- model_means(fit, "copy_c", data.frame(copy_c = c(-1, 0, 1)))

  lb <- c(paste0(substr(locus, 1, 3), " Deletion"), "Control",
          paste0(substr(locus, 1, 3), " Duplication"))
  s$grp  <- factor(as.character(s$Genotype),
                   levels = c(paste0(locus, "Del"), "Control", paste0(locus, "Dup")),
                   labels = lb)
  mm$grp <- factor(lb, levels = lb)

  yl <- robust_ylim(s[[adj_col]], c(mm$lower, mm$upper), mult = FIG_CLIP_IQR)

  ggplot(s, aes(x = grp, y = .data[[adj_col]])) +
    geom_jitter(aes(colour = GUID), width = 0.12, height = 0, size = 1.9, alpha = 0.6) +
    geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA,
                 colour = "black", linewidth = 0.7) +
    geom_line(data = mm, aes(x = grp, y = estimate, group = 1),
              inherit.aes = FALSE, linewidth = 0.9) +
    geom_errorbar(data = mm, aes(x = grp, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.8) +
    geom_point(data = mm, aes(x = grp, y = estimate), inherit.aes = FALSE,
               shape = 21, size = 4, fill = "white", colour = "black", stroke = 1) +
    scale_colour_manual(values = .participant_colors(s$GUID)) +
    coord_cartesian(ylim = yl) +
    labs(
      title = sprintf("%s copy number and %s, adjusted for sleep duration",
                      substr(locus, 1, 3), label),
      subtitle = sprintf("Mirror effect after removing total sleep time: p = %s",
                         fmt_p(co["copy_c", "Pr(>|t|)"])),
      x = "Copy number group",
      y = paste0("Adjusted ", label, "\n(age, age², sex, site and sleep duration removed)"),
      caption = paste(
        paste("Adjustment removes age, age², sex, site AND total sleep time (within- and between-participant),",
              "\nbut not genotype. Compare with the unadjusted panel: if the copy-number effect persists here, it is",
              "\nnot merely a consequence of the depth / duration trade-off.", sep = ""),
        clip_note(s[[adj_col]], yl), sep = "")
    ) +
    .base_theme()
}

save_fig <- function(p, file, w = FIG_WIDTH, h = FIG_HEIGHT, dpi = FIG_DPI) {
  if (is.null(p)) return(invisible(NULL))
  suppressWarnings(ggsave(file, p, width = w, height = h, dpi = dpi, bg = "white"))
  invisible(file)
}
