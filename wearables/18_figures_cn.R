# Report figures for the copy-number results.
# Adjusted values = partial residual w.r.t. genotype from ONE reference model per outcome, fitted on the
# full five-group cohort:  y ~ Genotype + Age_c + I(Age_c^2) + Sex + Site + (1 | GUID)  (REML).
# Covariate design columns are centred, so adjusted values keep the raw scale and mean; genotype is in the
# model but not subtracted.  Points = participant means of adjusted nightly values; diamonds = model means
# (covariates at sample means) with 95 % CI.  Annotation = mirror (linear copy-number) estimate from
# results_all_tests.csv.  Colours: deletion orange, control grey, duplication blue (dose order).
# Usage: Rscript 18_figures_cn.R <aim1_dir> <out_dir>
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(ggplot2) })
args <- commandArgs(TRUE); A1 <- args[1]; OUT <- args[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CM <- file.path(A1, "cn_models")
P <- read.csv(file.path(CM, "nights_paired.csv")); W <- read.csv(file.path(CM, "nights_wrist_all.csv"))
W <- W[W$valid %in% c(TRUE, "True"), ]
TT <- read.csv(file.path(CM, "results_all_tests.csv"))
COL <- c(Deletion = "#eb6834", Control = "#8a8a86", Duplication = "#2a78d6")
prep <- function(d) {
  d <- d[d$link_status == "linked" & d$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"), ]
  d$Age <- as.numeric(d$age_final); d <- d[!is.na(d$Age), ]; d$Age_c <- d$Age - mean(d$Age)
  d$Sex <- factor(d$sex_final); d$Site <- factor(d$redcap_data_access_group); d$GUID <- d$guid
  d <- d[!is.na(d$Sex) & !is.na(d$Site), ]
  d$Genotype <- relevel(factor(d$genotype), "Control"); d
}
adj_plot <- function(d, y, ylab, title, file, back = identity, logy = FALSE, unit = "", set = "paired") {
  d <- prep(d); d <- d[!is.na(d[[y]]), ]
  m <- lmer(as.formula(paste(y, "~ Genotype + Age_c + I(Age_c^2) + Sex + Site + (1 | GUID)")), data = d, REML = TRUE)
  X <- getME(m, "X"); b <- fixef(m); cc <- setdiff(colnames(X), c("(Intercept)", grep("^Genotype", colnames(X), value = TRUE)))
  Xc <- sweep(X[, cc, drop = FALSE], 2, colMeans(X[, cc, drop = FALSE]))
  mf <- model.frame(m); mf$adj <- mf[[y]] - as.numeric(Xc %*% b[cc]); mf$genotype <- as.character(mf$Genotype)
  pp <- aggregate(adj ~ GUID + genotype, data = mf, FUN = mean)
  # model means at sample-mean covariates
  g <- c("Control", "16pDel", "16pDup", "22qDel", "22qDup"); V <- vcov(m)
  xm <- colMeans(X); mm <- do.call(rbind, lapply(g, function(k) {
    x <- xm; x[grep("^Genotype", names(x))] <- 0; if (k != "Control") x[paste0("Genotype", k)] <- 1
    data.frame(genotype = k, est = sum(x * b), se = sqrt(as.numeric(t(x) %*% as.matrix(V) %*% x))) }))
  panel <- function(df) {
    df$locus <- ifelse(df$genotype == "Control", NA, substr(df$genotype, 1, 3))
    rbind(transform(df[df$genotype %in% c("16pDel", "16pDup", "Control"), ], locus = "16p11.2"),
          transform(df[df$genotype %in% c("22qDel", "22qDup", "Control"), ], locus = "22q11.2"))
  }
  lab <- function(df) { df$copy <- ifelse(grepl("Del", df$genotype), 1, ifelse(grepl("Dup", df$genotype), 3, 2))
    df$group <- factor(ifelse(df$copy == 1, "Deletion", ifelse(df$copy == 3, "Duplication", "Control")), names(COL)); df }
  pp <- lab(panel(pp)); mm <- lab(panel(mm))
  ann <- do.call(rbind, lapply(c("16p", "22q"), function(l) {
    r <- TT[TT$outcome == y & TT$night_set == set & TT$cohort == l & TT$model == "mirror_single" & TT$random == "(1 | GUID)", ]
    if (!nrow(r)) return(NULL)
    r <- r[1, ]
    est <- if (logy) sprintf("%+.0f%% per copy", 100 * (exp(r$estimate) - 1)) else sprintf("%+.1f%s per copy", r$estimate, unit)
    data.frame(locus = paste0(l, ifelse(l == "16p", "11.2", "11.2")), txt = sprintf("mirror: %s, p = %s", est, signif(r$p_value, 2)))
  }))
  pp$yy <- back(pp$adj); mm$yy <- back(mm$est); mm$lo <- back(mm$est - 1.96 * mm$se); mm$hi <- back(mm$est + 1.96 * mm$se)
  set.seed(1)
  pl <- ggplot(pp, aes(copy, yy, colour = group)) +
    geom_jitter(width = 0.12, height = 0, size = 1.9, alpha = 0.75) +
    geom_line(data = mm, aes(copy, yy, group = 1), colour = "#444444", linewidth = 0.6, inherit.aes = FALSE) +
    geom_errorbar(data = mm, aes(copy, ymin = lo, ymax = hi), width = 0.12, colour = "#222222", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(data = mm, aes(copy, yy), shape = 23, size = 3.2, fill = "white", colour = "#222222", inherit.aes = FALSE) +
    facet_wrap(~locus) + scale_colour_manual(values = COL, name = NULL) +
    scale_x_continuous(breaks = 1:3, labels = c("Deletion\n(1 copy)", "Control\n(2 copies)", "Duplication\n(3 copies)"), limits = c(0.6, 3.4)) +
    labs(x = NULL, y = ylab, title = title, subtitle = "Dot = one person (mean of their nights), adjusted for age, sex, site. Diamond = model mean, 95% CI") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      legend.position = "none", plot.subtitle = element_text(size = 8.5, colour = "#555555"), strip.text = element_text(face = "bold", size = 11))
  if (!is.null(ann)) pl <- pl + geom_text(data = ann, aes(x = 2, y = Inf, label = txt), vjust = 1.6, size = 3.3, colour = "#222222", inherit.aes = FALSE)
  if (logy) pl <- pl + scale_y_continuous(trans = "log1p", breaks = c(0, 5, 10, 20, 40, 80, 160))
  ggsave(file.path(OUT, file), pl, width = 7.5, height = 4.2, dpi = 220, bg = "white")
  write.csv(mm, file.path(OUT, sub(".png", "_model_means.csv", file)), row.names = FALSE)
  invisible(mm)
}
adj_plot(P, "dreem_tst", "Total sleep time (minutes)", "Sleep duration, EEG headband (same nights as wrist)", "cn_tst_dreem.png", unit = " min")
adj_plot(W, "wrist_tst", "Total sleep time (minutes)", "Sleep duration, wrist device, all valid nights", "cn_tst_wrist_all.png", unit = " min", set = "all wrist")
adj_plot(P, "sol_lightsout_log", "Time to fall asleep (minutes, log scale)", "Time to fall asleep after lights out (EEG headband + mat)",
         "cn_sol_lightsout.png", back = function(x) exp(x) - 1, logy = TRUE)

# forest plot of mirror estimates
F <- TT[TT$model == "mirror_single" & TT$random == "(1 | GUID)", ]
F$cohort <- paste0(F$cohort, "11.2")
F$lo <- F$estimate - 1.96 * F$se; F$hi <- F$estimate + 1.96 * F$se
tst <- F[grepl("TST", F$outcome_label), ]; sol <- F[grepl("SOL", F$outcome_label) & grepl("log", F$outcome_label), ]
sol[, c("estimate", "lo", "hi")] <- 100 * (exp(sol[, c("estimate", "lo", "hi")]) - 1)
DEVLAB <- c("Dreem TST (min)" = "Headband (EEG), same nights as wrist", "Wrist TST (min), same nights" = "Wrist band, same nights as headband",
            "Wrist TST (min), all valid nights" = "Wrist band, all valid nights",
            "log(1+SOL), from lights-out [PRIMARY]" = "Headband + mat: from lights-out (primary)",
            "log(1+SOL), from headband start" = "Headband: from headband on",
            "log(1+SOL), from mat in-bed" = "Headband + mat: from mat in-bed",
            "log(1+SOL), wrist onset from mat in-bed, all nights" = "Wrist band + mat: from mat in-bed")
fp <- function(df, xl, file, title) {
  df$outcome_label <- ifelse(df$outcome_label %in% names(DEVLAB), DEVLAB[df$outcome_label], df$outcome_label)
  df$label <- factor(paste0(df$outcome_label), levels = rev(unique(df$outcome_label)))
  pl <- ggplot(df, aes(estimate, label)) + geom_vline(xintercept = 0, colour = "#999999") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "#444444") + geom_point(size = 2.6, colour = "#2a78d6") +
    facet_wrap(~cohort, ncol = 2) + labs(x = xl, y = NULL, title = title) + theme_minimal(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT, file), pl, width = 9, height = 0.55 * nrow(df) / 2 + 1.8, dpi = 220, bg = "white")
}
fp(tst, "Change per extra copy (minutes of sleep), 95% CI", "forest_tst.png", "Sleep duration: effect of each extra copy")
fp(sol, "Change per extra copy (% time to fall asleep), 95% CI", "forest_sol.png", "Time to fall asleep: effect of each extra copy")
cat("figures written to", OUT, "\n")
