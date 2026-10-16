# Aim 4 report figures (WRIST BAND only).
#  aim4_overview.png     standardized mirror / U-shape estimates (per SD of the phenotype) for every phenotype, by locus;
#                        filled = q_domain < 0.05
#  aim4_cn_<var>.png     copy-number plots: partial residuals from one reference model per phenotype on the five-group
#                        cohort (genotype kept, covariates centred), participant means, model means +/- 95 % CI
#  aim4_hr_checks.png    heart rate in sleep, each carrier group vs controls, under each confounder check (95 % CI)
# Env SEASON=1 (primary since 2026-10-14): read aim4_results_season.csv and add sin/cos day-of-year to the plot models.
# Usage: Rscript 24_figures_aim4.R <analysis_root> <out_dir>
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(ggplot2) })
a <- commandArgs(TRUE); ROOT <- a[1]; OUT <- a[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); A4 <- file.path(ROOT, "aim4")
COL <- c(Deletion = "#eb6834", Control = "#8a8a86", Duplication = "#2a78d6")
SEASON <- Sys.getenv("SEASON", "0") == "1"
TT <- read.csv(file.path(A4, if (SEASON) "aim4_results_season.csv" else "aim4_results.csv"))
SUBT <- paste0("Age at recording, age², sex, site", if (SEASON) ", season" else "", "; person and family random effects; FDR within each domain")
th <- theme_minimal(base_size = 10.5) + theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 8.5, colour = "#555555"))

# ---- overview ---------------------------------------------------------------------------------------
O <- TT[TT$model %in% c("mirror_single", "ushape_single"), ]
O$lo <- (O$estimate - 1.96 * O$se) / O$sd_outcome; O$hi <- (O$estimate + 1.96 * O$se) / O$sd_outcome
O$model <- factor(ifelse(O$model == "mirror_single", "Mirror (per copy)", "U-shape (both CNVs vs control)"), c("Mirror (per copy)", "U-shape (both CNVs vs control)"))
O$locus <- paste0(O$cohort, "11.2"); O$sig <- ifelse(O$q_domain < 0.05, "q < 0.05 (within domain)", "not significant")
O$lab <- sub("^Wrist band: ", "", O$label)
lv <- unique(O[order(O$domain, O$lab), "lab"]); O$lab <- factor(O$lab, rev(lv))
pl <- ggplot(O, aes(std_estimate, lab, colour = locus, shape = sig)) + geom_vline(xintercept = 0, colour = "#999999") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, alpha = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.6), fill = "white") +
  scale_shape_manual(values = c("q < 0.05 (within domain)" = 16, "not significant" = 21), name = NULL) +
  scale_colour_manual(values = c("16p11.2" = "#4a3aa7", "22q11.2" = "#1baf7a"), name = NULL) + facet_wrap(~model) +
  coord_cartesian(xlim = c(-1.6, 1.6)) + labs(x = "Effect in SD units of the phenotype, 95% CI", y = NULL,
  title = "Wrist band: copy-number effects on every phenotype", subtitle = SUBT) +
  th + theme(legend.position = "bottom", axis.text.y = element_text(size = 8.5))
ggsave(file.path(OUT, "aim4_overview.png"), pl, width = 10, height = 8.6, dpi = 200, bg = "white")

# ---- copy-number plots -----------------------------------------------------------------------------
L <- read.csv(file.path(ROOT, "phase0", "linkage.csv")); COV <- read.csv(file.path(ROOT, "phase0", "covariates.csv"))
base <- merge(L[L$link_status == "linked" & L$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"),
                c("empatica_id", "guid", "genotype", "sex_final", "redcap_data_access_group", "rarecnv_family_id")], COV[, c("empatica_id", "age_at_recording")])
base <- base[!is.na(base$age_at_recording), ]; base$Age_c <- base$age_at_recording - mean(base$age_at_recording)
base$Sex <- factor(base$sex_final); base$Site <- factor(base$redcap_data_access_group); base$GUID <- base$guid
base$family <- ifelse(is.na(base$rarecnv_family_id) | base$rarecnv_family_id == "", base$guid, base$rarecnv_family_id)
base$Genotype <- relevel(factor(base$genotype), "Control")
PH <- list(night = read.csv(file.path(A4, "phenotypes_night.csv")), participant = read.csv(file.path(A4, "phenotypes_participant.csv")))
PH$participant <- PH$participant[PH$participant$n_valid_days >= 3, ]
doy <- function(x) as.numeric(format(as.Date(x), "%j")); DAY <- read.csv(file.path(A4, "phenotypes_day.csv"))
PH$night$doy <- doy(PH$night$night)
PH$participant <- merge(PH$participant, aggregate(list(doy = as.numeric(as.Date(DAY$date))), list(empatica_id = DAY$empatica_id),
                        function(z) doy(as.Date(median(z), origin = "1970-01-01"))), by = "empatica_id", all.x = TRUE)
for (k in names(PH)) { PH[[k]]$sin_doy <- sin(2 * pi * PH[[k]]$doy / 365.25); PH[[k]]$cos_doy <- cos(2 * pi * PH[[k]]$doy / 365.25) }
cnplot <- function(v, lev, ylab, title, unit) {
  d <- merge(PH[[lev]][, c("empatica_id", v, "sin_doy", "cos_doy")], base, by = "empatica_id"); d <- d[is.finite(d[[v]]), ]
  if (SEASON) d <- d[is.finite(d$sin_doy), ]
  re <- if (lev == "participant") "(1 | family)" else "(1 | GUID) + (1 | family)"
  m <- lmer(as.formula(paste(v, "~ Genotype + Age_c + I(Age_c^2) + Sex + Site +", if (SEASON) "sin_doy + cos_doy +", re)), data = d)
  X <- getME(m, "X"); b <- fixef(m); cc <- setdiff(colnames(X), c("(Intercept)", grep("^Genotype", colnames(X), value = TRUE)))
  Xc <- sweep(X[, cc, drop = FALSE], 2, colMeans(X[, cc, drop = FALSE])); mf <- model.frame(m)
  mf$adj <- mf[[v]] - as.numeric(Xc %*% b[cc]); mf$GUID <- d$GUID[as.numeric(rownames(mf))] %||% mf$GUID
  pp <- aggregate(adj ~ GUID + Genotype, data = cbind(mf, GUID = if ("GUID" %in% names(mf)) mf$GUID else d$GUID), FUN = mean)
  xm <- colMeans(X); V <- as.matrix(vcov(m))
  mm <- do.call(rbind, lapply(c("Control", "16pDel", "16pDup", "22qDel", "22qDup"), function(k) { x <- xm; x[grep("^Genotype", names(x))] <- 0
    if (k != "Control") x[paste0("Genotype", k)] <- 1; data.frame(Genotype = k, est = sum(x * b), se = sqrt(as.numeric(t(x) %*% V %*% x))) }))
  pan <- function(df) rbind(transform(df[df$Genotype %in% c("16pDel", "16pDup", "Control"), ], locus = "16p11.2"),
                            transform(df[df$Genotype %in% c("22qDel", "22qDup", "Control"), ], locus = "22q11.2"))
  lab <- function(df) { df$copy <- ifelse(grepl("Del", df$Genotype), 1, ifelse(grepl("Dup", df$Genotype), 3, 2))
    df$group <- factor(ifelse(df$copy == 1, "Deletion", ifelse(df$copy == 3, "Duplication", "Control")), names(COL)); df }
  pp <- lab(pan(pp)); mm <- lab(pan(mm))
  ann <- do.call(rbind, lapply(c("16p", "22q"), function(l) { r <- TT[TT$var == v & TT$cohort == l & TT$model %in% c("mirror_single", "ushape_single"), ]
    data.frame(locus = paste0(l, "11.2"), txt = sprintf("mirror %+.2g%s/copy (p=%.1g); U %+.2g%s (p=%.1g)", r$estimate[r$model == "mirror_single"], unit,
               r$p_value[r$model == "mirror_single"], r$estimate[r$model == "ushape_single"], unit, r$p_value[r$model == "ushape_single"])) }))
  set.seed(1)
  pl <- ggplot(pp, aes(copy, adj, colour = group)) + geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
    geom_line(data = mm, aes(copy, est, group = 1), colour = "#444444", linewidth = 0.6, inherit.aes = FALSE) +
    geom_errorbar(data = mm, aes(copy, ymin = est - 1.96 * se, ymax = est + 1.96 * se), width = 0.12, colour = "#222222", inherit.aes = FALSE) +
    geom_point(data = mm, aes(copy, est), shape = 23, size = 3, fill = "white", colour = "#222222", inherit.aes = FALSE) +
    geom_text(data = ann, aes(x = 2, y = Inf, label = txt), vjust = 1.5, size = 2.45, colour = "#222222", inherit.aes = FALSE) +
    facet_wrap(~locus) + scale_colour_manual(values = COL) +
    scale_x_continuous(breaks = 1:3, labels = c("Deletion\n(1 copy)", "Control\n(2 copies)", "Duplication\n(3 copies)"), limits = c(0.6, 3.4)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = ylab, title = title, subtitle = paste0("Wrist band. Dot = one person, adjusted for age, sex, site", if (SEASON) ", season", ". Diamond = model mean, 95% CI")) + th + theme(legend.position = "none")
  ggsave(file.path(OUT, paste0("aim4_cn_", v, ".png")), pl, width = 7.5, height = 4.2, dpi = 220, bg = "white")
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
cnplot("hr_sleep_mean", "night", "Heart rate in sleep (bpm)", "Heart rate during sleep", " bpm")
cnplot("tst_h", "night", "Total sleep time (h)", "Total sleep time", " h")
cnplot("IV", "participant", "Intradaily variability", "Rest-activity fragmentation (intradaily variability)", "")

# ---- HR confounder checks --------------------------------------------------------------------------
H <- read.csv(file.path(A4, "hr_checks_results.csv"))
H <- H[H$cohort == "All5" & H$outcome %in% c("hr_sleep_mean", "hr_still_mean"), ]
H$check <- factor(c(base = "Base result", still = "Still sleep epochs only", no_meds = "Medicated carriers removed", bmi = "+ BMI covariate",
                    no_chd = "Known heart malformation removed", all_checks = "All checks together")[H$variant],
                  rev(c("Base result", "Still sleep epochs only", "Medicated carriers removed", "+ BMI covariate", "Known heart malformation removed", "All checks together")))
H$group <- sub("Genotype", "", H$term); H$group <- factor(c("16pDel" = "16p deletion", "16pDup" = "16p duplication", "22qDel" = "22q deletion", "22qDup" = "22q duplication")[H$group],
                                                        c("16p deletion", "16p duplication", "22q deletion", "22q duplication"))
H$kind <- ifelse(grepl("del", H$group), "Deletion", "Duplication")
pl <- ggplot(H, aes(estimate, check, colour = kind)) + geom_vline(xintercept = 0, colour = "#999999") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), height = 0.2) + geom_point(size = 2.4) +
  facet_wrap(~group, nrow = 1) + scale_colour_manual(values = COL[c("Deletion", "Duplication")], guide = "none") +
  labs(x = "Heart rate in sleep vs controls (bpm), 95% CI", y = NULL, title = "Wrist band: heart rate during sleep under each confounder check") + th
ggsave(file.path(OUT, "aim4_hr_checks.png"), pl, width = 11, height = 3.6, dpi = 200, bg = "white")
cat("done\n")
