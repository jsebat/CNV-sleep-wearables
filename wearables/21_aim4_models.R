# Aim 4, step 2: copy-number models for the wrist-band phenotypes (20_phenotypes.py).
# Same specification as the Aim 1 / sleep-architecture models, with age at recording and a family random effect:
#   night/day level : y ~ term + Age_c + I(Age_c^2) + Sex + Site + (1 | GUID) + (1 | family)
#   participant level: y ~ term + Age_c + I(Age_c^2) + Sex + Site + (1 | family)
#   per locus: mirror copy_c (-1/0/+1), U-shape copy_q (1/0/1), joint (LRT, ML); five-group vs control.
# FDR (BH) within each domain over the primary mirror + U-shape single-term tests.
# Env: AGE_MIN / AGE_MAX (age at recording) and OUT_TAG.
#      SEASON=1 adds season as sin + cos of day of year (night/day rows: that night's/day's date; participant rows:
#      the participant's median recording date). All four sites are northern-hemisphere.
# Usage: Rscript 21_aim4_models.R <analysis_root>   (reads phase0/ and aim4/)
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
args <- commandArgs(TRUE); ROOT <- args[1]; A4 <- file.path(ROOT, "aim4")
AGE_MIN <- as.numeric(Sys.getenv("AGE_MIN", "-Inf")); AGE_MAX <- as.numeric(Sys.getenv("AGE_MAX", "Inf")); OUT_TAG <- Sys.getenv("OUT_TAG", "")
L <- read.csv(file.path(ROOT, "phase0", "linkage.csv")); COV <- read.csv(file.path(ROOT, "phase0", "covariates.csv"))
base <- merge(L[L$link_status == "linked" & L$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"),
                c("empatica_id", "guid", "genotype", "sex_final", "redcap_data_access_group", "rarecnv_family_id")],
              COV[, c("empatica_id", "age_at_recording")], by = "empatica_id")
base <- base[!is.na(base$age_at_recording) & base$age_at_recording >= AGE_MIN & base$age_at_recording <= AGE_MAX, ]
base$Age_c <- base$age_at_recording - mean(base$age_at_recording); base$Sex <- factor(base$sex_final)
base$Site <- factor(base$redcap_data_access_group); base$GUID <- base$guid
base$family <- ifelse(is.na(base$rarecnv_family_id) | base$rarecnv_family_id == "", base$guid, base$rarecnv_family_id)
base$Genotype <- relevel(factor(base$genotype), "Control")
PH <- list(participant = read.csv(file.path(A4, "phenotypes_participant.csv")), night = read.csv(file.path(A4, "phenotypes_night.csv")),
           day = read.csv(file.path(A4, "phenotypes_day.csv")))
PH$participant <- PH$participant[PH$participant$n_valid_days >= 3, ]
SEASON <- Sys.getenv("SEASON", "0") == "1"
doy <- function(x) as.numeric(format(as.Date(x), "%j"))
PH$night$doy <- doy(PH$night$night); PH$day$doy <- doy(PH$day$date)
PH$participant <- merge(PH$participant, aggregate(list(doy = as.numeric(as.Date(PH$day$date))), list(empatica_id = PH$day$empatica_id),
                        function(z) doy(as.Date(median(z), origin = "1970-01-01"))), by = "empatica_id", all.x = TRUE)
for (k in names(PH)) { PH[[k]]$sin_doy <- sin(2 * pi * PH[[k]]$doy / 365.25); PH[[k]]$cos_doy <- cos(2 * pi * PH[[k]]$doy / 365.25) }
DICT <- read.csv(text = "
var,level,domain,label
tst_h,night,Sleep duration & quality,Wrist band: total sleep time (h)
spt_h,night,Sleep duration & quality,Wrist band: sleep period length (h)
sleep_eff,night,Sleep duration & quality,Wrist band: sleep efficiency
waso_min,night,Sleep duration & quality,Wrist band: wake after sleep onset (min)
n_wake_bouts,night,Sleep duration & quality,Wrist band: night-time wake bouts (n)
onset_clock_h,night,Sleep timing,Wrist band: sleep onset clock time (h after 18:00)
offset_clock_h,night,Sleep timing,Wrist band: wake-up clock time (h after 18:00)
midpoint_clock_h,night,Sleep timing,Wrist band: sleep midpoint (h after 18:00)
SRI,participant,Sleep regularity,Wrist band: Sleep Regularity Index
tst_sd_h,participant,Sleep regularity,Wrist band: night-to-night SD of sleep time (h)
midpoint_sd_h,participant,Sleep regularity,Wrist band: night-to-night SD of sleep midpoint (h)
social_jetlag_h,participant,Sleep regularity,Wrist band: social jet lag (h)
IS,participant,Rest-activity rhythm,Wrist band: interdaily stability
IV,participant,Rest-activity rhythm,Wrist band: intradaily variability
RA,participant,Rest-activity rhythm,Wrist band: relative amplitude
cos_amp,participant,Rest-activity rhythm,Wrist band: cosinor amplitude (log10 mg)
cos_acrophase_h,participant,Rest-activity rhythm,Wrist band: cosinor peak time (clock h)
m10_onset_h,day,Rest-activity rhythm,Wrist band: most-active 10 h onset (clock h)
l5_onset_h,day,Rest-activity rhythm,Wrist band: least-active 5 h onset (clock h)
enmo_mg,day,Physical activity,Wrist band: mean acceleration ENMO (mg)
m10_mg,day,Physical activity,Wrist band: most-active 10 h level (mg)
l5_mg,day,Physical activity,Wrist band: least-active 5 h level (mg)
sed_h,day,Physical activity,Wrist band: sedentary time (h/day; cut-point)
light_h,day,Physical activity,Wrist band: light activity (h/day; cut-point)
mvpa_min,day,Physical activity,Wrist band: moderate-vigorous activity (min/day; cut-point)
hr_sleep_min,night,Physiology (exploratory),Wrist band: lowest heart rate in sleep (bpm)
hr_sleep_mean,night,Physiology (exploratory),Wrist band: mean heart rate in sleep (bpm)
rmssd_sleep,night,Physiology (exploratory),Wrist band: HRV in sleep RMSSD (ms)
temp_amp,day,Physiology (exploratory),Wrist band: wrist temperature range (C)
")
fit <- function(f, d, REML = TRUE) tryCatch(suppressMessages(lmer(f, data = d, REML = REML, control = lmerControl(check.conv.singular = "ignore"))), error = function(e) NULL)
cov_terms <- function(d) { k <- c("Age_c", "I(Age_c^2)"); if (nlevels(droplevels(d$Sex)) > 1) k <- c(k, "Sex"); if (nlevels(droplevels(d$Site)) > 1) k <- c(k, "Site"); if (SEASON) k <- c(k, "sin_doy", "cos_doy"); k }
res <- list(); desc <- list()
for (i in seq_len(nrow(DICT))) {
  v <- DICT$var[i]; lev <- DICT$level[i]
  d0 <- merge(PH[[lev]][, c("empatica_id", v, "sin_doy", "cos_doy")], base, by = "empatica_id"); d0 <- d0[is.finite(d0[[v]]), ]
  if (SEASON) d0 <- d0[is.finite(d0$sin_doy), ]
  if (nrow(d0) < 30) next
  re <- if (lev == "participant") "(1 | family)" else "(1 | GUID) + (1 | family)"
  fml <- function(g, d) as.formula(paste(v, "~", paste(c(g, cov_terms(d), re), collapse = " + ")))
  desc[[length(desc) + 1]] <- do.call(rbind, lapply(split(d0, d0$Genotype), function(x)
    data.frame(var = v, label = DICT$label[i], genotype = x$Genotype[1], people = length(unique(x$GUID)), obs = nrow(x), mean = mean(x[[v]]), sd = sd(x[[v]]))))
  out <- function(m, term, cohort, model, d) { if (is.null(m)) return(NULL); co <- summary(m)$coefficients; r <- grep(term, rownames(co), fixed = TRUE)
    if (!length(r)) return(NULL)
    data.frame(var = v, label = DICT$label[i], domain = DICT$domain[i], level = lev, cohort = cohort, model = model, term = rownames(co)[r],
               estimate = co[r, 1], se = co[r, 2], p_value = co[r, 5], sd_outcome = sd(d[[v]]), n_obs = nrow(model.frame(m)),
               n_people = length(unique(d$GUID)), singular = isSingular(m)) }
  m5 <- fit(fml("Genotype", d0), d0); res[[length(res) + 1]] <- out(m5, "Genotype", "All5", "group_vs_control", d0)
  for (locus in c("16p", "22q")) {
    d <- d0[d0$Genotype %in% c(paste0(locus, "Del"), paste0(locus, "Dup"), "Control"), ]
    d$copy_c <- ifelse(grepl("Del", d$Genotype), -1, ifelse(grepl("Dup", d$Genotype), 1, 0)); d$copy_q <- abs(d$copy_c)
    res[[length(res) + 1]] <- out(fit(fml("copy_c", d), d), "copy_c", locus, "mirror_single", d)
    res[[length(res) + 1]] <- out(fit(fml("copy_q", d), d), "copy_q", locus, "ushape_single", d)
    fj <- fit(fml(c("copy_c", "copy_q"), d), d, FALSE); fm <- fit(fml("copy_c", d), d, FALSE); fu <- fit(fml("copy_q", d), d, FALSE)
    if (!is.null(fj) && !is.null(fm) && !is.null(fu)) {
      a1 <- suppressMessages(anova(fm, fj)); a2 <- suppressMessages(anova(fu, fj))
      res[[length(res) + 1]] <- data.frame(var = v, label = DICT$label[i], domain = DICT$domain[i], level = lev, cohort = locus, model = "joint",
        term = c("U-shape beyond mirror", "mirror beyond U-shape"), estimate = NA, se = NA, p_value = c(a1$`Pr(>Chisq)`[2], a2$`Pr(>Chisq)`[2]),
        sd_outcome = sd(d[[v]]), n_obs = nrow(d), n_people = length(unique(d$GUID)), singular = isSingular(fj))
    }
  }
}
T <- do.call(rbind, res); D <- do.call(rbind, desc)
T$std_estimate <- T$estimate / T$sd_outcome
T$q_domain <- NA
for (dm in unique(T$domain)) { k <- T$domain == dm & T$model %in% c("mirror_single", "ushape_single"); T$q_domain[k] <- p.adjust(T$p_value[k], "BH") }
write.csv(T, file.path(A4, paste0("aim4_results", OUT_TAG, ".csv")), row.names = FALSE)
write.csv(D, file.path(A4, paste0("aim4_descriptives", OUT_TAG, ".csv")), row.names = FALSE)
options(width = 230)
S <- T[T$model %in% c("mirror_single", "ushape_single"), c("label", "cohort", "model", "estimate", "std_estimate", "p_value", "q_domain", "n_people")]
S <- S[order(S$p_value), ]; S$estimate <- signif(S$estimate, 3); S$std_estimate <- round(S$std_estimate, 2); S$p_value <- signif(S$p_value, 2); S$q_domain <- signif(S$q_domain, 2)
cat("AGE range", AGE_MIN, AGE_MAX, "| people in base", nrow(base), "\n"); print(head(S, 30), row.names = FALSE)
cat("\ntests with q_domain < 0.05:", sum(S$q_domain < 0.05, na.rm = TRUE), "of", nrow(S), "\n")
