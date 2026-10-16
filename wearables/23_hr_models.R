# Aim 4, heart-rate confounder checks - part 2: refit the sleep heart-rate effects (WRIST BAND) under each check.
# Model as 21_aim4_models.R (night level): y ~ term + Age_c + I(Age_c^2) + Sex + Site [+ extra] + (1 | GUID) + (1 | family)
# Reported per variant: five-group vs control, 16p and 22q mirror and U-shape terms.
# Variants
#   base        mean HR in sleep (as Aim 4)
#   still       HR from still sleep epochs only (hr_still_mean)
#   no_meds     carriers on heart-rate-active psychotropics removed (stimulant, non-stimulant ADHD, alpha-2 agonist,
#               antipsychotic, antidepressant, lithium); carriers with no medication form also removed; controls kept
#               (family members were mostly not asked - their medication status is unknown)
#   bmi         + BMI (enrolment) as covariate
#   no_chd      participants with a congenital heart malformation (med_birthdefects_3 = 1) removed
#   all_checks  no_meds + no_chd + BMI, on still-epoch HR
# Signal quality (hr_coverage, enmo_sleep_mg) is modelled as an outcome by genotype.
# Within carriers: HR ~ heart-rate-active medication + Genotype + covariates.
# Usage: Rscript 23_hr_models.R <analysis_root>
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
ROOT <- commandArgs(TRUE)[1]; A4 <- file.path(ROOT, "aim4")
L <- read.csv(file.path(ROOT, "phase0", "linkage.csv")); COV <- read.csv(file.path(ROOT, "phase0", "covariates.csv"))
HH <- read.csv(file.path(A4, "heart_history.csv"))
N <- merge(read.csv(file.path(A4, "phenotypes_night.csv"))[, c("empatica_id", "night", "hr_sleep_mean", "rmssd_sleep")],
           read.csv(file.path(A4, "hr_checks_night.csv")), by = c("empatica_id", "night"))
b <- L[L$link_status == "linked" & L$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"),
       c("empatica_id", "guid", "genotype", "sex_final", "redcap_data_access_group", "rarecnv_family_id")]
b <- merge(b, COV[, c("empatica_id", "age_at_recording", "bmi", "meds_form_present", "med_stimulant", "med_adhd_nonstim", "med_alpha2",
                      "med_antipsychotic", "med_antidepressant", "med_lithium")], by = "empatica_id")
b <- merge(b, HH[, c("empatica_id", "med_birthdefects_3")], by = "empatica_id", all.x = TRUE)
b$hr_med <- as.integer(rowSums(b[, c("med_stimulant", "med_adhd_nonstim", "med_alpha2", "med_antipsychotic", "med_antidepressant", "med_lithium")], na.rm = TRUE) > 0)
b$carrier <- b$genotype != "Control"; b$chd <- b$med_birthdefects_3 %in% c(1, "1")
b <- b[!is.na(b$age_at_recording), ]
b$Age_c <- b$age_at_recording - mean(b$age_at_recording); b$Sex <- factor(b$sex_final); b$Site <- factor(b$redcap_data_access_group)
b$GUID <- b$guid; b$family <- ifelse(is.na(b$rarecnv_family_id) | b$rarecnv_family_id == "", b$guid, b$rarecnv_family_id)
b$Genotype <- relevel(factor(b$genotype), "Control")
D <- merge(N, b, by = "empatica_id")
cat("carriers on HR-active meds:", sum(b$carrier & b$hr_med == 1), "of", sum(b$carrier), "| carriers without meds form:", sum(b$carrier & b$meds_form_present %in% 0),
    "| CHD:", paste(names(table(b$genotype[b$chd])), table(b$genotype[b$chd]), collapse = ", "), "\n")
fit <- function(f, d) tryCatch(suppressMessages(lmer(f, data = d, control = lmerControl(check.conv.singular = "ignore"))), error = function(e) NULL)
cvt <- function(d, extra) { k <- c("Age_c", "I(Age_c^2)", extra); if (nlevels(droplevels(d$Sex)) > 1) k <- c(k, "Sex"); if (nlevels(droplevels(d$Site)) > 1) k <- c(k, "Site"); k }
run <- function(d, y, variant, extra = character(0)) {
  d <- d[is.finite(d[[y]]) & complete.cases(d[, extra, drop = FALSE]), ]
  f <- function(g, dd) as.formula(paste(y, "~", paste(c(g, cvt(dd, extra), "(1 | GUID) + (1 | family)"), collapse = " + ")))
  pick <- function(m, term, cohort) { if (is.null(m)) return(NULL); co <- summary(m)$coefficients; r <- grep(term, rownames(co), fixed = TRUE)
    data.frame(variant = variant, outcome = y, cohort = cohort, term = rownames(co)[r], estimate = co[r, 1], se = co[r, 2], p = co[r, 5],
               people = length(unique(d$GUID[d$Genotype %in% unique(model.frame(m)$Genotype %||% d$Genotype)])), nights = nrow(model.frame(m))) }
  out <- list(pick(fit(f("Genotype", d), d), "Genotype", "All5"))
  for (l in c("16p", "22q")) { dd <- d[d$Genotype %in% c(paste0(l, "Del"), paste0(l, "Dup"), "Control"), ]
    dd$copy_c <- ifelse(grepl("Del", dd$Genotype), -1, ifelse(grepl("Dup", dd$Genotype), 1, 0)); dd$copy_q <- abs(dd$copy_c)
    out[[length(out) + 1]] <- pick(fit(f("copy_c", dd), dd), "copy_c", l); out[[length(out) + 1]] <- pick(fit(f("copy_q", dd), dd), "copy_q", l) }
  do.call(rbind, out)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
nomed <- D[!(D$carrier & (D$hr_med == 1 | D$meds_form_present %in% 0)), ]
R <- rbind(run(D, "hr_sleep_mean", "base"), run(D, "hr_still_mean", "still"), run(nomed, "hr_sleep_mean", "no_meds"),
           run(D, "hr_sleep_mean", "bmi", "bmi"), run(D[!D$chd, ], "hr_sleep_mean", "no_chd"),
           run(nomed[!nomed$chd, ], "hr_still_mean", "all_checks", "bmi"),
           run(D, "rmssd_sleep", "base"), run(nomed[!nomed$chd, ], "rmssd_still", "all_checks", "bmi"),
           run(D, "hr_coverage", "signal: share of sleep epochs with HR"), run(D, "enmo_sleep_mg", "signal: movement in sleep (mg)"))
write.csv(R, file.path(A4, "hr_checks_results.csv"), row.names = FALSE)
# within carriers: does HR-active medication relate to sleep HR?
C <- D[D$carrier & D$meds_form_present %in% 1, ]
mw <- fit(as.formula(paste("hr_sleep_mean ~ hr_med + Genotype +", paste(cvt(C, character(0)), collapse = " + "), "+ (1 | GUID) + (1 | family)")), C)
cw <- summary(mw)$coefficients["hr_med", ]
options(width = 200)
show <- R[R$cohort != "All5" | grepl("signal", R$variant), ]; show$estimate <- round(show$estimate, 2); show$p <- signif(show$p, 2); print(show[, c("variant", "outcome", "cohort", "term", "estimate", "p", "nights")], row.names = FALSE)
g <- R[R$cohort == "All5" & !grepl("signal", R$variant), ]; g$estimate <- round(g$estimate, 1); g$p <- signif(g$p, 2)
print(reshape(g[, c("variant", "term", "estimate")], idvar = "variant", timevar = "term", direction = "wide"), row.names = FALSE)
cat(sprintf("\nwithin carriers with a meds form (%d people): HR-active medication effect on sleep HR = %+.2f bpm (SE %.2f, p = %.2g)\n",
            length(unique(C$GUID)), cw[1], cw[2], cw[5]))
