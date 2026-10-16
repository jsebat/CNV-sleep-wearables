# Copy-number models for sleep timing/duration outcomes - same specification as the sleep-architecture
# analysis (SleepArchitectureResults/code/R/03_models.R):
#   y ~ <genotype term> + Age_c + I(Age_c^2) + Sex + Site + (1 | GUID)
#   per locus (copy number 1/2/3): mirror copy_c = -1/0/+1, U-shape copy_q = 1/0/1, joint = both (LRT, ML)
#   five-group: Genotype vs Control.  Wald tests REML + Satterthwaite (lmerTest).
# Sensitivity: + (1 | family) random effect.
# Inputs : <aim1>/cn_models/nights_paired.csv, nights_wrist_all.csv (16_night_table.py)
# Output : <aim1>/cn_models/results_all_tests.csv, group_descriptives.csv
# Usage  : Rscript 17_cn_models.R <aim1_dir>
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
args <- commandArgs(TRUE); R <- file.path(args[1], "cn_models")
# options (environment): AGE_VAR = age_final (default, age at enrolment - as the sleep paper) | age_at_recording
#                        AGE_MIN / AGE_MAX = restrict to an age range (on AGE_VAR);  OUT_TAG = suffix for output files
AGE_VAR <- Sys.getenv("AGE_VAR", "age_final"); AGE_MIN <- as.numeric(Sys.getenv("AGE_MIN", "-Inf")); AGE_MAX <- as.numeric(Sys.getenv("AGE_MAX", "Inf"))
OUT_TAG <- Sys.getenv("OUT_TAG", "")
COV <- read.csv(file.path(dirname(args[1]), "phase0", "covariates.csv"))[, c("empatica_id", "age_at_recording")]
P <- read.csv(file.path(R, "nights_paired.csv")); A <- read.csv(file.path(R, "nights_wrist_all.csv"))
A <- A[A$valid %in% c(TRUE, "True"), ]
prep <- function(d) {
  d <- d[d$link_status == "linked" & d$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"), ]
  d <- merge(d, COV, by = "empatica_id", all.x = TRUE)
  d$Age <- as.numeric(d[[AGE_VAR]]); d <- d[!is.na(d$Age) & d$Age >= AGE_MIN & d$Age <= AGE_MAX, ]
  d$Age_c <- d$Age - mean(d$Age); d$Sex <- factor(d$sex_final); d$Site <- factor(d$redcap_data_access_group)
  d$GUID <- d$guid; d$family <- ifelse(is.na(d$rarecnv_family_id) | d$rarecnv_family_id == "", d$guid, d$rarecnv_family_id)
  d$Genotype <- relevel(factor(d$genotype), "Control"); d
}
P <- prep(P); A <- prep(A)
outcomes <- list(
  list(d = "P", y = "dreem_tst",            lab = "Dreem TST (min)",                       set = "paired"),
  list(d = "P", y = "wrist_tst",            lab = "Wrist TST (min), same nights",          set = "paired"),
  list(d = "A", y = "wrist_tst",            lab = "Wrist TST (min), all valid nights",     set = "all wrist"),
  list(d = "P", y = "sol_lightsout_log",    lab = "log(1+SOL), from lights-out [PRIMARY]", set = "paired"),
  list(d = "P", y = "sol_headband_log",     lab = "log(1+SOL), from headband start",       set = "paired"),
  list(d = "P", y = "sol_matinbed_log",     lab = "log(1+SOL), from mat in-bed",           set = "paired"),
  list(d = "P", y = "sol_lightsout",        lab = "SOL (min), from lights-out",            set = "paired"),
  list(d = "A", y = "wrist_sol_lightsout_log", lab = "log(1+SOL), wrist onset from mat in-bed, all nights", set = "all wrist"),
  list(d = "P", y = "solp_lightsout_log",  lab = "log(1+latency to PERSISTENT sleep), from lights-out", set = "paired"),
  list(d = "P", y = "solp_headband_log",   lab = "log(1+latency to PERSISTENT sleep), from headband start", set = "paired"),
  list(d = "P", y = "solp_matinbed_log",   lab = "log(1+latency to PERSISTENT sleep), from mat in-bed", set = "paired"),
  list(d = "P", y = "solp_lightsout",      lab = "Latency to PERSISTENT sleep (min), from lights-out", set = "paired"))
fit <- function(f, d, REML = TRUE) tryCatch(suppressMessages(lmer(f, data = d, REML = REML,
           control = lmerControl(check.conv.singular = "ignore"))), error = function(e) NULL)
cov_terms <- function(d) { k <- c("Age_c", "I(Age_c^2)"); if (nlevels(droplevels(d$Sex)) > 1) k <- c(k, "Sex")
  if (nlevels(droplevels(d$Site)) > 1) k <- c(k, "Site"); k }
fml <- function(y, g, d, re = "(1 | GUID)") as.formula(paste(y, "~", paste(c(g, cov_terms(d), re), collapse = " + ")))
wald <- function(m, term, o, cohort, model, d, re) {
  if (is.null(m)) return(NULL); co <- summary(m)$coefficients; rows <- grep(term, rownames(co), fixed = TRUE)
  if (!length(rows)) return(NULL)
  data.frame(outcome = o$y, outcome_label = o$lab, night_set = o$set, cohort = cohort, model = model, random = re,
             term = rownames(co)[rows], estimate = co[rows, 1], se = co[rows, 2], df = co[rows, 3], p_value = co[rows, 5],
             n_nights = nrow(model.frame(m)), n_subjects = length(unique(model.frame(m)$GUID)), singular = isSingular(m))
}
lrt <- function(full, red, o, cohort, label, re) {
  if (is.null(full) || is.null(red)) return(NULL); a <- suppressMessages(anova(red, full))
  data.frame(outcome = o$y, outcome_label = o$lab, night_set = o$set, cohort = cohort, model = "joint", random = re,
             term = label, estimate = NA, se = NA, df = a$Df[2], p_value = a$`Pr(>Chisq)`[2],
             n_nights = nrow(model.frame(full)), n_subjects = length(unique(model.frame(full)$GUID)), singular = isSingular(full))
}
res <- list(); desc <- list()
for (o in outcomes) {
  d0 <- get(o$d); d0 <- d0[!is.na(d0[[o$y]]), ]
  desc[[length(desc) + 1]] <- do.call(rbind, lapply(split(d0, d0$Genotype), function(x)
    data.frame(outcome = o$y, night_set = o$set, genotype = x$Genotype[1], subjects = length(unique(x$GUID)), nights = nrow(x),
               mean = mean(x[[o$y]]), median = median(x[[o$y]]), age_mean = mean(tapply(x$Age, x$GUID, mean)))))
  for (re in c("(1 | GUID)", "(1 | GUID) + (1 | family)")) {
    m5 <- fit(fml(o$y, "Genotype", d0, re), d0)
    res[[length(res) + 1]] <- wald(m5, "Genotype", o, "All5", "group_vs_control", d0, re)
    for (locus in c("16p", "22q")) {
      d <- d0[d0$Genotype %in% c(paste0(locus, "Del"), paste0(locus, "Dup"), "Control"), ]
      d$copy_c <- ifelse(grepl("Del", d$Genotype), -1, ifelse(grepl("Dup", d$Genotype), 1, 0)); d$copy_q <- abs(d$copy_c)
      mm <- fit(fml(o$y, "copy_c", d, re), d); mu <- fit(fml(o$y, "copy_q", d, re), d)
      res[[length(res) + 1]] <- wald(mm, "copy_c", o, locus, "mirror_single", d, re)
      res[[length(res) + 1]] <- wald(mu, "copy_q", o, locus, "ushape_single", d, re)
      fj <- fit(fml(o$y, c("copy_c", "copy_q"), d, re), d, REML = FALSE)
      fm <- fit(fml(o$y, "copy_c", d, re), d, REML = FALSE); fu <- fit(fml(o$y, "copy_q", d, re), d, REML = FALSE)
      res[[length(res) + 1]] <- lrt(fj, fm, o, locus, "U-shape beyond mirror (1 df)", re)
      res[[length(res) + 1]] <- lrt(fj, fu, o, locus, "mirror beyond U-shape (1 df)", re)
    }
  }
}
T <- do.call(rbind, res); D <- do.call(rbind, desc)
prim <- T$model %in% c("mirror_single", "ushape_single") & T$random == "(1 | GUID)"
T$q_value <- NA; T$q_value[prim] <- p.adjust(T$p_value[prim], "BH")
write.csv(T, file.path(R, paste0("results_all_tests", OUT_TAG, ".csv")), row.names = FALSE); write.csv(D, file.path(R, paste0("group_descriptives", OUT_TAG, ".csv")), row.names = FALSE)
cat("AGE_VAR", AGE_VAR, "range", AGE_MIN, AGE_MAX, "
")
options(width = 220)
show <- T[T$random == "(1 | GUID)" & T$model %in% c("mirror_single", "ushape_single"), c("outcome_label", "cohort", "model", "estimate", "se", "p_value", "q_value", "n_nights", "n_subjects")]
show$estimate <- round(show$estimate, 3); show$se <- round(show$se, 3); show$p_value <- signif(show$p_value, 2); show$q_value <- signif(show$q_value, 2)
print(show, row.names = FALSE)
cat("\nfamily random-effect sensitivity (mirror terms):\n")
s2 <- T[T$random != "(1 | GUID)" & T$model == "mirror_single", c("outcome_label", "cohort", "estimate", "p_value")]
s2$estimate <- round(s2$estimate, 3); s2$p_value <- signif(s2$p_value, 2); print(s2, row.names = FALSE)
cat("\nfive-group vs control (primary RE):\n")
g5 <- T[T$random == "(1 | GUID)" & T$model == "group_vs_control", c("outcome_label", "term", "estimate", "p_value")]
g5$estimate <- round(g5$estimate, 2); g5$p_value <- signif(g5$p_value, 2); print(g5, row.names = FALSE)
cat("\ndescriptives:\n"); D$mean <- round(D$mean, 2); D$median <- round(D$median, 2); D$age_mean <- round(D$age_mean, 1); print(D, row.names = FALSE)
