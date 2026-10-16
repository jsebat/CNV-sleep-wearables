# Aim 3 association layer: do wearable phenotypes relate to clinical features in CNV carriers, beyond CNV group?
# Carriers only (visit 1); controls were not psychiatrically assessed.
#   base  : outcome ~ phenotype_z + Genotype (4 groups) + Age_c + Age_c^2 + Sex + Site            (logistic or linear)
#   meds  : base + any psychotropic medication (stimulant, non-stimulant ADHD, alpha-2, antipsychotic, antidepressant,
#           lithium, anxiolytic/hypnotic, antiepileptic) - carriers without a meds form dropped
#   nomeds: base, carriers on any psychotropic (or without a meds form) removed
#   interaction: LRT for phenotype x locus (16p vs 22q) added to base
# SEs: model-based. Cluster-robust (sandwich::vcovCL, env ROBUST=1) was tried and rejected: with ~50-80 carriers, ~12
#      parameters and mostly singleton families it returned implausibly small SEs (p ~1e-10 at n = 28).
#      Relatedness is handled instead by a sensitivity model keeping one carrier per family ('onefam', fixed seed).
# Phenotypes: WRIST core (primary family) and HEADBAND core (secondary family), z-scored over carriers.
# FDR (BH) within outcome x device family, over phenotypes, base model.
# Usage: Rscript 41_aim3_assoc.R <analysis_root>
suppressPackageStartupMessages({ library(stats) })
ROOT <- commandArgs(TRUE)[1]; OUT <- file.path(ROOT, "aim3")
ROB <- Sys.getenv("ROBUST") == "1" && requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)
cat("cluster-robust SEs:", ROB, "\n")
M <- read.csv(file.path(ROOT, "aim4", "master", "G2MH_wearable_phenotypes_per_participant.csv"))
M <- M[M$visit %in% 1 & M$link_status == "linked" & M$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup"), ]
O <- read.csv(file.path(OUT, "outcomes.csv"))
D <- merge(M, O[, setdiff(names(O), c("rarecnv_id", "genotype"))], by = "empatica_id")
D$Age_c <- D$age_at_recording - mean(D$age_at_recording, na.rm = TRUE); D$Sex <- factor(D$sex); D$Site <- factor(D$site)
D$Genotype <- factor(D$genotype); D$Locus <- factor(substr(D$genotype, 1, 3))
D$family <- ifelse(is.na(D$family) | D$family == "", D$study_id, D$family)
MEDS <- c("med_stimulant", "med_adhd_nonstim", "med_alpha2", "med_antipsychotic", "med_antidepressant", "med_lithium", "med_anxiolytic_hypnotic", "med_antiepileptic")
D$med_any <- ifelse(D$meds_form_present %in% 1, as.integer(rowSums(D[, MEDS], na.rm = TRUE) > 0), NA)
COUNTS <- c("wr_n_valid_days", "wr_n_valid_nights")
WR <- setdiff(grep("^wr_", names(D), value = TRUE), c(grep("_p2_", names(D), value = TRUE), COUNTS))
HB <- grep("^hb_", names(D), value = TRUE); HB <- HB[!grepl("^hb_(TRT|min_|frac_MT)", HB)]
set.seed(1); D1 <- D[sample(nrow(D)), ]; D1 <- D1[!duplicated(D1$family), ]
cat("families:", length(unique(D$family)), "| carriers in multi-carrier families:", sum(duplicated(D$family) | duplicated(D$family, fromLast = TRUE)), "
")
BIN <- c("asd", "adhd", "id", "psychosis", "anxiety", "depression", "ocd"); CONT <- c("fsiq", "cnb_acc", "cnb_spd")
PRIMARY <- c("asd", "adhd", "id", "psychosis", "fsiq")
cat("carriers with wrist data:", sum(is.finite(D$wr_tst_h)), "| with headband:", sum(is.finite(D$hb_TST)),
    "| meds form:", sum(!is.na(D$med_any)), "| on psychotropic:", sum(D$med_any %in% 1), "\n")

fit1 <- function(y, x, d, extra = character(0), inter = FALSE) {
  d <- d[is.finite(d[[y]]) & is.finite(d[[x]]), ]
  for (e in extra) d <- d[!is.na(d[[e]]), ]
  if (nrow(d) < 25) return(NULL)
  isb <- y %in% BIN; if (isb && (sum(d[[y]]) < 5 || sum(1 - d[[y]]) < 5)) return(NULL)
  d$X <- (d[[x]] - mean(d[[x]])) / sd(d[[x]])
  cv <- c("Genotype", "Age_c", "I(Age_c^2)", if (nlevels(droplevels(d$Sex)) > 1) "Sex", if (nlevels(droplevels(d$Site)) > 1) "Site", extra)
  f <- as.formula(paste(y, "~ X +", paste(cv, collapse = " + ")))
  fam <- if (isb) binomial() else gaussian()
  m <- tryCatch(glm(f, data = d, family = fam), error = function(e) NULL, warning = function(w) suppressWarnings(glm(f, data = d, family = fam)))
  if (is.null(m) || !"X" %in% names(coef(m)) || is.na(coef(m)["X"])) return(NULL)
  co <- if (ROB) lmtest::coeftest(m, vcov = sandwich::vcovCL(m, cluster = d$family))["X", ] else summary(m)$coefficients["X", ]
  pint <- NA
  if (inter && nlevels(droplevels(d$Locus)) > 1) {
    m1 <- tryCatch(suppressWarnings(glm(update(f, . ~ . + X:Locus), data = d, family = fam)), error = function(e) NULL)
    if (!is.null(m1)) pint <- anova(m, m1, test = if (isb) "LRT" else "F")[2, ncol(anova(m, m1, test = if (isb) "LRT" else "F"))]
  }
  data.frame(estimate = co[1], se = co[2], p = co[4], n = nrow(d), n_cases = if (isb) sum(d[[y]]) else NA,
             effect = if (isb) exp(co[1]) else co[1] / sd(d[[y]]), p_interaction_locus = pint,
             sep_warning = isb && any(fitted(m) < 1e-6 | fitted(m) > 1 - 1e-6))
}
res <- list()
for (y in c(BIN, CONT)) for (fam_ in c("wrist", "headband")) for (x in if (fam_ == "wrist") WR else HB) {
  for (mod in c("base", "meds", "nomeds", "onefam")) {
    d <- if (mod == "nomeds") D[D$med_any %in% 0, ] else if (mod == "onefam") D1 else D
    r <- fit1(y, x, d, extra = if (mod == "meds") "med_any" else character(0), inter = mod == "base")
    if (!is.null(r)) res[[length(res) + 1]] <- cbind(data.frame(outcome = y, primary = y %in% PRIMARY, device = fam_, phenotype = x, model = mod), r)
  }
}
T <- do.call(rbind, res)
T$q <- NA; k <- T$model == "base"
T$q[k] <- ave(T$p[k], paste(T$outcome[k], T$device[k]), FUN = function(p) p.adjust(p, "BH"))
T$q_all_primary <- NA; k2 <- k & T$primary & T$device == "wrist"; T$q_all_primary[k2] <- p.adjust(T$p[k2], "BH")
write.csv(T, file.path(OUT, "aim3_association.csv"), row.names = FALSE)
options(width = 230)
B <- T[T$model == "base", ]
cat("\nbase-model tests:", nrow(B), "| q < 0.05 within outcome x device:", sum(B$q < 0.05), "| primary wrist family q < 0.05:", sum(B$q_all_primary < 0.05, na.rm = TRUE), "\n")
top <- B[order(B$p), ][1:25, c("outcome", "device", "phenotype", "effect", "p", "q", "q_all_primary", "n", "n_cases", "p_interaction_locus")]
top[, c("effect")] <- round(top$effect, 2); top$p <- signif(top$p, 2); top$q <- signif(top$q, 2); top$q_all_primary <- signif(top$q_all_primary, 2); top$p_interaction_locus <- signif(top$p_interaction_locus, 2)
print(top, row.names = FALSE)
# how the top hits behave under the medication models
cat("\nmedication checks for the 10 strongest base associations:\n")
for (i in seq_len(10)) { r <- B[order(B$p), ][i, ]
  z <- T[T$outcome == r$outcome & T$phenotype == r$phenotype, c("model", "effect", "p", "n")]; z$effect <- round(z$effect, 2); z$p <- signif(z$p, 2)
  cat(r$outcome, r$phenotype, ":", paste(apply(z, 1, function(v) paste0(v["model"], " ", v["effect"], " (p ", v["p"], ", n ", v["n"], ")")), collapse = " | "), "\n") }
