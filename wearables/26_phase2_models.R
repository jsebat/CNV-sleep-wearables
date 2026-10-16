# Aim 4 SECOND PANEL models - kept separate from the core phenotypes (own outputs, own FDR families).
# Person-level outcomes: y ~ term + Age_c + I(Age_c^2) + Sex + Site + (1 | family); mirror / U / joint per locus;
# five-group vs control. FDR (BH) within each second-panel domain.
# Env: AGE_MIN / AGE_MAX / OUT_TAG.   Usage: Rscript 26_phase2_models.R <analysis_root> <phenotype_csv> <dictionary_csv> <out_prefix>
#   dictionary columns: var, domain, label  (label must name the device)
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
a <- commandArgs(TRUE); ROOT <- a[1]; PHF <- a[2]; DICF <- a[3]; PRE <- a[4]
AGE_MIN <- as.numeric(Sys.getenv("AGE_MIN", "-Inf")); AGE_MAX <- as.numeric(Sys.getenv("AGE_MAX", "Inf")); OUT_TAG <- Sys.getenv("OUT_TAG", "")
L <- read.csv(file.path(ROOT, "phase0", "linkage.csv")); COV <- read.csv(file.path(ROOT, "phase0", "covariates.csv"))
base <- merge(L[L$link_status == "linked" & L$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup", "Control"),
                c("empatica_id", "guid", "genotype", "sex_final", "redcap_data_access_group", "rarecnv_family_id")], COV[, c("empatica_id", "age_at_recording")])
base <- base[!is.na(base$age_at_recording) & base$age_at_recording >= AGE_MIN & base$age_at_recording <= AGE_MAX, ]
base$Age_c <- base$age_at_recording - mean(base$age_at_recording); base$Sex <- factor(base$sex_final); base$Site <- factor(base$redcap_data_access_group)
base$GUID <- base$guid; base$family <- ifelse(is.na(base$rarecnv_family_id) | base$rarecnv_family_id == "", base$guid, base$rarecnv_family_id)
base$Genotype <- relevel(factor(base$genotype), "Control")
PH <- read.csv(PHF); DICT <- read.csv(DICF)
fit <- function(f, d, REML = TRUE) tryCatch(suppressMessages(lmer(f, data = d, REML = REML, control = lmerControl(check.conv.singular = "ignore"))), error = function(e) NULL)
cvt <- function(d) { k <- c("Age_c", "I(Age_c^2)"); if (nlevels(droplevels(d$Sex)) > 1) k <- c(k, "Sex"); if (nlevels(droplevels(d$Site)) > 1) k <- c(k, "Site"); k }
res <- list()
for (i in seq_len(nrow(DICT))) {
  v <- DICT$var[i]; d0 <- merge(PH[, c("empatica_id", v)], base, by = "empatica_id"); d0 <- d0[is.finite(d0[[v]]), ]
  if (nrow(d0) < 30) next
  f <- function(g, d) as.formula(paste(v, "~", paste(c(g, cvt(d), "(1 | family)"), collapse = " + ")))
  out <- function(m, term, cohort, model, d) { if (is.null(m)) return(NULL); co <- summary(m)$coefficients; r <- grep(term, rownames(co), fixed = TRUE)
    data.frame(var = v, label = DICT$label[i], domain = DICT$domain[i], cohort = cohort, model = model, term = rownames(co)[r], estimate = co[r, 1],
               se = co[r, 2], p_value = co[r, 5], sd_outcome = sd(d[[v]]), n_people = length(unique(d$GUID))) }
  res[[length(res) + 1]] <- out(fit(f("Genotype", d0), d0), "Genotype", "All5", "group_vs_control", d0)
  for (l in c("16p", "22q")) { d <- d0[d0$Genotype %in% c(paste0(l, "Del"), paste0(l, "Dup"), "Control"), ]
    d$copy_c <- ifelse(grepl("Del", d$Genotype), -1, ifelse(grepl("Dup", d$Genotype), 1, 0)); d$copy_q <- abs(d$copy_c)
    res[[length(res) + 1]] <- out(fit(f("copy_c", d), d), "copy_c", l, "mirror_single", d)
    res[[length(res) + 1]] <- out(fit(f("copy_q", d), d), "copy_q", l, "ushape_single", d)
    fj <- fit(f(c("copy_c", "copy_q"), d), d, FALSE); fm <- fit(f("copy_c", d), d, FALSE); fu <- fit(f("copy_q", d), d, FALSE)
    if (!is.null(fj) && !is.null(fm) && !is.null(fu)) res[[length(res) + 1]] <- data.frame(var = v, label = DICT$label[i], domain = DICT$domain[i], cohort = l,
      model = "joint", term = c("U-shape beyond mirror", "mirror beyond U-shape"), estimate = NA, se = NA,
      p_value = c(suppressMessages(anova(fm, fj))$`Pr(>Chisq)`[2], suppressMessages(anova(fu, fj))$`Pr(>Chisq)`[2]), sd_outcome = sd(d[[v]]), n_people = length(unique(d$GUID))) }
}
T <- do.call(rbind, res); T$std_estimate <- T$estimate / T$sd_outcome; T$q_domain <- NA
for (dm in unique(T$domain)) { k <- T$domain == dm & T$model %in% c("mirror_single", "ushape_single"); T$q_domain[k] <- p.adjust(T$p_value[k], "BH") }
write.csv(T, paste0(PRE, OUT_TAG, ".csv"), row.names = FALSE)
options(width = 220)
S <- T[T$model %in% c("mirror_single", "ushape_single", "group_vs_control"), c("label", "cohort", "model", "term", "estimate", "p_value", "q_domain", "n_people")]
S$estimate <- signif(S$estimate, 3); S$p_value <- signif(S$p_value, 2); S$q_domain <- signif(S$q_domain, 2)
cat("AGE range", AGE_MIN, AGE_MAX, "\n"); print(S, row.names = FALSE)
