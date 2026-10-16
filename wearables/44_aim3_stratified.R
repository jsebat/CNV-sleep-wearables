# Aim 3 follow-up: HEADBAND sleep efficiency vs Penn CNB accuracy, stratified by CNV group.
# Controls cannot be included: no family control has CNB data (cnb_model_data.tsv holds carriers only).
# Pooled reference model (as in 41_aim3_assoc.R, base): cnb_acc ~ SE + Genotype + Age_c + Age_c^2 + Sex + Site.
# Panel A: partial residuals of CNB accuracy w.r.t. sleep efficiency from that ONE reference model (covariates and genotype
#          centred and removed), coloured by group; pooled line, plus each group's own slope (within-group model:
#          Age_c + Age_c^2 + Sex; site dropped, too few per group).
# Panel B: slope per 10 percentage points of sleep efficiency, per group and pooled, 95 % CI; heterogeneity = LRT of
#          SE x Genotype in the pooled model. 16p duplication (n = 3) is shown as points only.
# Usage: Rscript 44_aim3_stratified.R <analysis_root> <out_dir>
suppressPackageStartupMessages({ library(ggplot2) })
a <- commandArgs(TRUE); ROOT <- a[1]; OUT <- a[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M <- read.csv(file.path(ROOT, "aim4", "master", "G2MH_wearable_phenotypes_per_participant.csv"))
M <- M[M$visit %in% 1 & M$link_status == "linked" & M$genotype %in% c("16pDel", "16pDup", "22qDel", "22qDup"), ]
O <- read.csv(file.path(ROOT, "aim3", "outcomes.csv"))
D <- merge(M, O[, c("empatica_id", "cnb_acc")], by = "empatica_id"); D <- D[is.finite(D$cnb_acc) & is.finite(D$hb_SE), ]
if (max(D$hb_SE) > 1.5) D$hb_SE <- D$hb_SE / 100       # stored as % in some builds
cat("hb_SE range:", range(D$hb_SE), "\n")
D$SE10 <- D$hb_SE * 10                                  # fraction x 10: 1 unit = 10 percentage points
D$Age_c <- D$age_at_recording - mean(D$age_at_recording); D$Sex <- factor(D$sex); D$Site <- factor(D$site); D$Genotype <- factor(D$genotype)
cat("n =", nrow(D), "\n"); print(table(D$genotype))
m <- lm(cnb_acc ~ SE10 + Genotype + Age_c + I(Age_c^2) + Sex + Site, D)
mi <- lm(cnb_acc ~ SE10 * Genotype + Age_c + I(Age_c^2) + Sex + Site, D)
p_het <- anova(m, mi)$`Pr(>F)`[2]
X <- model.matrix(m); b <- coef(m); cc <- setdiff(colnames(X), c("(Intercept)", "SE10"))
D$adj <- D$cnb_acc - as.numeric(sweep(X[, cc], 2, colMeans(X[, cc])) %*% b[cc])
GN <- c("16pDel" = "16p deletion", "16pDup" = "16p duplication", "22qDel" = "22q deletion", "22qDup" = "22q duplication")
COL <- c("16p deletion" = "#d95f02", "16p duplication" = "#1b9e77", "22q deletion" = "#7570b3", "22q duplication" = "#2a78d6")
D$group <- factor(GN[D$genotype], GN)
S <- do.call(rbind, lapply(names(GN), function(g) { d <- D[D$genotype == g, ]
  if (nrow(d) < 8) return(data.frame(group = GN[g], n = nrow(d), est = NA, lo = NA, hi = NA, p = NA))
  f <- lm(cnb_acc ~ SE10 + Age_c + I(Age_c^2) + Sex, d); ci <- confint(f)["SE10", ]
  data.frame(group = GN[g], n = nrow(d), est = coef(f)["SE10"], lo = ci[1], hi = ci[2], p = summary(f)$coefficients["SE10", 4]) }))
ci <- confint(m)["SE10", ]
S <- rbind(S, data.frame(group = "All carriers (pooled)", n = nrow(D), est = b["SE10"], lo = ci[1], hi = ci[2], p = summary(m)$coefficients["SE10", 4]))
write.csv(S, file.path(OUT, "aim3_SE_cnb_by_group.csv"), row.names = FALSE); print(S); cat("heterogeneity p =", signif(p_het, 2), "\n")
# group lines in the partial-residual space: each group's own slope, through that group's mean point
L <- do.call(rbind, lapply(names(GN), function(g) { s <- S[S$group == GN[g], ]; d <- D[D$genotype == g, ]
  if (is.na(s$est)) return(NULL); x <- range(d$SE10); data.frame(group = GN[g], x = x, y = mean(d$adj) + s$est * (x - mean(d$SE10))) }))
L$group <- factor(L$group, GN)
xr <- range(D$SE10); pool <- data.frame(x = xr, y = mean(D$adj) + b["SE10"] * (xr - mean(D$SE10)))
th <- theme_minimal(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8.5, colour = "#555555"))
pa <- ggplot(D, aes(SE10 * 10, adj, colour = group)) + geom_point(size = 2.2, alpha = 0.85) +
  geom_line(data = L, aes(x * 10, y, colour = group), linewidth = 0.8) +
  geom_line(data = pool, aes(x * 10, y), colour = "black", linewidth = 1.1, linetype = "22", inherit.aes = FALSE) +
  scale_colour_manual(values = COL, name = NULL, drop = FALSE) +
  labs(x = "Headband sleep efficiency (%)", y = "Penn CNB accuracy composite (adjusted, z)", title = "Carriers: headband sleep efficiency and cognition",
       subtitle = "One pooled model (CNV group, age, age², sex, site).
Dashed = pooled slope; coloured = each group's own slope") +
  th + theme(legend.position = "bottom")
S$group <- factor(S$group, rev(c(GN, "All carriers (pooled)")))
pb <- ggplot(S[!is.na(S$est), ], aes(est, group)) + geom_vline(xintercept = 0, colour = "#999999") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) + geom_point(aes(colour = group), size = 3) +
  geom_text(aes(label = sprintf("n = %d", n)), nudge_y = 0.32, size = 3, colour = "#444444") +
  scale_colour_manual(values = c(COL, "All carriers (pooled)" = "black"), guide = "none") +
  labs(x = "Change in CNB accuracy (z) per +10 points of sleep efficiency, 95% CI", y = NULL, title = "Slope by CNV group",
       subtitle = sprintf("Heterogeneity across groups p = %.2f.
16p duplication (n = %d) not estimated; no controls have CNB data.", p_het, sum(D$genotype == "16pDup"))) + th
ggsave(file.path(OUT, "aim3_SE_cnb_scatter.png"), pa, width = 6.8, height = 5.2, dpi = 200, bg = "white")
ggsave(file.path(OUT, "aim3_SE_cnb_by_group.png"), pb, width = 6.6, height = 3.6, dpi = 200, bg = "white")
cat("done\n")
