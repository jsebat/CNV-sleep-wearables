# =====================================================================
# 06_tradeoff.R -- is there a depth / duration trade-off, and is it
#                  specific to CNV carriers?
#
# QUESTION
#   Do longer nights have proportionally less deep sleep, and does that
#   relationship differ between controls and 16p11.2 / 22q11.2 carriers?
#
# IS IT DEFINITIONAL?
#   DI = N3 / (TST - N3), so TST appears in the denominator.  But DI is
#   scale-invariant in TST: if N3 = c*TST then DI = c/(1-c), independent of
#   TST.  So a DI-TST association is NOT forced by the algebra -- it means
#   sleep *composition* changes with duration.  As an additional check that
#   is immune to any ratio artefact, `tradeoff_loglog()` regresses
#   log(N3 minutes) on log(TST minutes): a slope of 1 is exact
#   proportionality (no trade-off), < 1 is sub-linear (a trade-off).
#
# WITHIN vs BETWEEN
#   Nights are nested in participants, so the raw association mixes two
#   different questions.  Person-mean centring separates them:
#     TST_w = TST(night) - mean TST(participant)   -> within-participant
#     TST_b = mean TST(participant) - grand mean   -> between-participant
#   Both go in the same model (Hoffman & Stawski).  They can differ in sign,
#   so pooling them would be misleading.
#
#   Slopes are expressed per HOUR of total sleep time.
# =====================================================================

TRADEOFF_GROUPS <- c("Control", "16pDel", "16pDup", "22qDel", "22qDup")

#' Add person-mean-centred TST (and a generic within/between pair).
add_person_centered <- function(d, var = "TST", per = 60) {
  gm <- tapply(d[[var]], d$GUID, mean, na.rm = TRUE)
  pm <- as.numeric(gm[as.character(d$GUID)])
  d[[paste0(var, "_w")]] <- (d[[var]] - pm) / per
  d[[paste0(var, "_b")]] <- (pm - mean(d[[var]], na.rm = TRUE)) / per
  d[[paste0(var, "_pm")]] <- pm
  d
}

#' Fit the within/between model, preferring a random slope for the
#' within-participant term and falling back to a random intercept if that is
#' singular or fails to converge.  Returns the fit plus which was used.
fit_tradeoff <- function(formula_fixed, data, slope_term = "TST_w") {
  f_rs <- stats::as.formula(paste(deparse(formula_fixed, width.cutoff = 500),
                                  collapse = ""))
  rs <- fit_lmm(stats::update(f_rs, paste("~ . - (1 | GUID) + (1 +", slope_term, "| GUID)")),
                data)
  if (!is.null(rs) && !isTRUE(lme4::isSingular(rs)) &&
      is.null(rs@optinfo$conv$lme4$messages)) {
    return(list(fit = rs, re = "random slope"))
  }
  list(fit = fit_lmm(f_rs, data), re = "random intercept")
}

.co <- function(fit, term) {
  if (is.null(fit)) return(NULL)
  co <- suppressWarnings(summary(fit)$coefficients)
  if (!term %in% rownames(co)) return(NULL)
  co[term, ]
}

#' Per-group within- and between-participant slopes of `outcome` on TST.
tradeoff_by_group <- function(d, outcome = "DI", groups = TRADEOFF_GROUPS,
                              verbose = TRUE) {
  rows <- list()
  for (g in groups) {
    s <- d[as.character(d$Genotype) == g & !is.na(d[[outcome]]), , drop = FALSE]
    if (!nrow(s)) next
    s <- add_person_centered(s, "TST")
    s$GUID <- droplevels(factor(s$GUID))
    s$Sex  <- droplevels(factor(s$Sex))
    s$Site <- droplevels(factor(s$Site))

    f <- build_formula(outcome, c("TST_w", "TST_b"), s)
    r <- fit_tradeoff(f, s)
    if (is.null(r$fit)) next

    for (tm in c("TST_w", "TST_b")) {
      cc <- .co(r$fit, tm)
      if (is.null(cc)) next
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, group = g,
        effect = if (tm == "TST_w") "within-participant" else "between-participant",
        term = tm,
        estimate = unname(cc[["Estimate"]]), se = unname(cc[["Std. Error"]]),
        df = unname(cc[["df"]]), statistic = unname(cc[["t value"]]),
        p_value = unname(cc[["Pr(>|t|)"]]),
        ci_lo = unname(cc[["Estimate"]] - 1.96 * cc[["Std. Error"]]),
        ci_hi = unname(cc[["Estimate"]] + 1.96 * cc[["Std. Error"]]),
        n_nights = nrow(s), n_subjects = length(unique(as.character(s$GUID))),
        random_effects = r$re,
        stringsAsFactors = FALSE
      )
    }
    if (verbose) msg("  %-8s n=%3d nights / %2d subjects (%s)", g, nrow(s),
                     length(unique(as.character(s$GUID))), r$re)
  }
  do.call(rbind, rows)
}

#' Does the trade-off differ by genotype?  Fits the pooled model plus three
#' interaction contrasts, each tested by LRT against its own no-interaction
#' reduced model.
tradeoff_interactions <- function(d, outcome = "DI", verbose = TRUE) {
  s <- d[!is.na(d[[outcome]]), , drop = FALSE]
  s <- add_person_centered(s, "TST")
  s$Carrier <- factor(ifelse(as.character(s$Genotype) == "Control",
                             "Control", "Carrier"),
                      levels = c("Control", "Carrier"))
  rows <- list()

  # ---- (a) pooled trade-off, all participants ------------------------
  f0 <- build_formula(outcome, c("TST_w", "TST_b"), s)
  r0 <- fit_tradeoff(f0, s)
  for (tm in c("TST_w", "TST_b")) {
    cc <- .co(r0$fit, tm)
    if (is.null(cc)) next
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, comparison = "pooled (all participants)",
      term = tm,
      effect = if (tm == "TST_w") "within-participant" else "between-participant",
      estimate = unname(cc[["Estimate"]]), se = unname(cc[["Std. Error"]]),
      statistic = unname(cc[["t value"]]), df = unname(cc[["df"]]),
      p_value = unname(cc[["Pr(>|t|)"]]), test = "Wald",
      n_nights = nrow(s), n_subjects = length(unique(as.character(s$GUID))),
      stringsAsFactors = FALSE)
  }

  # ---- (b) does the within slope differ by group? --------------------
  spec <- list(
    list(lab = "carrier vs control", inter = "TST_w:Carrier", grp = "Carrier"),
    list(lab = "5 genotype groups",  inter = "TST_w:Genotype", grp = "Genotype")
  )
  for (sp in spec) {
    full <- fit_lmm(build_formula(outcome,
              c("TST_w", "TST_b", sp$grp, sp$inter), s), s)
    red  <- fit_lmm(build_formula(outcome,
              c("TST_w", "TST_b", sp$grp), s), s)
    lr <- lrt_row(full, red, outcome, "tradeoff", sp$lab,
                  paste("TST_w x", sp$grp, "interaction"), s)
    if (!is.null(lr)) {
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, comparison = sp$lab,
        term = paste("TST_w x", sp$grp), effect = "interaction",
        estimate = NA_real_, se = NA_real_,
        statistic = lr$statistic, df = lr$df, p_value = lr$p_value,
        test = "LRT", n_nights = nrow(s),
        n_subjects = length(unique(as.character(s$GUID))),
        stringsAsFactors = FALSE)
    }
  }

  # ---- (c) does the within slope vary with copy number, per locus? ---
  for (locus in c("16p", "22q")) {
    sl <- locus_subset(s, locus)
    sl <- add_person_centered(sl, "TST")
    full <- fit_lmm(build_formula(outcome,
              c("TST_w", "TST_b", "copy_c", "TST_w:copy_c"), sl), sl)
    red  <- fit_lmm(build_formula(outcome,
              c("TST_w", "TST_b", "copy_c"), sl), sl)
    cc <- .co(full, "TST_w:copy_c")
    if (!is.null(cc)) {
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, comparison = paste(locus, "copy number"),
        term = "TST_w x copy_c", effect = "interaction",
        estimate = unname(cc[["Estimate"]]), se = unname(cc[["Std. Error"]]),
        statistic = unname(cc[["t value"]]), df = unname(cc[["df"]]),
        p_value = unname(cc[["Pr(>|t|)"]]), test = "Wald",
        n_nights = nrow(sl), n_subjects = length(unique(as.character(sl$GUID))),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Ratio-free check: log(N3 minutes) ~ log(TST minutes), per group.
#' Slope 1 = exact proportionality (no trade-off); < 1 = sub-linear.
#' Reported with a Wald test of H0: slope = 1.
tradeoff_loglog <- function(d, groups = TRADEOFF_GROUPS, verbose = TRUE) {
  rows <- list()
  for (g in groups) {
    s <- d[as.character(d$Genotype) == g & d$min_N3 > 0 & d$TST > 0, , drop = FALSE]
    if (nrow(s) < 20) next
    s$log_N3  <- log(s$min_N3)
    s$log_TST <- log(s$TST)
    s <- add_person_centered(s, "log_TST", per = 1)
    s$GUID <- droplevels(factor(s$GUID))
    s$Sex  <- droplevels(factor(s$Sex)); s$Site <- droplevels(factor(s$Site))
    fit <- fit_lmm(build_formula("log_N3", c("log_TST_w", "log_TST_b"), s), s)
    cc <- .co(fit, "log_TST_w")
    if (is.null(cc)) next
    est <- unname(cc[["Estimate"]]); se <- unname(cc[["Std. Error"]])
    tval <- (est - 1) / se
    rows[[length(rows) + 1]] <- data.frame(
      group = g, n_nights = nrow(s),
      n_subjects = length(unique(as.character(s$GUID))),
      slope = est, se = se,
      ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
      t_vs_1 = tval,
      p_vs_1 = 2 * stats::pt(abs(tval), df = unname(cc[["df"]]), lower.tail = FALSE),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Plain correlations between DI and TST, for the combined sample and each
#' group, at three levels:
#'   night level      -- all nights pooled (inflated N: nights are not
#'                       independent, so the p-value is anti-conservative;
#'                       reported for completeness)
#'   participant level-- one point per participant (their mean DI, mean TST)
#'   within-participant-- person-mean-centred values, i.e. the night-to-night
#'                       trade-off with stable between-person differences
#'                       removed
tradeoff_correlations <- function(d, outcome = "DI",
                                  groups = TRADEOFF_GROUPS, verbose = TRUE) {
  one <- function(x, y, level, g, n_subj) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 4) return(NULL)
    pe <- suppressWarnings(stats::cor.test(x, y, method = "pearson"))
    sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman",
                                           exact = FALSE))
    data.frame(
      outcome = outcome, group = g, level = level, n = length(x),
      n_subjects = n_subj,
      pearson_r = unname(pe$estimate),
      pearson_ci_lo = if (!is.null(pe$conf.int)) pe$conf.int[1] else NA_real_,
      pearson_ci_hi = if (!is.null(pe$conf.int)) pe$conf.int[2] else NA_real_,
      pearson_p = pe$p.value,
      spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
      stringsAsFactors = FALSE)
  }

  rows <- list()
  for (g in c("ALL (combined sample)", groups)) {
    s <- if (startsWith(g, "ALL")) d else d[as.character(d$Genotype) == g, ]
    s <- s[!is.na(s[[outcome]]) & !is.na(s$TST), , drop = FALSE]
    if (nrow(s) < 5) next
    ns <- length(unique(as.character(s$GUID)))

    rows[[length(rows) + 1]] <- one(s[[outcome]], s$TST, "night", g, ns)

    sm <- stats::aggregate(cbind(y = s[[outcome]], TST = s$TST),
                           by = list(GUID = as.character(s$GUID)), FUN = mean)
    rows[[length(rows) + 1]] <- one(sm$y, sm$TST, "participant", g, nrow(sm))

    sc <- add_person_centered(s, "TST", per = 1)
    ymean <- tapply(s[[outcome]], as.character(s$GUID), mean)
    y_w <- s[[outcome]] - as.numeric(ymean[as.character(s$GUID)])
    multi <- as.character(s$GUID) %in% names(which(table(as.character(s$GUID)) > 1))
    rows[[length(rows) + 1]] <- one(y_w[multi], sc$TST_w[multi],
                                    "within-participant", g, sum(multi))
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (verbose && !is.null(out)) {
    a <- out[out$group == "ALL (combined sample)", ]
    c_ <- out[out$group == "Control", ]
    for (lv in c("night", "participant", "within-participant")) {
      ra <- a$pearson_r[a$level == lv]; rc <- c_$pearson_r[c_$level == lv]
      if (length(ra) && length(rc))
        msg("  r(%s, TST) %-19s combined %+.3f | controls %+.3f", outcome, lv, ra, rc)
    }
  }
  out
}

# ---------------------------------------------------------------------
# Controls for the fact that TST enters the definition of DI
# ---------------------------------------------------------------------

#' CONTROL 1 -- fixed-window depth index.
#'
#' Compute DI using only the FIRST `window_min` minutes of *sleep* in each
#' recording, restricted to nights that reach that much sleep.  Every night's
#' DI is then computed from exactly the same amount of sleep, so duration is
#' held constant by construction and cannot enter the ratio at all.
#'
#' If longer sleepers still show a lower fixed-window DI, the relationship is
#' about sleep composition, not about the algebra of the index.
#'
#' Note this asks a slightly different question -- it compares the *early*
#' night across people -- so it is a control, not a replacement.
fixed_window_di <- function(epochs, window_min = 240) {
  if (!"dur_sec" %in% names(epochs)) epochs$dur_sec <- 30
  epochs <- epochs[order(epochs$recording_id, epochs$epoch_idx), , drop = FALSE]
  sp_stage <- split(as.character(epochs$stage), epochs$recording_id)
  sp_dur   <- split(epochs$dur_sec / 60,        epochs$recording_id)

  res <- lapply(names(sp_stage), function(rid) {
    st <- sp_stage[[rid]]; dm <- sp_dur[[rid]]
    is_sleep <- st %in% SLEEP_STAGES
    cum <- cumsum(ifelse(is_sleep, dm, 0))
    if (max(cum) < window_min)
      return(data.frame(recording_id = rid, DI_fix = NA_real_,
                        sleep_in_window = max(cum), stringsAsFactors = FALSE))
    last <- which(cum >= window_min)[1]
    keep <- seq_len(last)
    st <- st[keep]; dm <- dm[keep]
    n3 <- sum(dm[st == "N3"])
    other <- sum(dm[st %in% setdiff(SLEEP_STAGES, "N3")])
    data.frame(recording_id = rid,
               DI_fix = if (other > 0) n3 / other else NA_real_,
               sleep_in_window = window_min, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

#' CONTROL 2 -- empirical null for the mechanical artefact.
#'
#' Under the null "sleep composition does not vary with duration", each night
#' is a draw of its own number of sleep epochs from ONE common composition.
#' Simulating that directly gives the distribution of DI-TST correlation
#' produced purely by length-dependent sampling noise in the ratio.
#'
#' Comparing the observed correlation with this null quantifies exactly how
#' much of it could be an artefact of TST appearing in the denominator.
tradeoff_null_calibration <- function(d, B = 500, seed = 1, verbose = TRUE) {
  st <- c("N1", "N2", "N3", "R")
  mins <- as.matrix(d[, paste0("min_", st)])
  p <- colSums(mins, na.rm = TRUE); p <- p / sum(p)          # pooled composition
  n_ep <- round(d$TST * 2)                                   # 30 s epochs
  ok <- is.finite(n_ep) & n_ep > 0 & is.finite(d$TST)
  n_ep <- n_ep[ok]; tst <- d$TST[ok]
  obs <- stats::cor(d$DI[ok], tst, use = "complete.obs")

  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  r_null <- vapply(seq_len(B), function(b) {
    cnt <- t(vapply(n_ep, function(n) stats::rmultinom(1, n, p)[, 1], numeric(4)))
    di <- cnt[, 3] / pmax(rowSums(cnt[, c(1, 2, 4), drop = FALSE]), 1e-9)
    stats::cor(di, tst)
  }, numeric(1))
  if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)

  out <- data.frame(
    observed_r = obs,
    null_mean_r = mean(r_null), null_sd_r = stats::sd(r_null),
    null_lo = unname(stats::quantile(r_null, 0.025)),
    null_hi = unname(stats::quantile(r_null, 0.975)),
    p_vs_null = (1 + sum(r_null <= obs)) / (B + 1),
    B = B, n_nights = length(tst), stringsAsFactors = FALSE)
  if (verbose) {
    msg("  mechanical null: observed r = %+.3f, null r = %+.3f [%+.3f, %+.3f]",
        obs, out$null_mean_r, out$null_lo, out$null_hi)
  }
  out
}

#' Does the copy-number effect on DI survive adjustment for sleep duration?
#'
#' This matters: 22q11.2 duplication carriers sleep LESS and have HIGHER depth
#' index, so the headline mirror effect could in principle be nothing more than
#' the duration trade-off re-expressed.  Adding TST (split into its within- and
#' between-participant parts) to the mirror model tests exactly that.  If the
#' copy-number term survives, the depth effect is not merely a duration effect.
#' Runs in either direction: `adjust_for` is split into within- and
#' between-participant parts and added to the mirror model for `outcome`.
#'
#' NOTE on the reciprocal direction (TST adjusted for DI): because the
#' trade-off means DI is partly a *function* of TST, DI is partly downstream of
#' the exposure-outcome pair.  Adjusting TST for DI is therefore not a clean
#' causal adjustment -- it answers "is the duration signal statistically
#' separable from the depth signal", not "what is the direct effect on TST".
tradeoff_adjust_main <- function(d, outcome = "DI", adjust_for = "TST",
                                 verbose = TRUE) {
  rows <- list()
  aw <- paste0(adjust_for, "_w"); ab <- paste0(adjust_for, "_b")
  for (locus in c("16p", "22q")) {
    s <- locus_subset(d[!is.na(d[[outcome]]) & !is.na(d[[adjust_for]]), ,
                        drop = FALSE], locus)
    s <- add_person_centered(s, adjust_for, per = if (adjust_for == "TST") 60 else 1)
    for (adj in c(FALSE, TRUE)) {
      terms <- if (adj) c("copy_c", aw, ab) else "copy_c"
      fit <- fit_lmm(build_formula(outcome, terms, s), s)
      cc <- .co(fit, "copy_c")
      if (is.null(cc)) next
      rows[[length(rows) + 1]] <- data.frame(
        outcome = outcome, adjusted_for = if (adj) adjust_for else NA_character_,
        locus = locus,
        model = if (adj) sprintf("mirror + %s (within & between)", adjust_for)
                else "mirror (primary)",
        term = "copy_c",
        estimate = unname(cc[["Estimate"]]), se = unname(cc[["Std. Error"]]),
        statistic = unname(cc[["t value"]]), df = unname(cc[["df"]]),
        p_value = unname(cc[["Pr(>|t|)"]]),
        n_nights = nrow(s), n_subjects = length(unique(as.character(s$GUID))),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (verbose && !is.null(out)) {
    for (locus in c("16p", "22q")) {
      r <- out[out$locus == locus, ]
      if (nrow(r) == 2)
        msg("  %s %s mirror: %+.4g (p=%.4g) -> adjusted for %s: %+.4g (p=%.4g)",
            locus, outcome, r$estimate[1], r$p_value[1], adjust_for,
            r$estimate[2], r$p_value[2])
    }
  }
  out
}

#' Add a TST-adjusted version of `outcome` for plotting.
#'
#' Same partial-residual logic as add_adjusted_outcomes(), but the reference
#' model additionally contains the within- and between-participant TST terms,
#' so sleep duration is removed alongside age, age^2, sex and site. Genotype is
#' in the model but excluded from the subtraction, so group differences remain.
#'
#' Result: `adjTST_<outcome>`, directly comparable to `adj_<outcome>`.
add_tst_adjusted_outcome <- function(d, outcome = "DI", verbose = TRUE) {
  d <- add_person_centered(d, "TST")
  sub <- d[!is.na(d[[outcome]]) & !is.na(d$TST_w), , drop = FALSE]
  f <- build_formula(outcome, c("Genotype", "TST_w", "TST_b"), sub)
  fit <- fit_lmm(f, sub)
  col <- paste0("adjTST_", outcome)
  d[[col]] <- NA_real_
  if (is.null(fit)) return(d)
  a <- covariate_adjusted(fit, "Genotype", outcome)
  idx <- match(a$.row, rownames(sub))
  d[[col]][match(rownames(sub)[idx], rownames(d))] <- a$adj
  if (verbose)
    msg("  TST-adjusted %s built from %d nights (removes age, age^2, sex, site, TST)",
        outcome, nrow(a))
  d
}

#' Run the whole trade-off battery.
run_tradeoff <- function(d, outcome = "DI", epochs = NULL,
                         window_min = 240, verbose = TRUE) {
  if (verbose) msg("Depth / duration trade-off:")
  corr     <- tradeoff_correlations(d, outcome, verbose = verbose)
  by_group <- tradeoff_by_group(d, outcome, verbose = verbose)
  inter    <- tradeoff_interactions(d, outcome, verbose = verbose)
  loglog   <- tradeoff_loglog(d, verbose = verbose)
  nullcal  <- tradeoff_null_calibration(d, verbose = verbose)
  # both directions: DI adjusted for TST, and TST adjusted for DI
  adj_main <- rbind(
    tradeoff_adjust_main(d, outcome, "TST", verbose = verbose),
    tradeoff_adjust_main(d, "TST", outcome, verbose = verbose))

  # ---- fixed-window control -----------------------------------------
  fixed_corr <- NULL; fixed_group <- NULL; d_fix <- NULL
  if (!is.null(epochs)) {
    fw <- fixed_window_di(epochs, window_min)
    d_fix <- merge(d, fw, by = "recording_id", all.x = TRUE)
    d_fix <- d_fix[!is.na(d_fix$DI_fix), , drop = FALSE]
    if (nrow(d_fix) > 30) {
      if (verbose) msg("  fixed-window (first %g min of sleep): %d of %d nights qualify",
                       window_min, nrow(d_fix), nrow(d))
      fixed_corr  <- tradeoff_correlations(d_fix, "DI_fix", verbose = FALSE)
      fixed_group <- tradeoff_by_group(d_fix, "DI_fix", verbose = FALSE)
      if (verbose) {
        a <- fixed_corr[fixed_corr$group == "ALL (combined sample)" &
                        fixed_corr$level == "night", ]
        cn <- fixed_corr[fixed_corr$group == "Control" &
                         fixed_corr$level == "night", ]
        if (nrow(a) && nrow(cn))
          msg("  r(DI_fix, TST) night-level       combined %+.3f | controls %+.3f",
              a$pearson_r, cn$pearson_r)
      }
    }
  }

  if (verbose) {
    msg("  per-group slopes: %d rows; interaction tests: %d rows",
        if (is.null(by_group)) 0 else nrow(by_group),
        if (is.null(inter)) 0 else nrow(inter))
  }
  list(correlations = corr, by_group = by_group, interactions = inter,
       loglog = loglog, null_calibration = nullcal, adjusted_main = adj_main,
       fixed_correlations = fixed_corr, fixed_by_group = fixed_group,
       data_fixed = d_fix, window_min = window_min)
}
