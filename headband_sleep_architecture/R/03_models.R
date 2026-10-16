# =====================================================================
# 03_models.R -- mixed-effects models and the covariate adjustment
#
# All models:   outcome ~ <genotype term> + Age + I(Age^2) + Sex + Site
#                                         + (1 | GUID)
#
# Genotype parameterisations, per locus (copy number 1 / 2 / 3):
#
#   mirror   copy_c   Del = -1, Control =  0, Dup = +1   (linear / dose)
#   ushape   copy_q   Del = +1, Control =  0, Dup = +1   (both tails vs middle)
#   joint    copy_c + copy_q                             (orthogonal in a
#                                                          balanced design)
#   group    Genotype factor, Control as reference       (Del and Dup each
#                                                          vs Control)
#
# Wald tests use REML with Satterthwaite d.f. (lmerTest); model comparisons
# use likelihood-ratio tests on ML refits.
# =====================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

# ---------------------------------------------------------------------
# Formula construction
# ---------------------------------------------------------------------

#' Build a model formula, dropping covariates that are constant in `data`.
#' (e.g. Site has a single level once the cohort is subset)
build_formula <- function(outcome, geno_terms, data,
                          covariates = COVARIATES, random = RANDOM_TERM) {
  keep <- character(0)
  for (cv in covariates) {
    v <- all.vars(stats::as.formula(paste("~", cv)))[1]
    if (!v %in% names(data)) next
    x <- data[[v]]
    nlev <- if (is.factor(x) || is.character(x)) length(unique(x[!is.na(x)]))
            else length(unique(x[!is.na(x)]))
    if (nlev >= 2) keep <- c(keep, cv)
  }
  rhs <- c(geno_terms, keep, random)
  rhs <- rhs[nzchar(rhs)]
  stats::as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

#' Fit a linear mixed model, returning NULL (with a warning) on failure.
fit_lmm <- function(formula, data, REML = TRUE) {
  ctrl <- lme4::lmerControl(check.conv.singular = "ignore")
  fit <- tryCatch(
    suppressMessages(lmerTest::lmer(formula, data = data, REML = REML,
                                    control = ctrl)),
    error = function(e) { warning("Model failed: ", deparse(formula),
                                  " -- ", conditionMessage(e), call. = FALSE); NULL }
  )
  fit
}

# ---------------------------------------------------------------------
# Covariate adjustment  -- the fix
# ---------------------------------------------------------------------
#
# For each outcome we fit ONE reference model on the full five-group cohort:
#
#     y ~ Genotype + Age + I(Age^2) + Sex + Site + (1 | GUID)
#
# The genotype term is saturated (a free mean per group), so the covariate
# coefficients are estimated *within* genotype and are not contaminated by
# group differences.  Adjusted values remove only the covariate part, each
# covariate column centred at its own sample mean:
#
#     y_adj[i] = y[i] - sum_k ( X[i,k] - mean(X[,k]) ) * beta[k]
#                       over covariate columns k only
#
# Properties this buys us:
#   * genotype differences, subject random effects and residual scatter are
#     all PRESERVED in the dots (unlike `y - predict(re.form = NA)`, which
#     removes the genotype effect);
#   * mean(y_adj) == mean(y), so the vertical scale is unchanged;
#   * the value is a property of the NIGHT, not of the panel, so a given
#     night is plotted at the same height in the single-genotype figure,
#     the mirror figure and the joint figure.
#
# Model-implied means are computed at the SAME covariate reference (each
# covariate column at its sample mean, i.e. proportional weighting over sex
# and site), so the overlaid line is on the same scale as the dots.
# ---------------------------------------------------------------------

#' Which fixed-effect columns encode genotype (as opposed to covariates)?
geno_cols_of <- function(fit, geno_terms) {
  cn <- colnames(lme4::getME(fit, "X"))
  hit <- rep(FALSE, length(cn))
  for (g in geno_terms) {
    v <- all.vars(stats::as.formula(paste("~", g)))[1]
    hit <- hit | (cn == g) | startsWith(cn, v)
  }
  cn[hit]
}

#' Covariate-adjusted values for every row used by `fit`.
#' Returns a data.frame: .row (row name in the model frame) + adj value.
covariate_adjusted <- function(fit, geno_terms, outcome) {
  X    <- lme4::getME(fit, "X")
  beta <- lme4::fixef(fit)
  mf   <- stats::model.frame(fit)

  gcols <- geno_cols_of(fit, geno_terms)
  covcols <- setdiff(colnames(X), c("(Intercept)", gcols))

  if (length(covcols)) {
    Xc  <- X[, covcols, drop = FALSE]
    Xc  <- sweep(Xc, 2, colMeans(Xc), "-")
    shift <- as.numeric(Xc %*% beta[covcols])
  } else {
    shift <- rep(0, nrow(X))
  }

  data.frame(
    .row = rownames(mf),
    adj  = mf[[outcome]] - shift,
    stringsAsFactors = FALSE
  )
}

#' Attach `adj_<outcome>` to `d` for every primary/secondary outcome, using a
#' single reference model per outcome fitted on the full five-group cohort.
add_adjusted_outcomes <- function(d, outcomes, verbose = TRUE) {
  for (o in outcomes) {
    y <- o$name
    if (!y %in% names(d)) next
    sub <- d[!is.na(d[[y]]), , drop = FALSE]
    if (nrow(sub) < 20 || length(unique(sub$Genotype)) < 2) next

    f   <- build_formula(y, "Genotype", sub)
    fit <- fit_lmm(f, sub, REML = TRUE)
    if (is.null(fit)) next

    a <- covariate_adjusted(fit, "Genotype", y)
    col <- paste0("adj_", y)
    d[[col]] <- NA_real_
    # align by the model frame's row names, which are the row names of `sub`
    idx <- match(a$.row, rownames(sub))
    d[[col]][match(rownames(sub)[idx], rownames(d))] <- a$adj

    if (verbose) {
      msg("  adjusted %-12s from reference model on %d nights / %d subjects",
          y, nrow(a), length(unique(as.character(sub$GUID))))
    }
    attr(d, paste0("refmodel_", y)) <- f
  }
  d
}

#' Model-implied means at the covariate reference (all covariate columns at
#' their sample means => proportional weighting over sex and site).
#'
#' `targets` is a data.frame whose columns are genotype design-matrix columns.
model_means <- function(fit, geno_terms, targets, level = 0.95) {
  X    <- lme4::getME(fit, "X")
  beta <- lme4::fixef(fit)
  V    <- as.matrix(stats::vcov(fit))

  base <- colMeans(X)
  base["(Intercept)"] <- 1

  gcols <- geno_cols_of(fit, geno_terms)
  unknown <- setdiff(names(targets), gcols)
  if (length(unknown)) {
    stop("model_means(): target column(s) not in the design matrix: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }

  out <- lapply(seq_len(nrow(targets)), function(i) {
    x <- base
    x[gcols] <- 0                       # start from the reference level
    for (nm in names(targets)) x[nm] <- targets[[nm]][i]
    est <- sum(x * beta)
    se  <- sqrt(as.numeric(t(x) %*% V %*% x))
    df  <- tryCatch(stats::df.residual(fit), error = function(e) Inf)
    if (is.null(df) || is.na(df)) df <- Inf
    crit <- stats::qt(1 - (1 - level) / 2, df)
    data.frame(estimate = est, se = se,
               lower = est - crit * se, upper = est + crit * se)
  })
  cbind(targets, do.call(rbind, out))
}

# ---------------------------------------------------------------------
# Test extraction
# ---------------------------------------------------------------------

#' Convergence / degeneracy flags.
#'
#' `singular = TRUE` usually means the random-intercept variance collapsed to
#' zero, i.e. this outcome shows no between-participant variance beyond
#' residual noise.  Fixed effects remain interpretable, but the row is
#' flagged so it can be checked.
.model_flags <- function(fit) {
  oi <- fit@optinfo
  list(
    singular  = isTRUE(lme4::isSingular(fit)),
    converged = is.null(oi$conv$lme4$messages) &&
                (is.null(oi$conv$opt) || identical(as.integer(oi$conv$opt), 0L))
  )
}

#' Pull Wald rows for the named terms out of a fitted model.
wald_rows <- function(fit, terms, outcome, cohort, model, data) {
  if (is.null(fit)) return(NULL)
  co <- suppressWarnings(summary(fit)$coefficients)
  fl <- .model_flags(fit)
  rows <- lapply(terms, function(tm) {
    if (!tm %in% rownames(co)) return(NULL)
    data.frame(
      outcome = outcome, cohort = cohort, model = model,
      term = tm, test = "Wald (Satterthwaite)",
      estimate = unname(co[tm, "Estimate"]),
      se       = unname(co[tm, "Std. Error"]),
      df       = unname(co[tm, "df"]),
      statistic = unname(co[tm, "t value"]),
      p_value  = unname(co[tm, "Pr(>|t|)"]),
      n_nights = nrow(data),
      n_subjects = length(unique(as.character(data$GUID))),
      singular = fl$singular, converged = fl$converged,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Likelihood-ratio test between two nested models (refit under ML).
lrt_row <- function(full, reduced, outcome, cohort, model, label, data) {
  if (is.null(full) || is.null(reduced)) return(NULL)
  a <- tryCatch(suppressMessages(stats::anova(reduced, full)),
                error = function(e) NULL)
  if (is.null(a) || nrow(a) < 2) return(NULL)
  data.frame(
    outcome = outcome, cohort = cohort, model = model,
    term = label, test = "LRT",
    estimate = NA_real_, se = NA_real_,
    df = unname(a$Df[2]),
    statistic = unname(a$Chisq[2]),
    p_value = unname(a$`Pr(>Chisq)`[2]),
    n_nights = nrow(data),
    n_subjects = length(unique(as.character(data$GUID))),
    singular = NA, converged = NA,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------
# The full test battery
# ---------------------------------------------------------------------

#' Cohort subsets used for the per-locus models: that locus's Del + Dup,
#' plus all Controls.
locus_subset <- function(d, locus) {
  keep <- as.character(d$Genotype) %in%
    c(paste0(locus, c("Del", "Dup")), "Control")
  s <- d[keep, , drop = FALSE]
  s$Genotype <- factor(as.character(s$Genotype),
                       levels = c("Control", paste0(locus, c("Del", "Dup"))))
  s$GUID <- droplevels(s$GUID)
  s$Sex  <- droplevels(s$Sex)
  s$Site <- droplevels(s$Site)
  s
}

#' Run every requested model for one outcome; returns a data.frame of tests
#' and (invisibly) the fitted objects needed for plotting.
run_outcome <- function(d, outcome, verbose = TRUE) {
  sub_all <- d[!is.na(d[[outcome]]), , drop = FALSE]
  if (nrow(sub_all) < 20) {
    warning("Too few observations for ", outcome, call. = FALSE)
    return(list(tests = NULL, fits = list()))
  }
  if (verbose) msg("  [%s] %d nights, %d subjects", outcome, nrow(sub_all),
                   length(unique(as.character(sub_all$GUID))))

  res  <- list()
  fits <- list()

  # ---- (A) all five groups vs Control, one model -------------------
  f_all  <- build_formula(outcome, "Genotype", sub_all)
  m_all  <- fit_lmm(f_all, sub_all)
  f_null <- build_formula(outcome, character(0), sub_all)
  m_null <- fit_lmm(f_null, sub_all)

  fits[["all_group"]] <- m_all
  fits[["all_data"]]  <- sub_all

  if (!is.null(m_all)) {
    tms <- grep("^Genotype", rownames(summary(m_all)$coefficients), value = TRUE)
    res[[length(res) + 1]] <- wald_rows(m_all, tms, outcome, "All5",
                                        "group_vs_control", sub_all)
    res[[length(res) + 1]] <- lrt_row(m_all, m_null, outcome, "All5",
                                      "group_vs_control", "omnibus (4 df)", sub_all)
  }

  # ---- (B) per-locus mirror / U-shape / joint / group ---------------
  for (locus in c("16p", "22q")) {
    s <- locus_subset(d, locus)
    s <- s[!is.na(s[[outcome]]), , drop = FALSE]
    if (nrow(s) < 20 || length(unique(s$CopyNum)) < 3) {
      warning("Skipping ", locus, " for ", outcome, ": insufficient data",
              call. = FALSE)
      next
    }

    m_mirror <- fit_lmm(build_formula(outcome, "copy_c", s), s)
    m_ushape <- fit_lmm(build_formula(outcome, "copy_q", s), s)
    m_joint  <- fit_lmm(build_formula(outcome, c("copy_c", "copy_q"), s), s)
    m_group  <- fit_lmm(build_formula(outcome, "Genotype", s), s)
    m_base   <- fit_lmm(build_formula(outcome, character(0), s), s)

    fits[[paste0(locus, "_mirror")]] <- m_mirror
    fits[[paste0(locus, "_ushape")]] <- m_ushape
    fits[[paste0(locus, "_joint")]]  <- m_joint
    fits[[paste0(locus, "_group")]]  <- m_group
    fits[[paste0(locus, "_data")]]   <- s

    res[[length(res) + 1]] <- wald_rows(m_mirror, "copy_c", outcome, locus,
                                        "mirror_single", s)
    res[[length(res) + 1]] <- wald_rows(m_ushape, "copy_q", outcome, locus,
                                        "ushape_single", s)
    res[[length(res) + 1]] <- wald_rows(m_joint, c("copy_c", "copy_q"),
                                        outcome, locus, "joint", s)
    tms <- if (!is.null(m_group))
      grep("^Genotype", rownames(summary(m_group)$coefficients), value = TRUE) else character(0)
    res[[length(res) + 1]] <- wald_rows(m_group, tms, outcome, locus,
                                        "group_vs_control", s)

    res[[length(res) + 1]] <- lrt_row(m_joint, m_mirror, outcome, locus, "joint",
                                      "U-shape beyond mirror (1 df)", s)
    res[[length(res) + 1]] <- lrt_row(m_joint, m_ushape, outcome, locus, "joint",
                                      "mirror beyond U-shape (1 df)", s)
    res[[length(res) + 1]] <- lrt_row(m_joint, m_base, outcome, locus, "joint",
                                      "omnibus copy number (2 df)", s)
    res[[length(res) + 1]] <- lrt_row(m_group, m_base, outcome, locus,
                                      "group_vs_control", "omnibus (2 df)", s)
  }

  tests <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  list(tests = tests, fits = fits)
}

#' Human-readable labels + FDR within test families.
finalise_results <- function(tests, outcomes, fdr_method = "BH") {
  lab <- setNames(vapply(outcomes, function(o) o$label, character(1)),
                  vapply(outcomes, function(o) o$name, character(1)))
  prim <- setNames(vapply(outcomes, function(o) isTRUE(o$primary), logical(1)),
                   vapply(outcomes, function(o) o$name, character(1)))

  tests$outcome_label  <- unname(lab[tests$outcome])
  tests$primary_outcome <- unname(prim[tests$outcome])

  tests$term_label <- tests$term
  tests$term_label[tests$term == "copy_c"] <- "Mirror (linear copy number)"
  tests$term_label[tests$term == "copy_q"] <- "U-shape (squared copy number)"
  tests$term_label <- sub("^Genotype", "", tests$term_label)
  tests$term_label[tests$model == "group_vs_control" & tests$test != "LRT"] <-
    paste0(tests$term_label[tests$model == "group_vs_control" & tests$test != "LRT"],
           " vs Control")

  # FDR families -- computed over the PRIMARY outcomes only
  tests$fdr_family <- NA_character_
  is_p <- tests$primary_outcome & tests$test != "LRT"
  tests$fdr_family[is_p & tests$model %in% c("mirror_single", "ushape_single")] <-
    "primary: mirror & U-shape (single-term models)"
  tests$fdr_family[is_p & tests$model == "joint"] <-
    "joint model terms"
  tests$fdr_family[is_p & tests$model == "group_vs_control"] <-
    "each group vs Control"

  tests$q_value <- NA_real_
  for (fam in unique(stats::na.omit(tests$fdr_family))) {
    i <- which(tests$fdr_family == fam)
    tests$q_value[i] <- stats::p.adjust(tests$p_value[i], method = fdr_method)
  }

  ord <- order(match(tests$outcome, names(lab)),
               match(tests$cohort, c("All5", "16p", "22q")),
               match(tests$model, c("mirror_single", "ushape_single",
                                    "joint", "group_vs_control")))
  tests <- tests[ord, , drop = FALSE]

  tests[, c("outcome", "outcome_label", "primary_outcome", "cohort", "model",
            "term", "term_label", "test", "estimate", "se", "df", "statistic",
            "p_value", "fdr_family", "q_value", "n_nights", "n_subjects",
            "singular", "converged")]
}
