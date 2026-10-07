# Kaplan-Meier plotting and MICE-pooled logistic / Cox regression.
# Packages survival, survminer and mice are called with `pkg::`, so notebooks
# only need library(mice) / library(survival) for their own direct calls.

# ---- Kaplan-Meier -------------------------------------------------------------

#' PFS Kaplan-Meier curves by group, in the paper's print style.
#'
#' Used for Fig1B, the Fig2 LN split and SFig6B. Time and event must already be
#' computed: pass load_sample_metadata() rows with km_include == TRUE (infusion
#' before 2024-12-01). pfs_time_km_days sets censored patients to 500 days, past
#' the 400-day window, so no censoring marks show.
#'
#' @param df         data.frame with one row per patient.
#' @param group      grouping column; its factor levels set the curve order.
#' @param palette    colours in curve order; if its names cover all levels, it is
#'                   matched by name instead.
#' @param time,event PFS time (days) and event (1 = progression) columns.
#' @param linetype   optional line types in curve order (Fig2: solid / twodash).
#' @param title      optional plot title.
#' @param major_grid passed to theme_cfrna_print_no_minor(): TRUE for Fig1,
#'                   FALSE for the Fig2 / SFig6B panels.
#' @param max_days   right end of the x axis (days).
#' @param risk_table add the survminer risk table (off in all paper panels).
#' @return a ggsurvplot object; save or compose its `$plot`.
km_by_group <- function(df, group, palette, time = "pfs_time_km_days", event = "pfs_event",
                        linetype = NULL, title = NULL, major_grid = TRUE,
                        max_days = 400, risk_table = FALSE) {
  km_df <- data.frame(time = df[[time]], event = df[[event]], group = factor(df[[group]]))
  if (all(levels(km_df$group) %in% names(palette))) palette <- palette[levels(km_df$group)]

  fit <- survminer::surv_fit(survival::Surv(time, event) ~ group, data = km_df)

  args <- list(fit, data = km_df,
               palette = unname(palette),
               legend = "none",
               xlab = "Months",
               ylab = "PFS Fraction",
               conf.int = FALSE,
               censor = TRUE,
               pval = FALSE,
               risk.table = risk_table,
               xlim = c(0, max_days),
               xscale = "d_m",
               break.x.by = 30.44 * 2,
               ggtheme = theme_cfrna_print_no_minor(major_grid))
  if (!is.null(linetype)) args$linetype <- linetype
  if (!is.null(title)) args$title <- title
  do.call(survminer::ggsurvplot, args)
}

# ---- Pooling across multiple imputations ----------------------------------------

# Column names of a tidied pooled fit: "OR"/"HR" (SFig9, Tables S3/S4) or
# "broom" (estimate/conf.low/conf.high/p.value; SFig2).
.pooled_cols <- function(col_style) {
  switch(col_style,
         OR    = c("OR", "OR_lo", "OR_hi", "p_value"),
         HR    = c("HR", "HR_lo", "HR_hi", "p_value"),
         broom = c("estimate", "conf.low", "conf.high", "p.value"))
}

# Flatten summary(<mipo>) into term + exponentiated estimate, 95% CI and p-value.
.tidy_pooled <- function(pooled, col_style) {
  sm <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)
  out <- tibble::tibble(term = as.character(sm$term), sm$estimate, sm$`2.5 %`,
                        sm$`97.5 %`, sm$p.value)
  colnames(out) <- c("term", .pooled_cols(col_style))
  out
}

# One all-NA row, returned when a model cannot be fitted.
.na_row <- function(term, col_style, extra) {
  cols <- .pooled_cols(col_style)
  row <- tibble::tibble(term = term, !!!setNames(as.list(rep(NA_real_, 4)), cols))
  dplyr::bind_cols(row, tibble::as_tibble(extra))
}

#' Tidy a pooled logistic fit (Rubin's rules); odds ratios, intercept dropped.
#' @param pooled   result of mice::pool().
#' @param col_style "OR" (OR, OR_lo, OR_hi, p_value) or "broom".
tidy_pooled_logistic <- function(pooled, col_style = "OR") {
  dplyr::filter(.tidy_pooled(pooled, col_style), term != "(Intercept)")
}

#' Tidy a pooled Cox fit (Rubin's rules); hazard ratios.
#' @param col_style "HR" (HR, HR_lo, HR_hi, p_value) or
#'   "broom" (estimate, conf.low, conf.high, p.value; SFig2).
tidy_pooled_cox <- function(pooled, col_style = "HR") {
  .tidy_pooled(pooled, col_style)
}

#' Fit `fitter(fml, data)` on each completed data set of a mids object.
#' @return a mira object, ready for mice::pool().
fit_across_imputations <- function(fml, imp, fitter) {
  fits <- lapply(seq_len(imp$m), function(i) fitter(fml, mice::complete(imp, i)))
  mice::as.mira(fits)
}

# Cox formula Surv(time, event) ~ terms [+ strata(strata)]. Surv() and strata()
# are resolved from survival even if it is not attached; "survival::strata"
# cannot be written in the formula because coxph would not treat it as strata.
.cox_formula <- function(terms, time, event, strata) {
  rhs <- paste(terms, collapse = " + ")
  if (!is.null(strata)) rhs <- paste0(rhs, " + strata(", strata, ")")
  fml <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", rhs))
  environment(fml) <- list2env(list(Surv = survival::Surv, strata = survival::strata),
                               parent = globalenv())
  fml
}

#' Pooled logistic regression of `outcome` on a score, adjusted for covariates.
#' @param score_var score column (per-SD if it was standardized).
#' @param adj_var   covariate columns, or NULL for the unadjusted model.
#' @param imp       mids object from mice::mice().
#' @return one row: term, OR, OR_lo, OR_hi, p_value, n (NA row if the fit fails).
fit_adjusted_logistic_mi <- function(score_var, adj_var, imp, outcome = "response") {
  fml <- reformulate(c(score_var, adj_var), outcome)
  tryCatch({
    fit <- fit_across_imputations(fml, imp, function(f, d) glm(f, data = d, family = binomial))
    tidy_pooled_logistic(mice::pool(fit)) %>%
      dplyr::filter(term == score_var) %>%
      dplyr::mutate(n = nrow(imp$data))
  }, error = function(e) .na_row(score_var, "OR", list(n = NA_integer_)))
}

#' Pooled Cox regression of PFS on a score, adjusted for covariates.
#' @param strata column to stratify on ("cohort" for SFig9), or NULL.
#' @return one row: term, HR, HR_lo, HR_hi, p_value, n, n_events.
fit_adjusted_cox_mi <- function(score_var, adj_var, imp, time = "pfs_time",
                                event = "pfs_event", strata = "cohort", col_style = "HR") {
  fml <- .cox_formula(c(score_var, adj_var), time, event, strata)
  tryCatch({
    fit <- fit_across_imputations(fml, imp, function(f, d) survival::coxph(f, data = d))
    tidy_pooled_cox(mice::pool(fit), col_style) %>%
      dplyr::filter(term == score_var) %>%
      dplyr::mutate(n = nrow(imp$data), n_events = sum(imp$data[[event]]))
  }, error = function(e) .na_row(score_var, col_style,
                                 list(n = NA_integer_, n_events = NA_integer_)))
}

#' Pooled univariate logistic regression of `outcome` on one variable.
#'
#' Used for SFig9B / Table S4. Same model as fit_adjusted_logistic_mi(variable, NULL, imp),
#' but all non-intercept terms are kept, so a factor gives one row per
#' non-reference level (e.g. "Bridging1"). For standardized scores the OR is per SD.
#' @return term, OR, OR_lo, OR_hi, p_value, variable, n (NA row if the fit fails).
fit_univariate_logistic_mi <- function(variable, imp, outcome = "response") {
  fml <- reformulate(variable, outcome)
  tryCatch({
    fit <- fit_across_imputations(fml, imp, function(f, d) glm(f, data = d, family = binomial))
    tidy_pooled_logistic(mice::pool(fit)) %>%
      dplyr::mutate(variable = variable, n = nrow(imp$data))
  }, error = function(e) .na_row(variable, "OR", list(variable = variable, n = NA_integer_)))
}

#' Pooled univariate Cox regression of PFS on one variable.
#'
#' Companion of fit_univariate_logistic_mi() (stratified by cohort like
#' fit_adjusted_cox_mi). With strata = NULL, time = "pd_rel_infusion",
#' event = "censor" and col_style = "broom" it gives the univariate Cox models of
#' SFig2. A failed fit returns an NA row.
#' @return term, HR, HR_lo, HR_hi, p_value, variable, n, n_events.
fit_univariate_cox_mi <- function(variable, imp, time = "pfs_time", event = "pfs_event",
                                  strata = "cohort", col_style = "HR") {
  fml <- .cox_formula(variable, time, event, strata)
  tryCatch({
    fit <- fit_across_imputations(fml, imp, function(f, d) survival::coxph(f, data = d))
    tidy_pooled_cox(mice::pool(fit), col_style) %>%
      dplyr::mutate(variable = variable, n = nrow(imp$data),
                    n_events = sum(imp$data[[event]]))
  }, error = function(e) .na_row(variable, col_style,
                                 list(variable = variable, n = NA_integer_,
                                      n_events = NA_integer_)))
}

# ---- Labels -------------------------------------------------------------------

#' Display labels for model terms in the SFig9 forest plots.
BASELINE_LABEL_MAP <- c(
  LN_set1 = "LN", FMAC_set1 = "FMAC", TEX_set1 = "TEX",
  IPI = "IPI", LDH = "LDH (log2)", LDH_D30 = "LDH Day-30 (log2)",
  Bridging1 = "Bridging Therapy", Stage1 = "Advanced Stage (III-IV)",
  Lymphoma_malignancy = "Lymphoma Malignancy"
)

#' Map model terms to display labels; unmapped terms are returned unchanged.
clean_baseline_term <- function(x, label_map = BASELINE_LABEL_MAP) {
  dplyr::coalesce(unname(label_map[as.character(x)]), as.character(x))
}
