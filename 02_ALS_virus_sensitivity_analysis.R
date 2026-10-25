# Copyright 2026 Dr Eleni Christoforidou
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# https://www.apache.org/licenses/LICENSE-2.0
# SPDX-License-Identifier: Apache-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ALS/MND respiratory-virus mortality: prespecified follow-up sensitivity analysis
# -----------------------------------------------------------------------------
# Run AFTER 01_ALS_virus_mortality_pilot.R. Uses only existing CDC exports and
# pilot_results files; does not download anything and does not overwrite pilot.
# From the project directory:
#   source("02_ALS_virus_sensitivity_analysis.R")
#
# Outcome: monthly US deaths with motor neuron disease (G12.2) as underlying
# cause, January 2011-December 2019. Ecological analysis, not individual risk.
# Exploratory model comparisons are NOT independent confirmatory tests.
# -----------------------------------------------------------------------------

START <- as.Date("2011-01-01")
END   <- as.Date("2019-12-01")
FIRST_REQUIRED <- as.Date("2010-12-01")
OUT <- "sensitivity_results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

packages <- c("readr", "dplyr", "tidyr", "ggplot2", "lubridate", "MASS", "splines", "sandwich")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Please install packages: install.packages(c(",
  paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")
need_file <- function(path) {
  if (!file.exists(path)) stop("Required file missing: ", path,
    ". Run 01_ALS_virus_mortality_pilot.R first, with the CDC exports in raw_data/.")
  path
}

# A. READ VERIFIED PILOT OUTPUTS AND ORIGINAL CDC RSV FILE
pilot <- readr::read_csv(need_file("pilot_results/pilot_analysis_dataset.csv"),
                        show_col_types = FALSE)
flu_month <- readr::read_csv(need_file("pilot_results/flu_monthly_coverage.csv"),
                            show_col_types = FALSE)
rsv_raw <- readr::read_csv(need_file("raw_data/cdc_nrevss_rsv_52kb-ccu2.csv"),
                          col_types = readr::cols(.default = "c"),
                          show_col_types = FALSE)

pilot$month <- as.Date(pilot$month)
flu_month$month <- as.Date(flu_month$month)
need_pilot <- c("month", "deaths", "rsv_pct", "flu_pct")
need_raw <- c("repweekdate", "hhs_region", "testtype", "rsvpos", "rsvtest", "outlier")
if (!all(need_pilot %in% names(pilot))) stop("Pilot dataset lacks: ",
  paste(setdiff(need_pilot, names(pilot)), collapse = ", "))
if (!all(need_raw %in% names(rsv_raw))) stop("RSV CDC file lacks: ",
  paste(setdiff(need_raw, names(rsv_raw)), collapse = ", "))
if (!all(c("month", "positivity", "coverage") %in% names(flu_month)))
  stop("Influenza monthly coverage file has unexpected columns")
expected <- seq(START, END, by = "month")
if (nrow(pilot) != length(expected) ||
    !identical(sort(pilot$month), expected) || anyDuplicated(pilot$month))
  stop("Pilot mortality data must contain precisely 108 months, 2011-2019")
if (anyNA(pilot$deaths) || any(pilot$deaths < 0)) stop("Invalid pilot death counts")

# B. PARSE RSV CDC WEEK-ENDING DATES (e.g. 10JUL2010) WITHOUT LOCALE DEPENDENCE
s <- toupper(trimws(rsv_raw$repweekdate))
valid_dates <- grepl("^[0-3][0-9][A-Z]{3}(19|20)[0-9]{2}$", s)
if (!all(valid_dates)) stop("Unexpected RSV report date format; first bad: ",
                          s[which(!valid_dates)[1]])
mon <- match(substr(s, 3, 5), toupper(month.abb))
if (anyNA(mon)) stop("Invalid month abbreviations in RSV date column")
week_ending <- as.Date(sprintf("%s-%02d-%s", substr(s, 6, 9), mon, substr(s, 1, 2)))
if (anyNA(week_ending)) stop("Unparseable RSV report dates")
rsv_raw$week_start <- week_ending - 6L  # NREVSS weeks are Sunday-Saturday
rsv_raw$method <- tolower(trimws(rsv_raw$testtype))
rsv_raw$region <- suppressWarnings(as.integer(rsv_raw$hhs_region))
rsv_raw$tests <- suppressWarnings(as.numeric(rsv_raw$rsvtest))
rsv_raw$positive <- suppressWarnings(as.numeric(rsv_raw$rsvpos))
rsv_raw$flag <- suppressWarnings(as.integer(rsv_raw$outlier))
if (anyNA(rsv_raw[, c("region", "tests", "positive", "flag")]) ||
    any(!rsv_raw$method %in% c("pcr", "antigen")) ||
    any(!rsv_raw$region %in% 1:10) ||
    any(!rsv_raw$flag %in% c(0L, 1L)) ||
    any(rsv_raw$tests < 0) || any(rsv_raw$positive < 0) ||
    any(rsv_raw$positive > rsv_raw$tests))
  stop("Invalid RSV test counts, method, outlier indicators or HHS region")

rsv_raw <- dplyr::filter(rsv_raw,
                         week_start >= as.Date("2010-11-01"),
                         week_start <= END + 30)
if (anyDuplicated(paste(rsv_raw$week_start, rsv_raw$region, rsv_raw$method)))
  stop("Duplicate RSV (week, region, test method) records")

# Monthly test-weighted positivity. The original pilot assumes that each
# weekly total is distributed uniformly over seven days, including weeks
# spanning two calendar months. Reuse THAT convention for comparability.
weekly_to_monthly <- function(rows, label) {
  w <- rows |>
    dplyr::group_by(week_start) |>
    dplyr::summarise(tests = sum(tests), positive = sum(positive), .groups = "drop")
  if (!nrow(w) || any(w$tests <= 0) || any(w$positive > w$tests))
    stop(label, " has no usable data or contains invalid weekly totals")
  daily <- dplyr::bind_rows(lapply(0:6, function(k) {
    dplyr::transmute(w, date = week_start + k,
                     tests = tests / 7, positive = positive / 7)
  }))
  m <- daily |>
    dplyr::mutate(month = as.Date(format(date, "%Y-%m-01"))) |>
    dplyr::group_by(month) |>
    dplyr::summarise(days_covered = dplyr::n_distinct(date),
                     tests = sum(tests), positive = sum(positive), .groups = "drop") |>
    dplyr::mutate(coverage = days_covered / lubridate::days_in_month(month),
                  positivity = 100 * positive / tests) |>
    dplyr::filter(month >= FIRST_REQUIRED, month <= END)
  if (!identical(sort(m$month), seq(FIRST_REQUIRED, END, by = "month")) ||
      anyNA(m$positivity) || any(m$coverage < 1) || any(m$tests <= 0))
    stop(label, " RSV data do not cover every day of every required month")
  m
}

rsv_variants <- list(
  pooled = rsv_raw,
  pcr = dplyr::filter(rsv_raw, method == "pcr"),
  antigen = dplyr::filter(rsv_raw, method == "antigen"),
  pcr_no_outliers = dplyr::filter(rsv_raw, method == "pcr", flag == 0),
  pooled_no_outliers = dplyr::filter(rsv_raw, flag == 0)
)
monthly_rsv <- lapply(names(rsv_variants), function(name) {
  m <- weekly_to_monthly(rsv_variants[[name]], name)
  names(m)[names(m) == "positivity"] <- paste0("rsv_", name)
  names(m)[names(m) == "tests"] <- paste0("tests_", name)
  names(m)[names(m) == "positive"] <- paste0("pos_", name)
  dplyr::select(m, month, dplyr::all_of(c(paste0("rsv_", name),
                                        paste0("tests_", name),
                                        paste0("pos_", name))))
})
monthly_rsv <- Reduce(function(a, b) dplyr::full_join(a, b, by = "month"), monthly_rsv)

# Essential reproducibility check: reconstruction of the original pooled RSV
# positivity MUST match the original pilot's values, to numerical tolerance.
check <- dplyr::left_join(dplyr::select(pilot, month, pilot_rsv = rsv_pct),
                          dplyr::select(monthly_rsv, month, rsv_pooled), by = "month")
max_difference <- max(abs(check$pilot_rsv - check$rsv_pooled))
if (!is.finite(max_difference) || max_difference > 1e-5)
  stop("Pooled RSV reconstruction disagrees with pilot (max diff ",
       signif(max_difference, 4), " percentage points). Check weekly parsing.")
say("Reconstructed RSV matches original pilot; max difference: ",
    signif(max_difference, 3), " percentage points")

# Annual change in relative test methods, for assessing measurement drift.
rsv_method_share <- rsv_raw |>
  dplyr::mutate(year = as.integer(format(week_start + 6, "%Y"))) |>
  dplyr::filter(year >= 2011, year <= 2019) |>
  dplyr::group_by(year, method) |>
  dplyr::summarise(tests = sum(tests), positives = sum(positive),
                   flagged_rows = sum(flag == 1), .groups = "drop") |>
  dplyr::group_by(year) |>
  dplyr::mutate(share_of_tests_pct = 100 * tests / sum(tests),
                positivity = 100 * positives / tests) |>
  dplyr::ungroup()
readr::write_csv(rsv_method_share, file.path(OUT, "rsv_method_mix_by_year.csv"))

# C. HARMONIZE MONTHS AND BUILD PREDEFINED LAGGED EXPOSURES
flu_series <- dplyr::select(flu_month, month,
                             flu_pct = positivity, flu_coverage = coverage)
monthly <- dplyr::left_join(monthly_rsv, flu_series, by = "month") |>
  dplyr::arrange(month)
if (anyNA(monthly) || any(monthly$flu_coverage < 1))
  stop("Missing/incomplete influenza or RSV monthly data")
rsv_cols <- paste0("rsv_", names(rsv_variants))
if (any(vapply(monthly[, c("flu_pct", rsv_cols)],
               function(z) any(z < 0 | z > 100), logical(1))))
  stop("Invalid percentage positivity outside 0-100")

# One sensitivity model needs 2 previous MONTHS. November 2010 is outside
# the prepared influenza monthly file. Therefore ALL lag-comparison models
# are fitted on the SAME 107-month sample: February 2011-December 2019.
lag_specs <- c("lag0", "lag1", "lag2", "mean01", "mean012")
add_lags <- function(d, field, prefix) {
  x0 <- d[[field]]
  x1 <- dplyr::lag(x0, 1)
  x2 <- dplyr::lag(x0, 2)
  d[[paste0(prefix, "_lag0_5pp")]] <- x0 / 5
  d[[paste0(prefix, "_lag1_5pp")]] <- x1 / 5
  d[[paste0(prefix, "_lag2_5pp")]] <- x2 / 5
  d[[paste0(prefix, "_mean01_5pp")]] <- (x0 + x1) / 10
  d[[paste0(prefix, "_mean012_5pp")]] <- (x0 + x1 + x2) / 15
  d
}
monthly <- add_lags(monthly, "flu_pct", "flu")
for (variant in names(rsv_variants))
  monthly <- add_lags(monthly, paste0("rsv_", variant), paste0("rsv_", variant))

analysis <- dplyr::left_join(pilot[, c("month", "deaths")], monthly, by = "month") |>
  dplyr::arrange(month) |>
  dplyr::mutate(t = seq_len(dplyr::n()),
                days = lubridate::days_in_month(month),
                calendar_month = factor(as.integer(format(month, "%m")), levels = 1:12),
                sin1 = sin(2 * pi * t / 12), cos1 = cos(2 * pi * t / 12),
                sin2 = sin(4 * pi * t / 12), cos2 = cos(4 * pi * t / 12),
                flu_lab_era = as.integer(month >= as.Date("2015-10-01")))
readr::write_csv(analysis, file.path(OUT, "sensitivity_analysis_dataset.csv"), na = "")

# D. NEGATIVE BINOMIAL MODELS: SAME OUTCOME AND DEFAULT BASELINE AS PILOT
# Conditional on the same period/specification, smaller AIC is better.
# HAC/Wald intervals are reported for coefficients. LRT is NON-HAC and only
# descriptive of model improvement; it is not a causal test.
base_rhs <- function(seasonality, dat) {
  era <- if (length(unique(dat$flu_lab_era)) > 1) " + flu_lab_era" else ""
  switch(seasonality,
    pilot = paste0("splines::ns(t, df=4) + sin1 + cos1 + sin2 + cos2", era,
                   " + offset(log(days))"),
    month_fixed = paste0("splines::ns(t, df=3) + calendar_month", era,
                         " + offset(log(days))"),
    one_harmonic = paste0("splines::ns(t, df=4) + sin1 + cos1", era,
                          " + offset(log(days))"),
    no_lab_era = "splines::ns(t, df=4) + sin1 + cos1 + sin2 + cos2 + offset(log(days))",
    short_period = "splines::ns(t, df=2) + sin1 + cos1 + offset(log(days))",
    stop("Unknown seasonality: ", seasonality))
}
fit_nb <- function(formula, data) {
  warnings <- character()
  m <- withCallingHandlers(
    tryCatch(MASS::glm.nb(formula, data = data,
                          control = stats::glm.control(maxit = 100)),
             error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(model = m, warnings = paste(unique(warnings), collapse = " | "))
}

results <- list()
model_errors <- list()
add_model <- function(label, variant = "pcr", exposure = "joint",
                      lag = "mean01", season = "pilot", from = START,
                      to = END) {
  dat <- dplyr::filter(analysis, month >= from, month <= to)
  if (length(unique(dat$month)) < 36) stop("Too few months for: ", label)
  flu <- paste0("flu_", lag, "_5pp")
  rsv <- paste0("rsv_", variant, "_", lag, "_5pp")
  predictors <- switch(exposure,
    joint = c(flu, rsv), flu_only = flu, rsv_only = rsv,
    stop("Invalid exposure model"))
  if (anyNA(dat[, predictors, drop=FALSE])) {
    model_errors[[length(model_errors) + 1]] <<- data.frame(
      model = label, issue = "Missing lagged predictor(s)")
    return(invisible(NULL))
  }
  base_formula <- stats::as.formula(paste("deaths ~", base_rhs(season, dat)))
  full_formula <- stats::as.formula(paste("deaths ~", base_rhs(season, dat), "+",
                                          paste(predictors, collapse = " + ")))
  fit0 <- fit_nb(base_formula, dat)
  fit1 <- fit_nb(full_formula, dat)
  if (inherits(fit0$model, "error") || inherits(fit1$model, "error")) {
    problems <- character()
    if (inherits(fit0$model, "error")) problems <- c(problems, conditionMessage(fit0$model))
    if (inherits(fit1$model, "error")) problems <- c(problems, conditionMessage(fit1$model))
    model_errors[[length(model_errors) + 1]] <<- data.frame(
      model = label, issue = paste(problems, collapse=" | "))
    return(invisible(NULL))
  }
  null_model <- fit0$model; full_model <- fit1$model
  beta <- stats::coef(full_model)
  if (!all(predictors %in% names(beta)) || anyNA(beta[predictors])) {
    model_errors[[length(model_errors) + 1]] <<- data.frame(
      model = label, issue = "Missing/aliased predictor coefficient(s)")
    return(invisible(NULL))
  }
  V <- tryCatch(sandwich::NeweyWest(full_model, lag = 3,
                                     prewhite = FALSE, adjust = TRUE),
                error = function(e) e)
  if (inherits(V, "error")) {
    model_errors[[length(model_errors) + 1]] <<- data.frame(
      model = label, issue = paste("HAC covariance:", conditionMessage(V)))
    return(invisible(NULL))
  }
  hac_se <- sqrt(diag(V)[predictors])
  LR <- 2 * (as.numeric(stats::logLik(full_model)) -
             as.numeric(stats::logLik(null_model)))
  summary <- data.frame(model=label, variant=variant, exposure=exposure,
    lag=lag, season=season, start=as.character(min(dat$month)),
    end=as.character(max(dat$month)), n_months=nrow(dat),
    predictor=predictors, coefficient=unname(beta[predictors]),
    HAC_SE=unname(hac_se), IRR=exp(unname(beta[predictors])),
    CI_low=exp(unname(beta[predictors]-1.96*hac_se)),
    CI_high=exp(unname(beta[predictors]+1.96*hac_se)),
    p_HAC=2*stats::pnorm(-abs(unname(beta[predictors]/hac_se))),
    AIC_baseline=stats::AIC(null_model), AIC_extended=stats::AIC(full_model),
    delta_AIC=stats::AIC(full_model)-stats::AIC(null_model),
    LR=LR, LR_df=length(predictors),
    LR_p=stats::pchisq(max(0, LR), df=length(predictors),lower.tail=FALSE),
    Ljung_Box_p=stats::Box.test(stats::residuals(full_model,type="pearson"),
      lag=12, type="Ljung-Box")$p.value,
    warnings=paste(c(fit0$warnings,fit1$warnings)[nzchar(c(fit0$warnings,fit1$warnings))],
                   collapse=" | "), stringsAsFactors=FALSE)
  results[[length(results)+1L]] <<- summary
  invisible(summary)
}

say("Fitting primary PCR-only and test-method sensitivity models...")
add_model("pilot_reproduction_pooled_joint", "pooled", "joint")
add_model("pooled_RSV_only", "pooled", "rsv_only")
add_model("influenza_only", "pooled", "flu_only")
add_model("PCR_joint", "pcr", "joint")
add_model("PCR_RSV_only", "pcr", "rsv_only")
add_model("PCR_no_outlier_joint", "pcr_no_outliers", "joint")
add_model("PCR_no_outlier_RSV_only", "pcr_no_outliers", "rsv_only")
add_model("pooled_no_outlier_joint", "pooled_no_outliers", "joint")
add_model("antigen_RSV_only", "antigen", "rsv_only")

say("Fitting lag sensitivity models on IDENTICAL Feb 2011-Dec 2019 months...")
for (variant in c("pcr", "pooled", "pcr_no_outliers")) {
  for (lag in lag_specs) {
    for (exposure in c("rsv_only", "joint")) {
      label <- paste("lag", variant, lag, exposure, sep="_")
      add_model(label, variant, exposure, lag, "pilot",
                from=as.Date("2011-02-01"))
    }
  }
}

say("Fitting alternative seasonality models...")
for (season in c("one_harmonic", "month_fixed", "no_lab_era")) {
  for (variant in c("pcr", "pcr_no_outliers")) {
    for (exposure in c("rsv_only", "joint")) {
      add_model(paste("season", season, variant, exposure, sep="_"),
                variant, exposure, "mean01", season)
    }
  }
}

say("Fitting restricted-calendar-period models...")
for (period in c("2013plus", "2016plus")) {
  from <- if (period == "2013plus") as.Date("2013-01-01") else as.Date("2016-01-01")
  season <- if (period == "2013plus") "pilot" else "short_period"
  for (variant in c("pcr", "pcr_no_outliers")) {
    for (exposure in c("rsv_only", "joint")) {
      add_model(paste(period, variant, exposure, sep="_"),
                variant, exposure, "mean01", season, from=from)
    }
  }
}

if (!length(results)) stop("All models failed; examine the available input files")
all_effects <- dplyr::bind_rows(results)
core_labels <- c("pilot_reproduction_pooled_joint", "pooled_RSV_only",
  "influenza_only", "PCR_joint", "PCR_RSV_only",
  "PCR_no_outlier_joint", "PCR_no_outlier_RSV_only",
  "pooled_no_outlier_joint", "antigen_RSV_only")
core <- dplyr::filter(all_effects, model %in% core_labels)
lags <- dplyr::filter(all_effects, grepl("^lag_", model))
alternatives <- dplyr::filter(all_effects,
              grepl("^(season_|2013plus_|2016plus_)", model))
readr::write_csv(core, file.path(OUT, "core_sensitivity_effects.csv"))
readr::write_csv(lags, file.path(OUT, "lag_sensitivity_effects.csv"))
readr::write_csv(alternatives, file.path(OUT, "season_period_sensitivity_effects.csv"))
readr::write_csv(all_effects, file.path(OUT, "all_model_effects.csv"))
if (length(model_errors)) {
  readr::write_csv(dplyr::bind_rows(model_errors), file.path(OUT, "model_errors.csv"))
}

# E. FIGURES: monthly test-method drift and coefficient robustness.
plot_data <- dplyr::select(monthly, month, rsv_pooled, rsv_pcr,
                            rsv_pcr_no_outliers) |>
  tidyr::pivot_longer(-month, names_to="series", values_to="positivity")
p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=month,y=positivity,colour=series)) +
  ggplot2::geom_line(linewidth=0.55) +
  ggplot2::scale_colour_manual(values=c("rsv_pooled"="#555555",
                                         "rsv_pcr"="#176AA4",
                                         "rsv_pcr_no_outliers"="#D07628"),
    breaks=c("rsv_pooled","rsv_pcr","rsv_pcr_no_outliers"),
    labels=c("PCR + antigen", "PCR only", "PCR excluding flagged outliers")) +
  ggplot2::labs(x="Month",y="RSV tests positive (%)",colour=NULL,
                title="RSV surveillance measures used in sensitivity analysis") +
  ggplot2::theme_classic(base_size=11) +
  ggplot2::theme(legend.position="bottom")
ggplot2::ggsave(file.path(OUT,"figure_RSV_test_method_comparison.png"),
                p1,width=10,height=4.5,dpi=300)

share <- dplyr::filter(rsv_method_share, method == "pcr")
p2 <- ggplot2::ggplot(share, ggplot2::aes(year,share_of_tests_pct)) +
  ggplot2::geom_line(linewidth=0.8) + ggplot2::geom_point(size=1.8) +
  ggplot2::scale_y_continuous(limits=c(0,100)) +
  ggplot2::scale_x_continuous(breaks=2011:2019) +
  ggplot2::labs(x="Calendar year", y="PCR proportion of RSV tests (%)",
                title="Changing composition of RSV testing, 2011-2019") +
  ggplot2::theme_classic(base_size=11)
ggplot2::ggsave(file.path(OUT,"figure_RSV_PCR_share.png"),
                p2,width=8,height=4,dpi=300)

forest <- dplyr::filter(core, grepl("^rsv_", predictor)) |>
  dplyr::mutate(model = factor(model, levels=rev(unique(model))))
p3 <- ggplot2::ggplot(forest, ggplot2::aes(x=IRR,y=model)) +
  ggplot2::geom_vline(xintercept=1,linetype=2,colour="grey60") +
  ggplot2::geom_segment(ggplot2::aes(x=CI_low,xend=CI_high,y=model,yend=model),
                        linewidth=0.65) +
  ggplot2::geom_point(size=2) +
  ggplot2::labs(x="Incident rate ratio per +5 percentage points RSV positivity",
                y=NULL,title="Primary RSV sensitivity estimates (95% HAC CI)") +
  ggplot2::theme_classic(base_size=10)
ggplot2::ggsave(file.path(OUT,"figure_RSV_sensitivity_forest.png"),
                p3,width=10,height=5,dpi=300)

# F. HUMAN-READABLE BRIEF. Do not identify one nominal p < .05 from dozens
# of comparisons as definitive: all follow-up analyses are exploratory.
reproduce <- dplyr::filter(core, model == "pilot_reproduction_pooled_joint")
pcr_joint <- dplyr::filter(core, model == "PCR_joint", grepl("^rsv_", predictor))
pcr_only <- dplyr::filter(core, model == "PCR_RSV_only")
pcr_no_out <- dplyr::filter(core, model == "PCR_no_outlier_joint", grepl("^rsv_", predictor))
format_est <- function(x) {
  if (nrow(x) != 1) return("Not available (see model_errors.csv)")
  sprintf("IRR %.4f (95%% CI %.4f-%.4f); HAC p=%.4f; delta AIC=%.3f; LRT p=%.4f",
          x$IRR,x$CI_low,x$CI_high,x$p_HAC,x$delta_AIC,x$LR_p)
}
pcr_shares <- dplyr::filter(share, year %in% c(2011,2019))
notes <- c(
  "ALS/MND respiratory-virus follow-up sensitivity analysis (exploratory)",
  "---------------------------------------------------------------------",
  "Underlying-cause mortality G12.2, national USA, 2011-2019.",
  "Each exposure effect is per 5 percentage-point increase in test positivity.",
  "Primary exposure: mean positivity in index and preceding month.",
  "Primary model: same 4-df spline, two Fourier harmonics, influenza",
  "laboratory-era indicator and days-per-month offset as first pilot.",
  "HAC covariance: Newey-West lag 3. Full model tests use conventional LRT.",
  paste("Pooled RSV reconstruction max difference:", signif(max_difference,4), "percentage points"),
  paste("RSV rows flagged by CDC across study input:", sum(rsv_raw$flag==1)),
  paste("PCR share of RSV tests, 2011/2019:", paste(sprintf("%s: %.1f%%",
    pcr_shares$year,pcr_shares$share_of_tests_pct),collapse="; ")),
  "", "Original pilot model reproduced (two predictors):")
for (i in seq_len(nrow(reproduce))) notes <- c(notes,
    paste0("  ",reproduce$predictor[i],": ",format_est(reproduce[i,,drop=FALSE])))
notes <- c(notes,
  "", paste("PCR-only (joint with influenza):", format_est(pcr_joint)),
  paste("PCR-only (RSV alone):", format_est(pcr_only)),
  paste("PCR-only excluding flagged outliers (joint):", format_est(pcr_no_out)),
  "", "Interpretation guidance:",
  "- Compare IRR magnitude/direction/uncertainty across reasonable definitions.",
  "- Changes in RSV testing practices can bias surveillance positivity over time.",
  "- The many lag/seasonality/period comparisons are EXPLORATORY, not independent", 
  "  tests that justify selecting the smallest p-value.",
  "- No individual ALS cases or RSV-attributable deaths can be identified here.",
  "- The smaller 2016-2019 analysis has lower precision and a simpler baseline.",
  "- A stable signal is necessary but not sufficient to justify a full study.",
  "",paste("Successful fitted model specifications:",length(results)),
  paste("Failed model specifications:",length(model_errors)),
  "", "Full results in CSV tables; figures in sensitivity_results/."
)
writeLines(notes, file.path(OUT,"sensitivity_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"R_session_info.txt"))
say("Done. Results saved in: ", normalizePath(OUT))
cat(paste(notes, collapse="\n"),"\n")
