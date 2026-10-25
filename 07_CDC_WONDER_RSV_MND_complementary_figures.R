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
# =============================================================================
# RSV and motor neuron disease mortality analysis: complementary figures
# These figures complement the analyses in scripts 01-06.
# =============================================================================
# FIGURE 1: National monthly MND deaths and RSV PCR test positivity.
# FIGURE 2: Contrasting national and regional MND associations, 2011-2019
#           and 2013-2019; joint model in a separate panel.
# FIGURE 3: Technical sensitivity analysis and regional comparator outcomes.
#
# OUTPUT: Three high-resolution PNG figures, plotting_data_*.csv files, and
#         figure_data_manifest.csv.
#
# HOW TO USE:
#  1. Place this script in the analysis project folder containing the results
#     directories or CSV files (subfolders are searched).
#  2. Set the R working directory to that folder.
#  3. Install packages once: install.packages(c("ggplot2", "patchwork", "ragg"))
#  4. Run source("07_CDC_WONDER_RSV_MND_complementary_figures.R").
#
# If only final validation outputs are available, Figure 1 also needs the
# sensitivity_analysis_dataset.csv or mortality_monthly.csv input.
#
# SCIENTIFIC INTEGRITY / PROVENANCE:
# - Final model tables from script 06 are preferred when a validated estimate
#   can be extracted unambiguously. Extracted values are checked against the
#   rounded reference estimates documented below; discrepant values are not used.
# - Where a final source table is missing or ambiguous, Figure 2 uses clearly
#   identified rounded reference estimates. The output manifest records sources.
# - Figure 3A requires the Poisson HAC technical-validation CSV.
# - Figure 3B shows comparator confidence intervals only when their bounds can
#   be extracted from the final model CSV; missing intervals are not estimated.
# - All inferential claims remain exploratory and ecological, not causal.
# =============================================================================
# --------------------------- User settings ------------------------------------
DATA_ROOT <- "."  # folder to search recursively for original CSV outputs
OUTPUT_DIR <- file.path(DATA_ROOT, "complementary_figures")
WIDTH_IN <- 13
DPI <- 200        # 13 inches x 200 dpi = 2,600 px wide
BASE_FONT_PT <- 18  # legible when figures are viewed at reduced size
# Set FALSE to prohibit summary-value fallback in Figure 2. In that mode,
# unrecognised/missing final result tables cause an informative error.
ALLOW_REVIEWED_SUMMARY_FALLBACK <- TRUE
# Prefer the harmonised PCR national exposure from 01 if its exact variable
# can be identified; otherwise use the earlier, explicitly labelled PCR series.
PREFER_HARMONISED_FIGURE1_EXPOSURE <- TRUE

# ------------------------------ Packages --------------------------------------
needed <- c("ggplot2", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these R packages first: ",
                          paste(missing, collapse = ", "),
                          "\nRun: install.packages(c(\"ggplot2\", \"patchwork\", \"ragg\"))")
library(ggplot2)
library(patchwork)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ------------------------- File & column utilities ----------------------------
clean_names <- function(x) {
  names(x) <- make.unique(gsub("(^_|_$)", "",
      gsub("[^a-z0-9]+", "_", tolower(names(x)))), sep = "_")
  x
}
find_file <- function(filename, required = FALSE) {
  locations <- c(DATA_ROOT,
                 file.path(DATA_ROOT, "final_method_validation_results"),
                 file.path(DATA_ROOT, "technical_check_results"),
                 file.path(DATA_ROOT, "sensitivity_results"),
                 file.path(DATA_ROOT, "pilot_results"),
                 file.path(DATA_ROOT, "regional_validation_results"))
  candidates <- file.path(locations, filename)
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    all_files <- list.files(DATA_ROOT, recursive = TRUE,
                            full.names = TRUE, include.dirs = FALSE)
    hits <- all_files[basename(all_files) == filename]
  }
  if (!length(hits)) {
    if (required) stop("Cannot find required input: ", filename,
                       "\nSet DATA_ROOT to the directory with your CSV results.")
    return(NULL)
  }
  if (length(hits) > 1L) {
    warning("Several copies of ", filename, "; using: ", hits[1])
  }
  message("Reading ", filename, " from ", normalizePath(hits[1], mustWork = TRUE))
  hits[1]
}
read_input <- function(filename, required = FALSE) {
  p <- find_file(filename, required)
  if (is.null(p)) return(NULL)
  clean_names(read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("", "NA", "NaN", "NULL")))
}
first_col <- function(d, candidates) {
  hit <- candidates[candidates %in% names(d)]
  if (length(hit)) hit[1] else NULL
}
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
as_month <- function(d) {
  datecol <- first_col(d, c("month", "year_month", "date", "month_date",
                                "month_start", "yearmonth", "ym"))
  if (!is.null(datecol)) {
    x <- as.character(d[[datecol]])
    x <- trimws(x)
    x <- sub("^([0-9]{4})/([0-9]{1,2})$", "\\1-\\2", x)
    x <- ifelse(grepl("^[0-9]{4}-[0-9]{1,2}$", x), paste0(x, "-01"), x)
    x <- ifelse(grepl("^[0-9]{6}$", x),
                paste0(substr(x,1,4), "-", substr(x,5,6), "-01"), x)
    out <- as.Date(x)
    if (sum(!is.na(out)) > 0) return(as.Date(format(out, "%Y-%m-01")))
  }
  yy <- first_col(d, c("year", "calendar_year"))
  mm <- first_col(d, c("month_number", "month_num", "calendar_month", "month"))
  if (!is.null(yy) && !is.null(mm)) {
    m <- suppressWarnings(as.integer(as.character(d[[mm]])))
    if (all(is.na(m))) m <- match(tolower(substr(as.character(d[[mm]]), 1, 3)),
                                   tolower(month.abb))
    return(as.Date(sprintf("%04d-%02d-01", as.integer(d[[yy]]), m)))
  }
  stop("Could not identify a monthly date. CSV columns: ", paste(names(d), collapse=", "))
}
period_label <- function(z) {
  s <- as.character(z)
  ifelse(grepl("2013", s) & grepl("2019", s), "2013–2019",
         ifelse(grepl("2011", s) & grepl("2019", s), "2011–2019", NA_character_))
}

THEME_INK <- "#213149"
THEME_MUTED <- "#627287"
COL_NATIONAL <- "#0873A5"
COL_REGIONAL <- "#CD6736"
COL_MND <- "#175A8A"
COL_RSV <- "#C36B35"
COL_COMPARATOR <- "#446A7A"

figure_theme <- function() {
  theme_minimal(base_size = BASE_FONT_PT, base_family = "sans") +
    theme(
      text = element_text(colour = THEME_INK),
      plot.title = element_text(size = BASE_FONT_PT + 4, face = "bold",
                                margin = margin(b = 7)),
      plot.subtitle = element_text(size = BASE_FONT_PT - 1,
                                   colour = THEME_MUTED, margin = margin(b = 15)),
      plot.caption = element_text(size = BASE_FONT_PT - 5,
                                  colour = THEME_MUTED, hjust = 0,
                                  margin = margin(t = 16)),
      axis.title = element_text(size = BASE_FONT_PT),
      axis.text = element_text(size = BASE_FONT_PT - 2),
      axis.text.y = element_text(colour = THEME_INK),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E6EAF0", linewidth = 0.5),
      plot.margin = margin(18, 24, 20, 16),
      legend.position = "bottom",
      legend.text = element_text(size = BASE_FONT_PT - 1),
      strip.text = element_text(size = BASE_FONT_PT, face = "bold")
    )
}
save_figure <- function(plot, filename, height_in) {
  full <- file.path(OUTPUT_DIR, filename)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(full, plot = plot, width = WIDTH_IN, height = height_in,
           units = "in", dpi = DPI, bg = "white", device = ragg::agg_png)
  } else {
    # PNG is widely supported; ragg improves antialiasing.
    ggsave(full, plot = plot, width = WIDTH_IN, height = height_in,
           units = "in", dpi = DPI, bg = "white", device = "png")
  }
  message("Saved: ", normalizePath(full, mustWork = TRUE), " [",
          round(WIDTH_IN * DPI), " pixels wide]")
}
manifest <- data.frame(figure = character(), series = character(),
                       period = character(), source = character(),
                       details = character(), stringsAsFactors = FALSE)
record_source <- function(figure, series, period, source, details) {
  manifest <<- rbind(manifest,
                     data.frame(figure, series, period, source, details))
}

# ------------------------------ FIGURE 1 --------------------------------------
# Plot descriptive PCR positivity (not pooled antigen/PCR positivity).
sens <- read_input("sensitivity_analysis_dataset.csv")
mort <- read_input("mortality_monthly.csv")
if (is.null(mort) && !is.null(sens)) {
  mort <- sens
}
if (is.null(mort)) stop("Figure 1 needs mortality_monthly.csv or sensitivity_analysis_dataset.csv")
if (!"deaths" %in% names(mort)) stop("Mortality CSV has no 'deaths' column.")
mort_d <- data.frame(month = as_month(mort), deaths = num(mort$deaths))
mort_d <- mort_d[!is.na(mort_d$month) & !is.na(mort_d$deaths), ]
mort_d <- mort_d[!duplicated(mort_d$month), ]

exposure <- NULL
exposure_note <- NULL
if (PREFER_HARMONISED_FIGURE1_EXPOSURE) {
  nat <- read_input("01_national_harmonised_exposure.csv")
  if (!is.null(nat)) {
    # These are deliberately restricted to identifiable national PCR lag-01
    # variables: a generic column called 'exposure' is NOT assumed to be PCR.
    candidate <- first_col(nat, c(
      "national_rsv_pcr_mean01_5pp", "national_rsv_pcr_lag01_5pp",
      "national_pcr_mean01_5pp", "national_pcr_lag01_5pp",
      "national_pcr_lag01", "national_rsv_pcr_mean01",
      "national_rsv_mean01_5pp", "national_exposure_5pp",
      "national_rsv_pcr_percent", "national_pcr_positivity",
      "rsv_pcr_mean01_5pp", "rsv_pcr_lag01_5pp"))
    if (!is.null(candidate)) {
      d <- data.frame(month = as_month(nat), val = num(nat[[candidate]]))
      d <- d[!is.na(d$month) & !is.na(d$val), ]
      # Repeated months must have identical exposure values.
      if (any(duplicated(d$month))) {
        sums <- tapply(d$val, d$month, function(z) diff(range(z)) < 1e-8)
        if (!all(sums)) stop("The harmonised exposure CSV contains conflicting ",
                              "values for the same month.")
        d <- d[!duplicated(d$month), ]
      }
      # _5pp means the variable has been divided by five percentage points.
      if (grepl("5pp", candidate)) d$val <- d$val * 5
      if (!grepl("5pp", candidate) && max(d$val, na.rm=TRUE) <= 1)
        d$val <- d$val * 100
      if (any(d$val < 0 | d$val > 100)) {
        stop("Converted harmonised RSV positivity outside 0–100%. ",
             "Inspect field: ", candidate)
      }
      exposure <- data.frame(month = d$month, pct = d$val)
      exposure_note <- paste0("Harmonised national PCR lag-01 positivity (",
                              candidate, ")")
      record_source("1", "RSV positivity", "2011–2019",
                    "01_national_harmonised_exposure.csv", candidate)
    } else {
      warning("Cannot unambiguously identify harmonised national lag-01 PCR ",
              "column in file 01. Columns: ", paste(names(nat), collapse=", "),
              "\nFalling back to the correctly labelled ORIGINAL national PCR series.")
    }
  }
}
if (is.null(exposure)) {
  if (is.null(sens) || !"rsv_pcr" %in% names(sens)) {
    stop("Cannot create Figure 1: need 'rsv_pcr' in ",
         "sensitivity_analysis_dataset.csv, or a recognised PCR field in file 01.")
  }
  exposure <- data.frame(month = as_month(sens), pct = num(sens$rsv_pcr))
  exposure_note <- "Original national PCR-only RSV positivity (descriptive series)"
  record_source("1", "RSV positivity", "2011–2019",
                "sensitivity_analysis_dataset.csv", "rsv_pcr (%)")
}
record_source("1", "MND deaths", "2011–2019",
              if ("deaths" %in% names(mort) && !is.null(find_file("mortality_monthly.csv")))
                "mortality_monthly.csv" else "sensitivity_analysis_dataset.csv",
              "underlying-cause G12.2, national monthly deaths")
figure1_data <- merge(mort_d, exposure, by = "month", all = FALSE)
figure1_data <- figure1_data[figure1_data$month >= as.Date("2011-01-01") &
                             figure1_data$month <= as.Date("2019-12-01"), ]
figure1_data <- figure1_data[order(figure1_data$month), ]
if (nrow(figure1_data) < 80) stop("Figure 1 has unexpectedly sparse data: ",
                                  nrow(figure1_data), " months")
if (anyDuplicated(figure1_data$month)) stop("Duplicate months in Figure 1")
write.csv(figure1_data, file.path(OUTPUT_DIR, "plotting_data_figure1.csv"), row.names=FALSE)

p_deaths <- ggplot(figure1_data, aes(month, deaths)) +
  geom_line(linewidth = 1.1, colour = COL_MND, lineend = "round") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(title = "A  Monthly MND deaths", y = "Deaths (count)", x = NULL) +
  figure_theme() + theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
p_rsv <- ggplot(figure1_data, aes(month, pct)) +
  geom_line(linewidth = 1.1, colour = COL_RSV, lineend = "round") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = function(z) paste0(z, "%")) +
  labs(title = "B  RSV PCR test positivity", y = "PCR positivity", x = "Calendar year") +
  figure_theme()
fig1 <- (p_deaths / p_rsv) + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Respiratory virus circulation and MND mortality | United States",
    subtitle = "Monthly aggregate observations, 2011–2019",
    caption = paste0(
      "Sources: CDC WONDER (underlying-cause ICD-10 G12.2); CDC NREVSS. ",
      exposure_note, ".\n",
      "Descriptive time series only; the harmonised regional models used 107 complete calendar months."
    ),
    theme = theme(plot.title = element_text(size = BASE_FONT_PT + 5,
                                           face = "bold", colour = THEME_INK),
                  plot.subtitle = element_text(size = BASE_FONT_PT - 1,
                                               colour = THEME_MUTED),
                  plot.caption = element_text(size = BASE_FONT_PT - 6,
                                              colour = THEME_MUTED))
  )
save_figure(fig1, "Figure_1_MND_mortality_and_RSV_circulation.png", 9.2)

# ------------------------ Model-table extraction ------------------------------
# Figure 2's vetted reference data were transcribed from the final methodological
# validation summary. They are rounded to three decimals, not fabricated.
# More-precise original CSV estimates replace them ONLY if the matching row can
# be identified uniquely and checks pass.
reference <- data.frame(
  key = c("nat11","reg11","nat13","reg13","joint_nat11","joint_reg11"),
  period = c("2011–2019","2011–2019","2013–2019","2013–2019",
             "2011–2019","2011–2019"),
  model = c("National temporal", "Within-region", "National temporal",
            "Within-region", "Joint: national component", "Joint: regional deviation"),
  IRR = c(1.051, 0.970, 1.056, 0.962, 1.047, 0.973),
  low = c(1.029, 0.944, 1.026, 0.919, 1.018, 0.949),
  high = c(1.074, 0.998, 1.088, 1.008, 1.077, 0.997),
  source = rep("Reviewed results (rounded to 3 decimals)", 6),
  stringsAsFactors = FALSE
)

# Generic extraction: source table must have an unambiguous IRR and CI columns;
# the resulting row still has to match the validated prior summary.
parse_effects <- function(d) {
  if (is.null(d)) return(NULL)
  n <- nrow(d)
  irrcol <- first_col(d, c("irr", "incidence_rate_ratio", "rate_ratio",
                            "exp_beta", "exp_coef", "exp_estimate"))
  lowcol <- first_col(d, c("ci_low", "ci_lower", "ci_lo", "low_95", "lower_95",
                            "lower_95_ci", "ci95_low", "lcl", "conf_low",
                            "conf_lo", "conf_lower", "lower", "lo95"))
  highcol <- first_col(d, c("ci_high", "ci_upper", "ci_hi", "high_95",
                             "upper_95", "upper_95_ci", "ci95_high", "ucl",
                             "conf_high", "conf_hi", "conf_upper", "upper", "hi95"))
  bcol <- first_col(d, c("beta", "coefficient", "estimate", "coef"))
  secol <- first_col(d, c("se", "std_error", "standard_error", "hac_se", "cluster_se"))
  if (is.null(irrcol) && !is.null(bcol)) {
    vals <- num(d[[bcol]])
    # Only transform when a coefficient field is explicitly on log scale.
    if (bcol %in% c("beta", "coefficient", "coef")) {
      irr <- exp(vals)
    } else return(NULL)
  } else if (!is.null(irrcol)) {
    irr <- num(d[[irrcol]])
  } else return(NULL)
  low <- if (!is.null(lowcol)) num(d[[lowcol]]) else rep(NA_real_, n)
  high <- if (!is.null(highcol)) num(d[[highcol]]) else rep(NA_real_, n)
  if ((all(is.na(low)) || all(is.na(high))) &&
      !is.null(bcol) && !is.null(secol) &&
      bcol %in% c("beta", "coefficient", "coef")) {
    beta <- num(d[[bcol]]); se <- num(d[[secol]])
    low <- exp(beta - 1.96 * se)
    high <- exp(beta + 1.96 * se)
  }
  percol <- first_col(d, c("period", "time_period", "sample", "study_period",
                             "window", "years", "date_range"))
  if (is.null(percol)) return(NULL)
  outcome_col <- first_col(d, c("outcome", "cause", "cause_of_death",
                                  "death_cause", "outcome_name", "disease"))
  meta <- intersect(c("model", "analysis", "analysis_type", "model_type", "term",
                      "predictor", "component", "effect", "specification",
                      "contrast", "exposure", "estimator", "source", "scope",
                      "covariance", "se_type", "hac_lag", "hac_lags", "hac_lag_months",
                      "se_lag", "lag_months"), names(d))
  labels <- if (length(meta)) apply(d[meta], 1, function(z)
    paste(tolower(as.character(z)), collapse=" ")) else rep("", n)
  out <- data.frame(period=period_label(d[[percol]]),
                    outcome=if (!is.null(outcome_col))
                      tolower(as.character(d[[outcome_col]])) else rep("", n),
                    meta=labels, IRR=irr, low=low, high=high,
                    stringsAsFactors=FALSE)
  out <- out[is.finite(out$IRR) & out$IRR > 0 &
               !is.na(out$period), ]
  out
}
validate_row <- function(row, reference_irr, tolerance = 0.012) {
  if (nrow(row) != 1) return(FALSE)
  if (!is.finite(row$low) || !is.finite(row$high)) return(FALSE)
  if (!(row$low <= row$IRR && row$IRR <= row$high)) return(FALSE)
  abs(row$IRR - reference_irr) <= tolerance
}

# Directly read original final model tables. They have precedence over any
# figure-data export for the national and within-region comparisons.
raw4 <- read_input("04_harmonised_national_models.csv")
raw5 <- read_input("05_harmonised_regional_models.csv")
raw11 <- read_input("11_figure_data_harmonised_IRRs.csv")
p4 <- parse_effects(raw4)
p5 <- parse_effects(raw5)
p11 <- parse_effects(raw11)

# Select *only* a unique validated estimate for each expected key. If one
# table has several uncertainty methods, prefer metadata explicitly saying
# HAC(6) for national models, and 'within' for regional models.
update_from_table <- function(ref, key, tab, origin) {
  if (is.null(tab) || !nrow(tab)) return(ref)
  idx <- match(key, ref$key)
  goal <- ref[idx, ]
  candidates <- tab[tab$period == goal$period &
      (tab$outcome == "" | grepl("mnd|motor.neuron|g12|als", tab$outcome)), , drop=FALSE]
  if (nrow(candidates) == 0) return(ref)
  # Rows must be genuinely comparable to the validated reference estimate.
  candidates <- candidates[is.finite(candidates$low) &
                             is.finite(candidates$high) &
                             abs(candidates$IRR - goal$IRR) <= 0.012, , drop=FALSE]
  if (nrow(candidates) > 1 && key %in% c("nat11","nat13")) {
    ix <- grepl("hac.?6|hac.?lag.?6|lag.?6", candidates$meta)
    if (sum(ix) == 1) candidates <- candidates[ix, , drop=FALSE]
  }
  if (nrow(candidates) > 1 && key %in% c("reg11","reg13")) {
    ix <- grepl("within|national.year.month|year.month.fe|region.year", candidates$meta) &
          !grepl("joint", candidates$meta)
    if (sum(ix) == 1) candidates <- candidates[ix, , drop=FALSE]
  }
  # If duplicates are numerically identical they are harmless.
  if (nrow(candidates) > 1) {
    eq <- all(abs(candidates$IRR - candidates$IRR[1]) < 1e-10) &&
          all(abs(candidates$low - candidates$low[1]) < 1e-10) &&
          all(abs(candidates$high - candidates$high[1]) < 1e-10)
    if (eq) candidates <- candidates[1, , drop=FALSE]
  }
  if (validate_row(candidates, goal$IRR)) {
    ref[idx, c("IRR", "low", "high")] <- candidates[1, c("IRR", "low", "high")]
    ref$source[idx] <- origin
  }
  ref
}
for (key in c("nat11", "nat13")) {
  reference <- update_from_table(reference, key, p4, "04_harmonised_national_models.csv")
}
for (key in c("reg11", "reg13")) {
  reference <- update_from_table(reference, key, p5, "05_harmonised_regional_models.csv")
}
# Fallback to file 11 when the direct result tables were not parsed.
for (key in c("nat11", "reg11", "nat13", "reg13")) {
  idx <- match(key, reference$key)
  if (grepl("Reviewed", reference$source[idx])) {
    reference <- update_from_table(reference, key, p11,
                                    "11_figure_data_harmonised_IRRs.csv")
  }
}
# The joint model component rows are used ONLY if labelled adequately. Never
# conflate main national/within estimates with joint component estimates.
if (!is.null(p11)) {
  for (key in c("joint_nat11", "joint_reg11")) {
    idx <- match(key, reference$key)
    matches <- p11[p11$period == "2011–2019" &
                     (p11$outcome == "" | grepl("mnd|motor.neuron|g12",p11$outcome)) &
                     grepl("joint", p11$meta), , drop=FALSE]
    if (key == "joint_nat11") matches <- matches[grepl("national|common",matches$meta), , drop=FALSE]
    if (key == "joint_reg11") matches <- matches[grepl("regional|local|deviation",matches$meta), , drop=FALSE]
    matches <- matches[abs(matches$IRR - reference$IRR[idx]) <= 0.012, , drop=FALSE]
    if (validate_row(matches, reference$IRR[idx])) {
      reference[idx,c("IRR","low","high")] <- matches[1,c("IRR","low","high")]
      reference$source[idx] <- "11_figure_data_harmonised_IRRs.csv"
    }
  }
}
if (!ALLOW_REVIEWED_SUMMARY_FALLBACK && any(grepl("Reviewed", reference$source))) {
  stop("Could not unambiguously identify all Figure 2 estimates in final CSVs. ",
       "Inspect their headers or enable ALLOW_REVIEWED_SUMMARY_FALLBACK.")
}
for (i in seq_len(nrow(reference))) {
  record_source("2", reference$model[i], reference$period[i],
                reference$source[i], paste("IRR",round(reference$IRR[i],5)))
}
write.csv(reference, file.path(OUTPUT_DIR,"plotting_data_figure2.csv"), row.names=FALSE)

# ------------------------------ FIGURE 2 --------------------------------------
forest <- function(df, group_colours, x_limits = c(0.89, 1.16), x_step = 0.05,
                   heading = NULL, caption = NULL) {
  df$label <- factor(df$label, levels = rev(df$label))
  ggplot(df, aes(x=IRR, y=label, colour=group)) +
    geom_vline(xintercept=1, colour="#81909F", linewidth=0.7,
               linetype="dashed") +
    geom_errorbar(aes(xmin=low,xmax=high), orientation="y",
                  width=0.18, linewidth=1) +
    geom_point(size=4.0) +
    scale_color_manual(values=group_colours, guide="none") +
    scale_x_continuous(breaks=seq(x_limits[1], x_limits[2], by=x_step),
                       limits=x_limits) +
    labs(title=heading, x="Incidence rate ratio per 5-percentage-point increase in RSV positivity",
         y=NULL, caption=caption) +
    figure_theme()
}
main <- reference[reference$key %in% c("nat11","reg11","nat13","reg13"), ]
main <- main[match(c("nat11","reg11","nat13","reg13"),main$key), ]
main$group <- ifelse(main$model=="National temporal","national","regional")
main$label <- paste0(main$period, "   |   ", main$model)
p_2a <- forest(main, c(national=COL_NATIONAL,regional=COL_REGIONAL),
               x_limits=c(0.89,1.16),
               heading="A  National temporal versus within-region estimates")
joint <- reference[reference$key %in% c("joint_nat11","joint_reg11"), ]
joint <- joint[match(c("joint_nat11","joint_reg11"),joint$key), ]
joint$group <- ifelse(joint$key=="joint_nat11","national","regional")
joint$label <- c("National component", "Local-minus-national deviation")
p_2b <- forest(joint, c(national=COL_NATIONAL,regional=COL_REGIONAL),
               x_limits=c(0.89,1.16),
               heading="B  Joint model components (2011–2019)")
fig2 <- (p_2a / p_2b) + plot_layout(heights = c(1.12, 0.88)) +
  plot_annotation(
    title="Different geographical contrasts, different associations",
    subtitle="Final harmonised models | MND as underlying cause of death (G12.2)",
    caption=paste0(
      "Points and whiskers: model IRR and reported 95% CI. Dashed line: IRR = 1. ",
      "National and regional coefficients represent different estimands.\n",
      "Regional cluster-based inference is limited by 10 HHS regions; ",
      "2011–2019 MND wild-cluster score P = 0.0586. No causal interpretation."
    ),
    theme=theme(plot.title=element_text(size=BASE_FONT_PT+5,face="bold",colour=THEME_INK),
                plot.subtitle=element_text(size=BASE_FONT_PT-1,colour=THEME_MUTED),
                plot.caption=element_text(size=BASE_FONT_PT-6,colour=THEME_MUTED))
  )
save_figure(fig2, "Figure_2_national_vs_regional_IRRs.png", 9.6)

# ------------------------------ FIGURE 3A -------------------------------------
tech <- read_input("poisson_HAC_technical_results.csv", required=TRUE)
tech_required <- c("specification","irr","ci_low","ci_high","p_hac")
if (!all(tech_required %in% names(tech))) {
  stop("The technical model CSV has unexpected columns. Missing: ",
       paste(setdiff(tech_required,names(tech)),collapse=", "))
}
# These rows have been verified against the technical-validation summary.
tech_keys <- c("all_2011_2019__two_harmonics__rsv_only",
               "all_2011_2019__month_fixed__rsv_only",
               "all_2011_2019__one_harmonic__rsv_only",
               "later_2013_2019__two_harmonics__rsv_only",
               "all_2011_2019__exclude_outliers__rsv_only")
tech_labels <- c("2011–2019 | Two harmonics",
                 "2011–2019 | Calendar-month effects",
                 "2011–2019 | One harmonic",
                 "2013–2019 | Two harmonics",
                 "2011–2019 | Exclude flagged tests")
idx <- match(tech_keys, tech$specification)
if (anyNA(idx)) stop("Technical results missing required model specifications: ",
                       paste(tech_keys[is.na(idx)], collapse=", "))
tech_sel <- tech[idx,]
# Restrict plotted estimates to PCR-only unadjusted-for-influenza comparisons.
if ("rsv_measure" %in% names(tech_sel) && any(tech_sel$rsv_measure !=
                                             c("pcr","pcr","pcr","pcr",
                                               "pcr_excluding_outliers"))) {
  warning("Check RSV_measure metadata for the technical sensitivity figure.")
}
tech_df <- data.frame(label=tech_labels,
                      IRR=num(tech_sel$irr), low=num(tech_sel$ci_low),
                      high=num(tech_sel$ci_high),
                      p=num(tech_sel$p_hac),
                      group=rep("national",length(tech_labels)))
if (any(!is.finite(as.matrix(tech_df[c("IRR","low","high","p")]))))
  stop("Invalid technical-validation estimates; inspect CSV.")
write.csv(tech_df,file.path(OUTPUT_DIR,"plotting_data_figure3A.csv"),row.names=FALSE)
record_source("3A", "National model specification checks", "2011–2019/2013–2019",
              "poisson_HAC_technical_results.csv", "Five prespecified plotted rows")
p_3a <- forest(tech_df, c(national=COL_NATIONAL),
               x_limits=c(0.89,1.16),
               heading="A  Earlier national-model sensitivity checks")

# ------------------------------ FIGURE 3B -------------------------------------
# Comparator wild-score P values below come from the original 1024-pattern
# enumeration log (analysis_log(2).txt), not from conventional Poisson z-tests.
compare <- data.frame(
  outcome=c("MND","Parkinson's disease","Transport accidents",
            "Influenza/pneumonia"),
  irr=c(0.970,0.993,0.998,0.954),
  p_wild=c(0.058594,0.29688,0.72266,0.015625),
  low=rep(NA_real_,4), high=rep(NA_real_,4),
  source=rep("Reviewed regional results (rounded); CI not embedded",4),
  p_source=rep("Reviewed 1024-pattern enumeration log",4),
  stringsAsFactors=FALSE
)
# Use comparator CI from the final regional-model table, if it can be matched
# uniquely. Without actual CI fields, render a point estimate ONLY.
if (!is.null(p5) && nrow(p5)) {
  classify_outcome <- function(x) {
    ifelse(grepl("parkins|g20",x), "Parkinson's disease",
      ifelse(grepl("transport|accident",x), "Transport accidents",
        ifelse(grepl("pneumonia|respiratory|j09|j18|influenza",x),
               "Influenza/pneumonia",
          ifelse(grepl("mnd|motor.neuron|g12|als",x),"MND",NA_character_))))
  }
  p5$named_outcome <- classify_outcome(p5$outcome)
  for (i in seq_len(nrow(compare))) {
    rows <- p5[p5$period=="2011–2019" &
                 p5$named_outcome==compare$outcome[i] &
                 abs(p5$IRR-compare$irr[i]) <= 0.012 &
                 is.finite(p5$low) & is.finite(p5$high), ,drop=FALSE]
    if (nrow(rows) > 1) {
      ix <- grepl("within|national.year.month|year.month.fe",rows$meta) &
            !grepl("joint",rows$meta)
      if (sum(ix)==1) rows <- rows[ix,,drop=FALSE]
    }
    if (nrow(rows)==1 && rows$low <= rows$IRR && rows$IRR <= rows$high) {
      compare$irr[i] <- rows$IRR
      compare$low[i] <- rows$low
      compare$high[i] <- rows$high
      compare$source[i] <- "05_harmonised_regional_models.csv"
    }
  }
}
# Prefer P values from the actual 07 test output when its schema permits a
# unique matched lookup. Otherwise retain verified log values and say so.
raw7 <- read_input("07_exactly_enumerated_wild_cluster_score_tests.csv")
if (!is.null(raw7)) {
  test_period <- first_col(raw7, c("period", "study_period", "window", "sample"))
  test_outcome <- first_col(raw7, c("outcome", "cause", "outcome_name", "disease"))
  test_p <- first_col(raw7, c("p_value", "p_wild", "p_score", "p_cluster_score",
                              "p_wild_score", "p_enumerated", "p_two_sided",
                              "score_p", "wild_cluster_p", "wild_p"))
  if (!is.null(test_period) && !is.null(test_outcome) && !is.null(test_p)) {
    table_period <- period_label(raw7[[test_period]])
    table_outcome <- tolower(as.character(raw7[[test_outcome]]))
    source_names <- c("mnd", "parkinson", "transport", "respiratory")
    for (i in seq_len(nrow(compare))) {
      rows <- which(table_period == "2011–2019" &
        grepl(source_names[i], table_outcome))
      if (length(rows) == 1 && is.finite(num(raw7[[test_p]][rows])) &&
          abs(num(raw7[[test_p]][rows]) - compare$p_wild[i]) < 0.004) {
        compare$p_wild[i] <- num(raw7[[test_p]][rows])
        compare$p_source[i] <- "07_exactly_enumerated_wild_cluster_score_tests.csv"
      }
    }
  } else {
    warning("The score-test CSV columns were not recognised; using the reviewed",
            " P values from analysis_log(2).txt, without changing them.")
  }
}
compare$label <- compare$outcome
compare$IRR <- compare$irr
compare$group <- "regional"
compare$label <- factor(compare$label, levels=rev(compare$label))
compare$p_label <- paste0("P = ", formatC(compare$p_wild, digits=3,format="f"))
write.csv(compare, file.path(OUTPUT_DIR,"plotting_data_figure3B.csv"), row.names=FALSE)
for (i in seq_len(nrow(compare))) {
  record_source("3B", compare$outcome[i], "2011–2019", compare$source[i],
                paste0("score-test P source: ", compare$p_source[i],
                       "; P = ", compare$p_wild[i]))
}
p_3b <- ggplot(compare, aes(IRR,label)) +
  geom_vline(xintercept=1,colour="#81909F",linewidth=0.7,linetype="dashed") +
  geom_errorbar(data=compare[is.finite(compare$low) &
                               is.finite(compare$high), ],
                aes(xmin=low,xmax=high),orientation="y",width=0.18,
                linewidth=1.0,colour=COL_REGIONAL) +
  geom_point(size=4.0,colour=COL_REGIONAL) +
  geom_text(aes(x=1.105,label=p_label),colour=THEME_INK,
            size=5.3,hjust=0) +
  scale_x_continuous(breaks=seq(0.88,1.16,by=0.04),limits=c(0.88,1.17)) +
  coord_cartesian(clip="off") +
  labs(title="B  Within-region comparator outcomes (2011–2019)",
       x="Incidence rate ratio per 5-percentage-point increase in RSV positivity",
       y=NULL) +
  figure_theme()
fig3 <- (p_3a / p_3b) + plot_layout(heights=c(1.15,1)) +
  plot_annotation(
    title="Sensitivity analyses and disease comparators",
    subtitle="Exploratory analyses; models and inference methods differ between panels",
    caption=paste0(
      "Panel A: original national PCR-only sensitivity models, Poisson HAC(3) 95% CIs. ",
      "Panel B: harmonised regional IRRs and wild-cluster score P values (10 regions).\n",
      "Panel B whiskers are drawn ONLY if actual CI bounds can be read from the final ",
      "regional model table; no missing intervals are reconstructed."
    ),
    theme=theme(plot.title=element_text(size=BASE_FONT_PT+5,face="bold",colour=THEME_INK),
                plot.subtitle=element_text(size=BASE_FONT_PT-1,colour=THEME_MUTED),
                plot.caption=element_text(size=BASE_FONT_PT-6,colour=THEME_MUTED))
  )
save_figure(fig3, "Figure_3_sensitivity_and_comparator_outcomes.png", 10.5)

# ------------------------------- Audit trail ----------------------------------
write.csv(manifest, file.path(OUTPUT_DIR, "figure_data_manifest.csv"), row.names=FALSE)
writeLines(c(
  "RSV–MND complementary figures | provenance notes",
  "=================================================",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Figure 1 plots descriptive monthly counts and RSV PCR positivity.",
  "If the matched harmonised exposure table could not be parsed, the figure",
  "explicitly labels the original national PCR positivity series.",
  "",
  "Figure 2 data sources are documented per-row in plotting_data_figure2.csv.",
  "Rows marked 'Reviewed results' are rounded values from the final analysis",
  "and MUST be checked against the corresponding original tables before",
  "publishing them as exact effect estimates.",
  "",
  "Figure 3A uses the actual technical-validation model output.",
  "Figure 3B uses reviewed comparator IRRs and exact enumerated cluster-score",
  "P values; intervals are drawn only where verified original CSV CI fields",
  "were available.",
  "",
  "Exact enumeration of wild-cluster sign flips is NOT an exact finite-sample test.",
  "All analyses are exploratory, ecological, and not evidence of causality.",
  "",
  "These figures complement the analyses in scripts 01-06.",
  "Preserve the original aspect ratio if resizing.",
  "",
  "Recommended minimum checks before publishing:",
  "  1. Inspect figure_data_manifest.csv for summary-value fallbacks.",
  "  2. Compare Figure 2/3 plotted_data CSVs with files 04, 05, 07 and 11.",
  "  3. Confirm source and exposure definition written under Figure 1.",
  "  4. Update manuscript figure captions to match the final rendering."
), file.path(OUTPUT_DIR,"README_figures.txt"))
message("\nDONE — complementary figures and provenance tables saved in:\n",
        normalizePath(OUTPUT_DIR, mustWork=TRUE),"\n")
