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
# ALS/MND and respiratory-virus mortality: feasibility pilot (2011-2019)
# ----------------------------------------------------------------------
# Study endpoint: US deaths with ICD-10 G12.2 as the UNDERLYING cause.
# G12.2 denotes motor neuron disease, not exclusively ALS.
# Predictors: national influenza and RSV laboratory percent positivity.
# This is an ecological association analysis, NOT a causal risk estimate.
#
# Run in RStudio with the working directory set to your project folder:
#   source("01_ALS_virus_mortality_pilot.R")
# Offline use: place the 4 official exports in raw_data/. See README.md.
# ----------------------------------------------------------------------

options(timeout = max(180, getOption("timeout")))
START_DATE <- as.Date("2011-01-01")
END_DATE   <- as.Date("2019-12-01")
RAW <- "raw_data"
OUT <- "pilot_results"
dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

needed <- c("dplyr", "tidyr", "readr", "lubridate", "ggplot2", "MASS",
            "splines", "sandwich", "remotes")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  if (!requireNamespace(p, quietly = TRUE)) stop("Cannot install required package: ", p)
}

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")
clean_names <- function(x) {
  names(x) <- gsub("(^_|_$)", "", gsub("_+", "_",
    gsub("[^a-z0-9]+", "_", tolower(names(x)))))
  x
}
num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))
first_col <- function(df, candidates, label) {
  z <- intersect(candidates, names(df))
  if (!length(z)) stop("Missing ", label, "; columns were: ", paste(names(df), collapse = ", "))
  df[[z[[1]]]]
}

# Read CDC date formats including 2012/04, 2012-04, April 2012 and Apr, 2012.
parse_month <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(as.Date(NA), length(x))
  for (i in seq_along(x)) {
    s <- x[i]
    if (is.na(s) || !nzchar(s)) next
    r <- regmatches(s, regexec("(19|20)[0-9]{2}[/.-](0?[1-9]|1[012])", s))[[1]]
    if (length(r)) {
      bits <- strsplit(gsub("[./]", "-", r[1]), "-", fixed = TRUE)[[1]]
      out[i] <- as.Date(sprintf("%s-%02d-01", bits[1], as.integer(bits[2])))
    } else {
      d <- suppressWarnings(lubridate::parse_date_time(s,
        orders = c("b Y", "B Y", "Y b", "Y B", "m Y", "Y m", "Y-m-d", "m/d/Y"),
        quiet = TRUE))
      if (!is.na(d)) out[i] <- as.Date(format(d, "%Y-%m-01"))
    }
  }
  out
}

parse_mortality <- function(x) {
  x <- clean_names(x)
  # If both label and code are present, prefer the machine-readable code.
  date_raw <- first_col(x,
    c("year_month_code", "year_and_month_code", "year_month", "year_and_month",
      "month_code", "month", "date"), "Year/Month")
  death_raw <- first_col(x, c("deaths", "death_count", "number_of_deaths"), "Deaths")
  d <- parse_month(date_raw)
  y <- num(death_raw)
  ok <- !is.na(d) & d >= START_DATE & d <= END_DATE
  z <- data.frame(month = d[ok], deaths = y[ok])
  if (anyDuplicated(z$month)) stop("Duplicate month in mortality file/API response; check grouping.")
  expected <- seq(START_DATE, END_DATE, by = "month")
  if (!identical(sort(z$month), expected)) {
    stop("Mortality missing months. Expected 108 months (2011-2019); found ", nrow(z),
         ". Re-export Year/Month, ALL US, G12.2 only, no additional grouping.")
  }
  if (anyNA(z$deaths) || any(z$deaths < 0) || any(z$deaths > 5000)) {
    stop("Mortality deaths invalid/suppressed or implausibly high. Confirm underlying cause G12.2 ONLY.")
  }
  z[order(z$month), ]
}

# 1) MORTALITY: prefer a previously verified manual export; otherwise try CDC WONDER API.
# The WONDER API can intermittently fail. We never fabricate or substitute counts.
mortality_file <- file.path(RAW, "wonder_G12_2_monthly.txt")
mortality_csv <- file.path(RAW, "wonder_G12_2_monthly.csv")
read_wonder_export <- function(path) {
  lines <- readr::read_lines(path, progress = FALSE)
  header <- grep("(?i)(year.?/?month|year and month|month).*(deaths)", lines, perl = TRUE)
  if (!length(header)) stop("Cannot locate Year/Month and Deaths header in ", path)
  # CDC WONDER exports append a blank line, annual totals, then narrative footnotes.
  # Stop at the first blank line *after* the contiguous tab-delimited table.
  h <- header[1]
  after <- if (h < length(lines)) seq.int(h + 1L, length(lines)) else integer()
  blanks <- after[!nzchar(trimws(lines[after]))]
  last <- if (length(blanks)) blanks[1] - 1L else length(lines)
  if (last <= h) stop("No mortality table rows found in ", path)
  content <- lines[h:last]
  tab <- readr::read_tsv(I(paste(content, collapse = "\n")),
     show_col_types = FALSE, progress = FALSE, col_types = readr::cols(.default = "c"))
  parse_mortality(tab)
}

get_mortality <- function() {
  if (file.exists(mortality_file)) return(read_wonder_export(mortality_file))
  if (file.exists(mortality_csv)) {
    x <- readr::read_csv(mortality_csv, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
    return(parse_mortality(x))
  }
  say("Trying automatic national-level CDC WONDER mortality request...")
  try_auto <- tryCatch({
    if (!requireNamespace("wonderapi", quietly = TRUE)) {
      remotes::install_github("socdataR/wonderapi", upgrade = "never", quiet = TRUE,
                              dependencies = NA, build_vignettes = FALSE)
    }
    # B_1 D76.V1-level2 groups by year/month; D76.V2 is underlying-cause ICD-10.
    x <- wonderapi::getData("Detailed Mortality", list(
      list("B_1", "D76.V1-level2"),
      list("V_D76.V2", "G12.2"),
      list("O_ucd", "D76.V2"),
      list("V_D76.V5", "*All*")
    ))
    z <- parse_mortality(x)
    readr::write_csv(z, file.path(RAW, "mortality_API_verified_2011_2019.csv"))
    z
  }, error = function(e) {
    say("Automated WONDER request unavailable: ", conditionMessage(e))
    NULL
  })
  if (is.null(try_auto)) stop(paste0(
    "CDC WONDER automated retrieval did not work. Export a national Year/Month query ",
    "for 2011-2019, UNDERLYING cause G12.2, and save Export Results as ",
    mortality_file, ". See README.md. Then rerun; influenza/RSV still download automatically."))
  try_auto
}
mortality <- get_mortality()
readr::write_csv(mortality, file.path(OUT, "mortality_monthly.csv"))

# 2) RSV: official CDC Socrata CSV (July 2010-June 2020).
rsv_path <- file.path(RAW, "cdc_nrevss_rsv_52kb-ccu2.csv")
if (!file.exists(rsv_path)) {
  say("Downloading CDC historical NREVSS RSV surveillance...")
  url <- "https://data.cdc.gov/resource/52kb-ccu2.csv?%24limit=50000"
  tryCatch(utils::download.file(url, rsv_path, mode = "wb", quiet = TRUE),
    error = function(e) stop("RSV download failed: ", conditionMessage(e),
                            ". Save dataset 52kb-ccu2 as ", rsv_path))
}
rsv_raw <- clean_names(readr::read_csv(rsv_path, show_col_types = FALSE,
                                        col_types = readr::cols(.default = "c")))
need_rsv <- c("repweekdate", "hhs_region", "testtype", "rsvpos", "rsvtest")
if (!all(need_rsv %in% names(rsv_raw))) stop("RSV schema changed: ", paste(names(rsv_raw), collapse=", "))
# Actual historical CDC export has repweekdate like 10JUL2010 (week ending
# Saturday); other downloads sometimes use ISO or US dates. Parse all three.
# Explicit English month lookup avoids locale-specific %b failures on Windows.
rsv_dates <- trimws(as.character(rsv_raw$repweekdate))
rsv_end <- suppressWarnings(as.Date(substr(rsv_dates, 1, 10), format="%Y-%m-%d"))
is_cdc_date <- is.na(rsv_end) & grepl("^[0-3][0-9][A-Za-z]{3}(19|20)[0-9]{2}$", rsv_dates)
if (any(is_cdc_date)) {
  v <- toupper(rsv_dates[is_cdc_date])
  mo <- match(substr(v, 3, 5), toupper(month.abb))
  if (anyNA(mo)) stop("Unrecognised month abbreviation in RSV file")
  rsv_end[is_cdc_date] <- as.Date(sprintf("%s-%02d-%s", substr(v, 6, 9), mo, substr(v, 1, 2)))
}
bad <- is.na(rsv_end)
if (any(bad)) rsv_end[bad] <- as.Date(lubridate::mdy(rsv_dates[bad]))
if (anyNA(rsv_end)) stop("RSV has invalid/unrecognised report dates: ",
                         paste(head(unique(rsv_dates[is.na(rsv_end)]), 5), collapse=", "))
rsv_raw$week_start <- rsv_end - 6
rsv_raw$tests <- num(rsv_raw$rsvtest)
rsv_raw$positive <- num(rsv_raw$rsvpos)
rsv_raw$region_n <- as.integer(num(rsv_raw$hhs_region))
rsv_raw$method <- gsub("[^a-z]", "", tolower(rsv_raw$testtype))
# Include all reported diagnostic types, weighted by tests. PCR-only is a
# planned sensitivity analysis; cross-year changes in test methods remain a limitation.
rsv_weekly <- rsv_raw |>
  dplyr::filter(.data$region_n %in% 1:10, .data$method %in% c("pcr", "antigen"),
                .data$week_start >= as.Date("2010-11-01"),
                .data$week_start <= as.Date("2019-12-31")) |>
  dplyr::group_by(.data$week_start) |>
  dplyr::summarise(tests = sum(.data$tests, na.rm=TRUE),
                   positive = sum(.data$positive, na.rm=TRUE),
                   bad_rows = sum(is.na(.data$tests) | is.na(.data$positive)),
                   .groups = "drop")
if (!nrow(rsv_weekly) || any(rsv_weekly$bad_rows > 0) ||
    any(rsv_weekly$tests <= 0) || any(rsv_weekly$positive > rsv_weekly$tests)) {
  stop("RSV invalid weekly surveillance data. Inspect input and CDC coding.")
}

# 3) INFLUENZA: national FluView WHO/NREVSS data via cdcfluview.
# Before 2015/16 data are combined; 2015/16 onward use clinical laboratories.
flu_cache <- file.path(RAW, "cdc_fluview_weekly_harmonized.csv")
normalize_flu <- function(x, source) {
  if (is.null(x) || nrow(x) == 0) return(data.frame())
  x <- clean_names(x)
  if (!"wk_date" %in% names(x)) {
    if (!all(c("year", "week") %in% names(x))) stop("FluView missing year/week/date columns")
    # CDC MMWR week starts on Sunday; week 1 contains January 4.
    yr <- as.integer(num(x$year)); wk <- as.integer(num(x$week))
    day4 <- as.Date(sprintf("%d-01-04", yr))
    x$wk_date <- lubridate::floor_date(day4, "week", week_start=7) + 7 * (wk - 1L)
  }
  date <- as.Date(x$wk_date)
  if (anyNA(date)) stop("FluView has unparseable dates")
  tests <- num(first_col(x, c("total_specimens", "totalspecimens"), "flu specimens"))
  if (all(c("total_a", "total_b") %in% names(x))) {
    pos <- num(x$total_a) + num(x$total_b)
  } else {
    pct <- num(first_col(x, c("percent_positive", "percentpositive"), "flu positivity"))
    pos <- tests * pct / 100
  }
  z <- data.frame(week_start=date, tests=tests, positive=pos, source=source)
  if (source == "pre_2015_combined") z <- subset(z, week_start < as.Date("2015-10-04"))
  if (source == "post_2015_clinical") z <- subset(z, week_start >= as.Date("2015-10-04"))
  z
}
read_flu_export <- function(path) {
  # Actual CDC FluView CSVs start with a quoted, single-cell explanatory note
  # containing the word "weekly"; matching "WEEK" anywhere in line 1
  # incorrectly treats that note as a header (and causes thousands of readr
  # parsing problems). Find the *field names* YEAR and WEEK instead.
  lines <- readLines(path, n=40, warn=FALSE)
  valid_header <- vapply(lines, function(line) {
    fields <- trimws(toupper(strsplit(gsub('"', '', line), ",", fixed=TRUE)[[1]]))
    all(c("YEAR", "WEEK", "TOTAL SPECIMENS") %in% fields)
  }, logical(1))
  h <- which(valid_header)
  if (!length(h)) stop("Cannot find FluView CSV column header in ",path,
                       ". File must have YEAR, WEEK, TOTAL SPECIMENS fields.")
  x <- readr::read_csv(path, skip=h[1]-1L, show_col_types=FALSE,
                       col_types=readr::cols(.default="c"), progress=FALSE)
  x <- clean_names(x)
  # In the original exports REGION TYPE is 'National' and REGION is 'X'.
  if ("region_type" %in% names(x)) {
    x <- x[tolower(trimws(x$region_type)) %in% c("national", "united states"), , drop=FALSE]
  } else if ("region" %in% names(x)) {
    x <- x[tolower(trimws(x$region)) %in% c("national", "united states"), , drop=FALSE]
  }
  if (!nrow(x)) stop("Manual FluView CSV has no national surveillance rows: ",path)
  if (anyDuplicated(paste(x$year, x$week, sep="/"))) {
    stop("Duplicate national YEAR/WEEK rows in FluView export: ",path)
  }
  x
}
get_flu <- function() {
  if (file.exists(flu_cache)) return(readr::read_csv(flu_cache, show_col_types=FALSE))
  pre_path <- file.path(RAW, "WHO_NREVSS_Combined_prior_to_2015_16.csv")
  post_path <- file.path(RAW, "WHO_NREVSS_Clinical_Labs.csv")
  if (file.exists(pre_path) && file.exists(post_path)) {
    say("Reading national FluView CSV exports supplied locally...")
    res <- list(combined_prior_to_2015_16=read_flu_export(pre_path),
                clinical_labs=read_flu_export(post_path))
  } else {
    say("Downloading national influenza surveillance via cdcfluview...")
    res <- tryCatch({
      if (!requireNamespace("cdcfluview", quietly=TRUE)) {
        install.packages("cdcfluview", repos="https://cloud.r-project.org")
      }
      if (!requireNamespace("cdcfluview", quietly=TRUE)) stop("cdcfluview unavailable")
      cdcfluview::who_nrevss("national", years=2010:2019)
    }, error=function(e) stop("FluView automatic download failed: ",conditionMessage(e),
         ". Place official National WHO_NREVSS_Combined_prior_to_2015_16.csv and ",
         "WHO_NREVSS_Clinical_Labs.csv in raw_data/; see README."))
  }
  pre <- normalize_flu(res$combined_prior_to_2015_16, "pre_2015_combined")
  post <- normalize_flu(res$clinical_labs, "post_2015_clinical")
  z <- dplyr::bind_rows(pre, post) |>
    dplyr::filter(.data$week_start >= as.Date("2010-11-01"),
                  .data$week_start <= as.Date("2019-12-31")) |>
    dplyr::arrange(.data$week_start)
  if (anyDuplicated(z$week_start) || anyNA(z$tests) || anyNA(z$positive) ||
      any(z$tests <= 0) || any(z$positive < 0) || any(z$positive > z$tests)) {
    stop("FluView dates are duplicated or lab counts are invalid; inspect raw FluView download")
  }
  readr::write_csv(z, flu_cache)
  z
}
flu_weekly <- get_flu()
flu_weekly$week_start <- as.Date(flu_weekly$week_start)
readr::write_csv(rsv_weekly, file.path(OUT,"rsv_weekly.csv"))
readr::write_csv(flu_weekly, file.path(OUT,"flu_weekly.csv"))

# Split weeks across their 7 days, then compute monthly test-weighted positivity.
# This avoids assigning entire boundary weeks to the wrong calendar month.
weekly_to_monthly <- function(df, name) {
  daily <- dplyr::bind_rows(lapply(0:6, function(k) {
    dplyr::transmute(df, date=.data$week_start+k,
                     tests=.data$tests/7, positive=.data$positive/7)
  })) |>
    dplyr::mutate(month=as.Date(format(.data$date, "%Y-%m-01")))
  d <- daily |>
    dplyr::group_by(.data$month) |>
    dplyr::summarise(days_covered=dplyr::n_distinct(.data$date),
                     tests=sum(.data$tests), positive=sum(.data$positive), .groups="drop") |>
    dplyr::mutate(days_in_month=lubridate::days_in_month(.data$month),
                  coverage=.data$days_covered / .data$days_in_month,
                  positivity=100 * .data$positive / .data$tests)
  report <- dplyr::filter(d, .data$month >= as.Date("2010-12-01"), .data$month <= END_DATE)
  readr::write_csv(report, file.path(OUT, paste0(name,"_monthly_coverage.csv")))
  expected <- seq(as.Date("2010-12-01"), END_DATE, by="month")
  if (!identical(sort(report$month), expected) || any(report$coverage < 1) ||
      anyNA(report$positivity) || any(report$tests <= 0)) {
    stop(name," weekly surveillance is missing whole/partial months during the pilot. ",
         "See ",file.path(OUT,paste0(name,"_monthly_coverage.csv")),
         ". No missing values have been replaced with zero.")
  }
  out <- report[,c("month","positivity")]
  names(out)[2] <- paste0(name,"_pct")
  out
}
flu_monthly <- weekly_to_monthly(flu_weekly,"flu")
rsv_monthly <- weekly_to_monthly(rsv_weekly,"rsv")

# 4) JOIN AND PREPARE LAGGED PREDICTORS
analysis_data <- dplyr::full_join(flu_monthly,rsv_monthly,by="month") |>
  dplyr::arrange(.data$month) |>
  dplyr::mutate(flu_lag1=dplyr::lag(.data$flu_pct),
                rsv_lag1=dplyr::lag(.data$rsv_pct),
                flu_5pp=(.data$flu_pct + .data$flu_lag1)/10,
                rsv_5pp=(.data$rsv_pct + .data$rsv_lag1)/10) |>
  dplyr::filter(.data$month >= START_DATE,.data$month <= END_DATE) |>
  dplyr::left_join(mortality,by="month") |>
  dplyr::mutate(t=seq_len(dplyr::n()),
                days=lubridate::days_in_month(.data$month),
                sin1=sin(2*pi*.data$t/12), cos1=cos(2*pi*.data$t/12),
                sin2=sin(4*pi*.data$t/12), cos2=cos(4*pi*.data$t/12),
                # Flu surveillance switched reporting conventions at 2015-16 season.
                flu_lab_era=as.integer(.data$month>=as.Date("2015-10-01")))
if (nrow(analysis_data) != 108 || anyNA(analysis_data)) stop("Pilot has missing/invalid monthly data")
readr::write_csv(analysis_data,file.path(OUT,"pilot_analysis_dataset.csv"))

# 5) PRIMARY PILOT MODELS: negative-binomial time-series regression.
# Offset month length. Time spline and two Fourier harmonics control secular
# trends/ordinary seasonality. Predictors = average of current and prior
# month's percent positivity, in 5 percentage-point units.
# NB addresses overdispersion; HAC inference is a preliminary serial-correlation
# safeguard. Check residual serial correlation before scientific interpretation.
base_formula <- deaths ~ splines::ns(t,df=4) + sin1+cos1+sin2+cos2 +
  flu_lab_era + offset(log(days))
full_formula <- update(base_formula, . ~ . + flu_5pp + rsv_5pp)
base <- MASS::glm.nb(base_formula, data=analysis_data, control=stats::glm.control(maxit=100))
full <- MASS::glm.nb(full_formula, data=analysis_data, control=stats::glm.control(maxit=100))

# Comparisons are exploratory. AIC alone does not establish a causal effect.
ll0 <- as.numeric(logLik(base)); ll1 <- as.numeric(logLik(full))
LR <- 2*(ll1-ll0)
LRT_p <- stats::pchisq(max(0,LR), df=2, lower.tail=FALSE)
model_compare <- data.frame(model=c("seasonal baseline","baseline + influenza + RSV"),
  AIC=c(AIC(base),AIC(full)), logLik=c(ll0,ll1), LRT_2df=c(NA,LR),
  LRT_p=c(NA,LRT_p))
readr::write_csv(model_compare,file.path(OUT,"model_comparison.csv"))

V <- sandwich::NeweyWest(full,lag=3,prewhite=FALSE,adjust=TRUE)
b <- stats::coef(full)
se <- sqrt(diag(V))
terms <- c("flu_5pp","rsv_5pp")
results <- data.frame(predictor=terms, log_IRR=unname(b[terms]), HAC_SE=unname(se[terms]))
results$IRR <- exp(results$log_IRR)
results$CI_low <- exp(results$log_IRR-1.96*results$HAC_SE)
results$CI_high <- exp(results$log_IRR+1.96*results$HAC_SE)
results$p_HAC <- 2*stats::pnorm(-abs(results$log_IRR/results$HAC_SE))
readr::write_csv(results,file.path(OUT,"effect_estimates_HAC.csv"))

analysis_data$expected_baseline <- as.numeric(predict(base, type="response"))
analysis_data$expected_full <- as.numeric(predict(full, type="response"))
readr::write_csv(analysis_data,file.path(OUT,"pilot_fitted_monthly.csv"))

# Diagnostics and compact visual outputs
pdf(file.path(OUT,"model_diagnostics.pdf"),width=9,height=7)
par(mfrow=c(2,2))
plot(residuals(full,type="pearson"),type="l",xlab="Month index",ylab="Pearson residual")
acf(residuals(full,type="pearson"), main="Residual autocorrelation")
plot(fitted(full),residuals(full,type="pearson"),xlab="Fitted deaths",ylab="Pearson residual")
abline(h=0,lty=2)
qqnorm(residuals(full,type="pearson")); qqline(residuals(full,type="pearson"))
dev.off()
box_p <- Box.test(residuals(full,type="pearson"),lag=12,type="Ljung-Box")$p.value

p1 <- ggplot2::ggplot(analysis_data,ggplot2::aes(x=month)) +
  ggplot2::geom_line(ggplot2::aes(y=deaths),linewidth=.55) +
  ggplot2::geom_line(ggplot2::aes(y=expected_baseline),linetype=2,linewidth=.75) +
  ggplot2::labs(x="Month",y="Deaths",title="Motor neuron disease deaths (G12.2)",
                subtitle="Solid: observed; dashed: expected under seasonal baseline") +
  ggplot2::theme_classic(base_size=12)
ggplot2::ggsave(file.path(OUT,"figure_1_mnd_deaths.png"),p1,width=10,height=4,dpi=300)

vlong <- analysis_data |>
  dplyr::select(.data$month,.data$flu_pct,.data$rsv_pct) |>
  tidyr::pivot_longer(cols=c("flu_pct","rsv_pct"),names_to="virus",values_to="percent_positive")
p2 <- ggplot2::ggplot(vlong,ggplot2::aes(month,percent_positive)) +
  ggplot2::geom_line(linewidth=.65) + ggplot2::facet_wrap(~virus,ncol=1,scales="free_y") +
  ggplot2::labs(x="Month",y="Tests positive (%)",title="US respiratory-virus surveillance") +
  ggplot2::theme_classic(base_size=12)
ggplot2::ggsave(file.path(OUT,"figure_2_virus_activity.png"),p2,width=10,height=6,dpi=300)

capture.output({
 cat("ALS/MND respiratory virus pilot 2011-2019\n")
 cat("Endpoint: National UCD G12.2 (motor neuron disease).\n")
 cat("Exposures: monthly surveillance positivity, current + one prior month mean.\n")
 cat("Flu series changes lab reporting at 2015/16. RSV includes PCR+antigen.\n")
 cat("Models adjusted for 4-df secular spline, 2 Fourier harmonics, lab era and days/month.\n\n")
 print(model_compare)
 cat("\nIncident rate ratios for each +5 percentage points in positivity:\n")
 print(results)
 cat("\nPearson residual Ljung-Box lag 12 p = ",box_p,"\n")
 cat("\nPearson correlation of virus predictors = ",
     cor(analysis_data$flu_5pp,analysis_data$rsv_5pp),"\n")
 cat("\nCAUTION: Ecological associations, no ALS prevalence denominator, no causality.\n")
 cat("Not publication-ready: investigate surveillance methods, seasonal controls,\n")
 cat("autocorrelation, residual diagnostics, comparator disease and sensitivity models.\n")
},file=file.path(OUT,"pilot_summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_session_info.txt"))
say("Done. Results in: ",normalizePath(OUT))
print(model_compare)
print(results)
