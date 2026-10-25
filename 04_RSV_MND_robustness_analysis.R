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
# 04_RSV_MND_robustness_analysis.R
# RSV surveillance and US motor neuron disease mortality: comparator and
# regional robustness study, 2011-2019.
# v2: accepts complete ten-region HHS WONDER files in *_national_monthly.txt
# filenames, aggregates them to national counts only if all 1080 cells exist,
# and reuses the same originals for regional models.
# Plan drafted 2026-10-24, AFTER exploratory national pilot results.
# Do not call these analyses confirmatory or infer individual infection risk.
#
# Run from the SAME project folder as raw_data/:
#   install.packages(c("dplyr", "tidyr", "sandwich", "fixest", "ggplot2"))
#   source("04_RSV_MND_robustness_analysis.R")
#
# REQUIRED raw_data:
#   cdc_nrevss_rsv_52kb-ccu2.csv [already used in pilot]
#   Parkinson_G20_national_monthly.txt (national or COMPLETE HHS regional)
#   Transport_national_monthly.txt (national or COMPLETE HHS regional)
#   Respiratory_national_monthly.txt (national or COMPLETE HHS regional)
#   and original G12.2 national file or pilot_results/pilot_analysis_dataset.csv
# OPTIONAL (all four needed for panel stage):
#   MND_G12_2_HHS_monthly_2011_2019.txt
#   Parkinson_G20_HHS_monthly_2011_2019.txt
#   Transport_V01_V99_HHS_monthly_2011_2019.txt
#   Respiratory_J09_J18_HHS_monthly_2011_2019.txt
#
# A separate analysis_results/ directory is created; previous pilot results
# and raw_data/ are NEVER overwritten. No online data download is performed.

# ========== 0. Settings =======================================================
START <- as.Date("2011-01-01")
END <- as.Date("2019-12-01")
FIRST_RSV <- as.Date("2010-12-01")  # for December-2010 exposure lag
OUT <- "analysis_results"
B_CLUSTER <- 299L  # diagnostic pairs-bootstrap, NOT exact few-cluster inference
BOOT_SEED <- 20261024L
MIN_PANEL_COVERAGE <- 0.80
MIN_MONTHS_PER_REGION <- 60L

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(OUT, "analysis_log.txt")
if (file.exists(logfile)) file.remove(logfile)
log_msg <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n")
  cat(s, "\n", file = logfile, append = TRUE)
}
fail <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
write_out <- function(x, filename) {
  utils::write.csv(x, file.path(OUT, filename), na = "", row.names = FALSE)
}
need_pkgs <- c("dplyr", "tidyr", "sandwich", "fixest", "ggplot2")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) fail("Install missing packages: install.packages(c(",
                               paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
expected_months <- seq(START, END, by = "month")
get_days <- function(m) {
  y <- as.integer(format(m, "%Y")); mo <- as.integer(format(m, "%m"))
  ny <- y + as.integer(mo == 12L); nm <- mo %% 12L + 1L
  as.integer(as.Date(sprintf("%04d-%02d-01", ny, nm)) - m)
}

# ========== 1. CDC WONDER import with metadata and suppression checks =========
# Read exact tabular section. Never compute suppressed counts using residuals
# from national totals; this would undermine WONDER's disclosure controls.
NAMES <- c(mnd="MND", parkinson="Parkinson", transport="Transport", respiratory="Respiratory")
CODES <- c(mnd="G12.2", parkinson="G20", transport="V01-V99", respiratory="J09-J18")
raw_files <- list(
  national = c(
    mnd="raw_data/wonder_G12_2_monthly.txt",
    parkinson="raw_data/Parkinson_G20_national_monthly.txt",
    transport="raw_data/Transport_national_monthly.txt",
    respiratory="raw_data/Respiratory_national_monthly.txt"),
  regional = c(
    mnd="raw_data/MND_G12_2_HHS_monthly_2011_2019.txt",
    parkinson="raw_data/Parkinson_G20_HHS_monthly_2011_2019.txt",
    transport="raw_data/Transport_V01_V99_HHS_monthly_2011_2019.txt",
    respiratory="raw_data/Respiratory_J09_J18_HHS_monthly_2011_2019.txt")
)

check_wonder_metadata <- function(lines, expected_code, filename) {
  meta <- lines[grepl("^\\\"?ICD-10 Codes:", trimws(lines), ignore.case = TRUE)]
  dataset <- lines[grepl("^\\\"?Dataset:", trimws(lines), ignore.case = TRUE)]
  if (!length(meta)) fail("Missing ICD-10 query metadata in ", filename,
                          ". Export the full CDC WONDER text table WITH footnotes.")
  met <- toupper(gsub("[–—]", "-", meta[[1]]))
  pattern <- switch(expected_code,
                    "G12.2" = "G12\\.2(?![0-9])",
                    "G20" = "G20(?![0-9])",
                    "V01-V99" = "V01\\s*[-:]\\s*V99",
                    "J09-J18" = "J09\\s*[-:]\\s*J18")
  if (!grepl(pattern, met, perl=TRUE))
    fail("Wrong/incomplete ICD selection in ", filename, ". Expected ",
         expected_code, "; export says: ", meta[[1]])
  if (!length(dataset) || !grepl("Underlying Cause of Death", dataset[[1]], ignore.case=TRUE))
    fail("Wrong CDC WONDER database (or missing Dataset metadata) in ", filename)
  if (!any(grepl("2011", lines, fixed=TRUE)) ||
      !any(grepl("2019", lines, fixed=TRUE)))
    fail("Cannot verify 2011-2019 coverage in ", filename)
}

read_wonder <- function(path, code, regional=FALSE) {
  if (!file.exists(path)) fail("Missing required file: ", path)
  lines <- readLines(path, warn=FALSE, encoding="UTF-8")
  check_wonder_metadata(lines, code, path)
  header_candidates <- which(grepl("Month Code", lines, fixed=TRUE) &
                             grepl("Deaths", lines, fixed=TRUE) &
                             grepl("\t", lines, fixed=TRUE))
  if (!length(header_candidates))
    fail("Missing WONDER data header (Month Code and Deaths): ", path)
  hdr <- header_candidates[[1]]
  tail_lines <- lines[hdr:length(lines)]
  end <- which(grepl('^\\s*"?---"?\\s*$', tail_lines) |
               grepl('^"?Dataset:', tail_lines) |
               grepl('^"?Query Parameters:', tail_lines))
  if (length(end)) tail_lines <- tail_lines[seq_len(end[[1]]-1L)]
  if (length(tail_lines) < 2) fail("No mortality data rows in ", path)
  tab <- tryCatch(utils::read.delim(text=paste(tail_lines, collapse="\n"),
                    sep="\t", quote='"', comment.char="", check.names=FALSE,
                    stringsAsFactors=FALSE, colClasses="character", fill=TRUE,
                    blank.lines.skip=TRUE), error=function(e) fail("Parse error in ", path,
                                                               ": ", conditionMessage(e)))
  required <- c("Month Code", "Deaths")
  if (!all(required %in% names(tab))) fail("Missing WONDER columns in ", path)
  tab <- tab[!is.na(tab[["Month Code"]]) &
             grepl("^20[0-9]{2}/(0[1-9]|1[012])$", tab[["Month Code"]]), , drop=FALSE]
  if (!nrow(tab)) fail("No year/month death rows found: ", path)
  tab$month <- as.Date(paste0(gsub("/", "-", tab[["Month Code"]]), "-01"))
  tab <- tab[tab$month >= START & tab$month <= END, , drop=FALSE]
  if (!nrow(tab)) fail("No observations for 2011-2019: ", path)
  original_deaths <- trimws(tab[["Deaths"]])
  tab$deaths <- suppressWarnings(as.numeric(gsub(",", "", original_deaths, fixed=TRUE)))
  if (any(!is.na(tab$deaths) &
          (tab$deaths < 0 | tab$deaths != round(tab$deaths))))
    fail("Invalid death counts: ", path)
  # Explicitly suppressed cells stay NA. Absent region/month rows ALSO stay NA.
  # Do NOT turn absent, blank, 'Suppressed', 'Unreliable' etc into 0.
  if (!regional) {
    if (anyDuplicated(tab$month)) {
      # These three CDC exports were actually GROUPED BY HHS REGION, YEAR,
      # MONTH, despite being named '*_national_monthly.txt'. This is usable:
      # aggregate ONLY fully observed, unrestricted 10-region months. NEVER
      # infer suppressed counts from totals or treat missing cells as zero.
      if (!"HHS Region Code" %in% names(tab))
        fail("Duplicate national month rows without HHS region codes in: ", path,
             ". Inspect the original WONDER grouping.")
      rr <- read_wonder(path, code, regional=TRUE)
      if (nrow(rr) != 1080L || anyNA(rr$deaths) ||
          !setequal(unique(rr$region), 1:10))
        fail("Cannot safely aggregate this HHS export to national counts: ",path,
             ". Expected 10 complete regions x 108 months, no suppression.")
      regional_n <- stats::aggregate(region ~ month, rr, FUN=length)
      if (!identical(sort(regional_n$month),expected_months) ||
          any(regional_n$region != 10L))
        fail("HHS region export missing region-months: ",path)
      national_counts <- stats::aggregate(deaths ~ month, rr, FUN=sum)
      national_counts <- national_counts[order(national_counts$month),,drop=FALSE]
      log_msg("National series derived from all 10 complete HHS regions: ", path,
              "; total deaths=",sum(national_counts$deaths),
              "; this is an aggregate, NOT a separate national WONDER export")
      return(national_counts)
    }
    if (!identical(sort(tab$month), expected_months))
      fail("National file must contain all 108 months, 2011-2019: ", path)
    if (anyNA(tab$deaths)) fail("National file has missing/suppressed deaths: ", path)
    return(data.frame(month=tab$month, deaths=tab$deaths))
  }
  reg_cols <- grep("^HHS Region Code$", names(tab), value=TRUE, ignore.case=TRUE)
  if (!length(reg_cols)) reg_cols <- grep("^HHS Region$", names(tab), value=TRUE,
                                         ignore.case=TRUE)
  if (!length(reg_cols))
    fail("Regional file is NOT grouped by HHS Region: ", path,
         ". Re-export with HHS Region AND Year/Month.")
  s <- trimws(tab[[reg_cols[[1]]]])
  # Supports '1', '01', 'Region 1', 'HHS Region 10', or labels containing 'Region 01'.
  m <- regexpr("(?<![0-9])(10|0?[1-9])(?![0-9])", s, perl=TRUE)
  tab$region <- as.integer(regmatches(s, m))
  if (anyNA(tab$region) || any(!tab$region %in% 1:10))
    fail("Cannot interpret the HHS region codes in ", path,
         ". Check that grouping is HHS, NOT Census region.")
  if (anyDuplicated(paste(tab$region, tab$month)))
    fail("Duplicate HHS region-month rows in ", path, ". Check grouping options.")
  if (any(tab$deaths >= 1 & tab$deaths <= 9, na.rm=TRUE))
    fail("Regional file contains explicit counts from 1 to 9: ",path,
         ". CDC WONDER restricted counts must not be released; check source.")
  data.frame(region=tab$region, month=tab$month, deaths=tab$deaths,
             suppression_label=ifelse(is.na(tab$deaths), original_deaths, ""))
}

# ========== 2. Rebuild PCR surveillance using original week/day allocation ====
rsv_path <- "raw_data/cdc_nrevss_rsv_52kb-ccu2.csv"
if (!file.exists(rsv_path)) fail("Missing original CDC RSV data: ", rsv_path)
rsv <- utils::read.csv(rsv_path, stringsAsFactors=FALSE, check.names=FALSE)
fields <- c("repweekdate","hhs_region","testtype","rsvpos","rsvtest","outlier")
if (!all(fields %in% names(rsv))) fail("RSV file lacks required columns: ",
                                     paste(setdiff(fields,names(rsv)),collapse=", "))
s <- toupper(trimws(as.character(rsv$repweekdate)))
if (any(!grepl("^[0-9]{2}[A-Z]{3}[0-9]{4}$", s))) fail("Invalid RSV week-ending date")
mo <- match(substr(s,3,5), toupper(month.abb))
if (anyNA(mo)) fail("Unknown month abbreviation in RSV week-ending date")
week_end <- as.Date(sprintf("%s-%02d-%s",substr(s,6,9),mo,substr(s,1,2)))
if (anyNA(week_end)) fail("Cannot parse RSV report weeks")
rsv$week_start <- week_end - 6L
rsv$region <- suppressWarnings(as.integer(rsv$hhs_region))
rsv$method <- tolower(trimws(as.character(rsv$testtype)))
rsv$tests <- suppressWarnings(as.numeric(rsv$rsvtest))
rsv$positives <- suppressWarnings(as.numeric(rsv$rsvpos))
rsv$flag <- suppressWarnings(as.integer(rsv$outlier))
if (anyNA(rsv[,c("region","method","tests","positives","flag")]) ||
    any(!rsv$region %in% 1:10) || any(!rsv$method %in% c("pcr","antigen")) ||
    any(!rsv$flag %in% 0:1) || any(rsv$tests < 0) ||
    any(rsv$positives < 0 | rsv$positives > rsv$tests))
  fail("Invalid RSV data types, counts, methods, outlier flags or regions")
if (anyDuplicated(paste(rsv$week_start,rsv$region,rsv$method)))
  fail("Duplicate RSV region-week-test_method records")

# Split every Sunday-Saturday weekly count equally among its seven calendar
# days, EXACTLY as in the earlier sensitivity analysis. Positivity is then
# based on SUM(positives)/SUM(tests), never an unweighted mean of percentages.
pcr <- rsv[rsv$method == "pcr", , drop=FALSE]
make_monthly_pcr <- function(z, name) {
  expanded <- dplyr::bind_rows(lapply(0:6, function(k) {
    data.frame(region=z$region, date=z$week_start+k,
               tests=z$tests/7, positives=z$positives/7)
  }))
  expanded$month <- as.Date(format(expanded$date,"%Y-%m-01"))
  expanded <- expanded[expanded$month>=FIRST_RSV & expanded$month<=END,,drop=FALSE]
  grouped <- expanded |>
    dplyr::group_by(region,month) |>
    dplyr::summarise(tests=sum(tests), positives=sum(positives),
                     days_covered=dplyr::n_distinct(date), .groups="drop")
  names(grouped)[names(grouped)=="tests"] <- paste0("tests_", name)
  names(grouped)[names(grouped)=="positives"] <- paste0("positives_", name)
  names(grouped)[names(grouped)=="days_covered"] <- paste0("coverage_days_", name)
  grouped[[paste0("pct_",name)]] <- ifelse(grouped[[paste0("tests_",name)]]>0,
              100*grouped[[paste0("positives_",name)]]/
                  grouped[[paste0("tests_",name)]], NA_real_)
  grouped
}
mon_pcr <- make_monthly_pcr(pcr, "pcr")
mon_clean <- make_monthly_pcr(pcr[pcr$flag==0,,drop=FALSE], "pcr_clean")
region_grid <- expand.grid(region=1:10,month=seq(FIRST_RSV,END,by="month"))
region_grid$month <- as.Date(region_grid$month, origin="1970-01-01")
reg_rsv <- region_grid |>
  dplyr::left_join(mon_pcr, by=c("region","month")) |>
  dplyr::left_join(mon_clean,by=c("region","month")) |>
  dplyr::arrange(region,month) |>
  dplyr::group_by(region) |>
  dplyr::mutate(rsv_mean01_5pp=(pct_pcr+dplyr::lag(pct_pcr))/10,
                rsv_clean_mean01_5pp=(pct_pcr_clean+dplyr::lag(pct_pcr_clean))/10,
                rsv_lag0_5pp=pct_pcr/5,
                rsv_lag1_5pp=dplyr::lag(pct_pcr)/5,
                rsv_lead2_5pp=dplyr::lead(pct_pcr,2)/5) |>
  dplyr::ungroup()
reg_rsv <- reg_rsv[reg_rsv$month>=START,,drop=FALSE]
# Reconstruct NATIONAL series by summing weighted counts across regions BEFORE
# calculating positivity; the lag is a mean of two MONTHLY percentages.
nt <- mon_pcr |>
  dplyr::group_by(month) |>
  dplyr::summarise(tests_pcr=sum(tests_pcr),positives_pcr=sum(positives_pcr),
                   .groups="drop")
nc <- mon_clean |>
  dplyr::group_by(month) |>
  dplyr::summarise(tests_clean=sum(tests_pcr_clean),
                   positives_clean=sum(positives_pcr_clean),.groups="drop")
nat_rsv <- dplyr::full_join(nt,nc,by="month") |>
  dplyr::arrange(month) |>
  dplyr::mutate(pct_pcr=100*positives_pcr/tests_pcr,
                pct_pcr_clean=100*positives_clean/tests_clean,
                rsv_mean01_5pp=(pct_pcr+dplyr::lag(pct_pcr))/10,
                rsv_clean_mean01_5pp=(pct_pcr_clean+dplyr::lag(pct_pcr_clean))/10)
nat_rsv <- nat_rsv[nat_rsv$month >= START & nat_rsv$month <= END,,drop=FALSE]
if (!identical(nat_rsv$month, expected_months) ||
    any(!is.finite(nat_rsv$rsv_mean01_5pp)))
  fail("National RSV PCR reconstruction has gaps or invalid exposure")
write_out(reg_rsv,"01_regional_rsv_monthly.csv")
write_out(nat_rsv,"02_national_rsv_monthly.csv")
method_annual <- rsv |>
  dplyr::mutate(year=as.integer(format(week_end,"%Y"))) |>
  dplyr::filter(year>=2011,year<=2019) |>
  dplyr::group_by(year,method) |>
  dplyr::summarise(tests=sum(tests),positives=sum(positives),
                   n_flagged=sum(flag==1),.groups="drop") |>
  dplyr::group_by(year) |>
  dplyr::mutate(share_of_all_tests=100*tests/sum(tests),
                positivity_pct=100*positives/tests) |>
  dplyr::ungroup()
write_out(method_annual,"03_RSV_test_method_drift.csv")

# Reconcile PCR reconstruction against the original exploratory CSV when found.
checkfiles <- c("sensitivity_results/sensitivity_analysis_dataset.csv",
                "input_data/sensitivity_analysis_dataset.csv",
                "sensitivity_analysis_dataset.csv")
checkfiles <- checkfiles[file.exists(checkfiles)]
if (length(checkfiles)) {
  old <- utils::read.csv(checkfiles[[1]],stringsAsFactors=FALSE)
  old$month <- as.Date(old$month)
  chk <- dplyr::inner_join(nat_rsv,old,by="month")
  for (v in c("rsv_pcr", "rsv_pcr_no_outliers")) {
    src <- if (v=="rsv_pcr") "pct_pcr" else "pct_pcr_clean"
    if (v %in% names(chk)) {
      delta <- max(abs(chk[[src]]-chk[[v]]))
      if (!is.finite(delta) || delta>1e-5)
        fail("RSV reconstruction DOES NOT match prior sensitivity dataset: ",
             v,"; max difference=",delta," percentage points")
      log_msg("Reconstructed RSV ",v," agrees with pilot; max diff=",signif(delta,3))
    }
  }
}

# ========== 3. National mortality import =====================================
# Fallback to original verified pilot MND file when G12 national WONDER not found.
get_mnd_national <- function() {
  if (file.exists(raw_files$national[["mnd"]]))
    return(read_wonder(raw_files$national[["mnd"]], CODES[["mnd"]]))
  oldfiles <- c("pilot_results/pilot_analysis_dataset.csv",
                "sensitivity_results/sensitivity_analysis_dataset.csv",
                "input_data/sensitivity_analysis_dataset.csv")
  oldfiles <- oldfiles[file.exists(oldfiles)]
  if (!length(oldfiles)) fail("Cannot find old national G12.2 mortality: provide ",
                              raw_files$national[["mnd"]],
                              " OR pilot_results/pilot_analysis_dataset.csv")
  d <- utils::read.csv(oldfiles[[1]],stringsAsFactors=FALSE)
  if (!all(c("month","deaths") %in% names(d))) fail("Invalid original MND pilot CSV")
  d$month <- as.Date(d$month)
  if (nrow(d)!=108L || !identical(d$month,expected_months) ||
      anyNA(d$deaths) || any(d$deaths<0)) fail("Invalid national MND pilot months")
  d[,c("month","deaths")]
}

nat_deaths <- list(mnd=get_mnd_national())
for (nm in setdiff(names(NAMES),"mnd")) {
  nat_deaths[[nm]] <- read_wonder(raw_files$national[[nm]],CODES[[nm]])
}
national <- data.frame(month=expected_months)
for (nm in names(nat_deaths)) {
  names(nat_deaths[[nm]])[2] <- paste0("deaths_",nm)
  national <- dplyr::left_join(national,nat_deaths[[nm]],by="month")
}
national <- dplyr::left_join(national,nat_rsv,by="month")
national$t <- seq_len(nrow(national))
national$year <- as.integer(format(national$month,"%Y"))
national$moy <- as.integer(format(national$month,"%m"))
national$month_factor <- factor(national$moy,levels=1:12)
national$days <- get_days(national$month)
national$sin1 <- sin(2*pi*national$t/12)
national$cos1 <- cos(2*pi*national$t/12)
national$sin2 <- sin(4*pi*national$t/12)
national$cos2 <- cos(4*pi*national$t/12)
national$lab_era <- as.integer(national$month>=as.Date("2015-10-01"))
if (anyNA(national[,paste0("deaths_",names(NAMES))]))
  fail("Missing national deaths in one or more comparator files")
write_out(national,"04_national_matched_monthly.csv")

# ========== 4. NATIONAL comparator regressions ===============================
# Exactly matching original technical Poisson/HAC framework: 4df natural spline
# in time, two Fourier harmonics, lab-era step, offset(log calendar days),
# Newey-West lag=3 prewhite=FALSE. 2013-2019 results are pilot-informed.
fit_national <- function(dat, outcome, period, season="two_harmonics") {
  y <- paste0("deaths_",outcome)
  ss <- switch(season,
      two_harmonics="sin1+cos1+sin2+cos2",
      month_fixed="month_factor",
      one_harmonic="sin1+cos1")
  form <- stats::as.formula(paste(y,"~ rsv_mean01_5pp + splines::ns(t,df=4) +",
                                  ss,"+ lab_era + offset(log(days))"))
  fit <- stats::glm(form, data=dat, family=stats::poisson(link="log"),
                    control=stats::glm.control(maxit=100))
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))))
    fail("National Poisson failed: ",outcome," / ",period," / ",season)
  v <- sandwich::NeweyWest(fit,lag=3,prewhite=FALSE,adjust=TRUE)
  b <- unname(stats::coef(fit)[["rsv_mean01_5pp"]])
  se <- sqrt(unname(v["rsv_mean01_5pp","rsv_mean01_5pp"]))
  if (!is.finite(se) || se <= 0) fail("Bad HAC SE: ",outcome)
  resid <- stats::residuals(fit,type="pearson")
  data.frame(outcome=outcome,period=period,season=season,
             n_months=nrow(dat),n_deaths=sum(dat[[y]]),beta=b,se_HAC=se,
             IRR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),
             p_HAC=2*stats::pnorm(-abs(b/se)),
             dispersion=sum(resid^2)/stats::df.residual(fit),
             LjungBox_p12=stats::Box.test(resid,lag=12,type="Ljung-Box")$p.value)
}
national_res <- list()
for (period in c("all_2011_2019","later_2013_2019")) {
  d <- if (period=="all_2011_2019") national else national[national$year>=2013,,drop=FALSE]
  for (outcome in names(NAMES))
    for (season in c("two_harmonics","month_fixed","one_harmonic"))
      national_res[[length(national_res)+1L]] <- fit_national(d,outcome,period,season)
}
national_res <- dplyr::bind_rows(national_res)
write_out(national_res,"05_national_comparator_results.csv")
log_msg("National comparator stage complete: ",nrow(national_res)," fitted models")

# ========== 5. REGIONAL input completeness / common sample ===================
# If a comparator file was named 'national' but is actually HHS-region data,
# use that SAME ORIGINAL FILE for the regional panel rather than requiring
# a second download. For national-only tables, no fallback is applied.
is_hhs_export <- function(path) {
  if (!file.exists(path)) return(FALSE)
  first <- readLines(path,n=1L,warn=FALSE)
  length(first)==1L && grepl("HHS Region Code",first,fixed=TRUE)
}
for (nm in setdiff(names(NAMES),"mnd")) {
  if (!file.exists(raw_files$regional[[nm]]) &&
      is_hhs_export(raw_files$national[[nm]])) {
    raw_files$regional[[nm]] <- raw_files$national[[nm]]
    log_msg("Using mislabeled national file as the HHS regional export for ",
            nm,": ",raw_files$regional[[nm]])
  }
}
reg_present <- file.exists(raw_files$regional)
file_audit <- data.frame(type=c(rep("national",4),rep("regional",4)),
                         outcome=rep(names(NAMES),2),
                         path=c(unname(raw_files$national),unname(raw_files$regional)),
                         exists=c(file.exists(raw_files$national),reg_present))
write_out(file_audit,"00_file_audit.csv")
panel_ok <- FALSE
panel_status <- "NOT RUN: regional mortality files absent"
error_rows <- list()
safe_panel <- function(label,fun) {
  tryCatch(fun(),error=function(e) {
    error_rows[[length(error_rows)+1L]] <<-
      data.frame(model=label,error=conditionMessage(e))
    NULL
  })
}

if (!all(reg_present)) {
  log_msg("Regional analysis SKIPPED. Missing: ",
          paste(raw_files$regional[!reg_present],collapse="; "))
} else {
  log_msg("Four regional mortality files found; checking coverage and metadata")
  regional_deaths <- lapply(names(NAMES),function(nm) {
    x <- read_wonder(raw_files$regional[[nm]],CODES[[nm]],regional=TRUE)
    x$outcome <- nm
    x
  })
  names(regional_deaths) <- names(NAMES)
  fullgrid <- expand.grid(region=1:10,month=expected_months)
  fullgrid$month <- as.Date(fullgrid$month,origin="1970-01-01")
  panel <- fullgrid
  coverage <- list()
  for (nm in names(NAMES)) {
    d <- regional_deaths[[nm]]
    cc <- dplyr::left_join(fullgrid,d,by=c("region","month"))
    coverage[[nm]] <- cc |>
      dplyr::group_by(region) |>
      dplyr::summarise(outcome=nm,n_months=dplyr::n(),
          n_reported=sum(!is.na(deaths)),n_missing_or_suppressed=sum(is.na(deaths)),
          n_explicit_suppression=sum(grepl("suppres",suppression_label,
                                        ignore.case=TRUE),na.rm=TRUE),
          .groups="drop")
    cc <- cc[,c("region","month","deaths")]
    names(cc)[3] <- paste0("deaths_",nm)
    panel <- dplyr::left_join(panel,cc,by=c("region","month"))
  }
  coverage <- dplyr::bind_rows(coverage)
  write_out(coverage,"06_regional_death_coverage.csv")
  # Reconcile totals ONLY for months with all TEN original regional counts
  # observed. NEVER calculate the difference from national totals for a month
  # containing a suppressed or absent region; this could disclose small counts.
  reconciliation <- list()
  for (nm in names(NAMES)) {
    dd <- regional_deaths[[nm]]
    dd <- dplyr::right_join(dd, fullgrid, by=c("region","month"))
    full <- dd |>
      dplyr::group_by(month) |>
      dplyr::summarise(n_available=sum(!is.na(deaths)),
                       total_if_complete=if (all(!is.na(deaths))) sum(deaths)
                                         else NA_real_, .groups="drop")
    national_ref <- nat_deaths[[nm]]
    national_ref$month <- as.Date(national_ref$month)
    full <- dplyr::left_join(full,national_ref,by="month")
    nat_col <- paste0("deaths_",nm)
    full$identical_when_complete <- ifelse(full$n_available==10,
                      full$total_if_complete==full[[nat_col]],NA)
    # Safe output: do not expose an inferred suppressed count by subtraction.
    reconciliation[[nm]] <- data.frame(outcome=nm,month=full$month,
                         n_regions_reported=full$n_available,
                         complete_reconciles=full$identical_when_complete)
  }
  reconciliation <- dplyr::bind_rows(reconciliation)
  write_out(reconciliation,"06b_regional_national_reconciliation.csv")
  if (any(reconciliation$complete_reconciles==FALSE,na.rm=TRUE))
    fail("One or more FULLY OBSERVED HHS-region monthly totals do not agree with national mortality exports. Check CDC WONDER export settings.")
  panel <- dplyr::left_join(panel, reg_rsv, by=c("region","month"))
  panel$year <- as.integer(format(panel$month,"%Y"))
  panel$moy <- as.integer(format(panel$month,"%m"))
  panel$days <- get_days(panel$month)
  deathcols <- paste0("deaths_",names(NAMES))
  panel$primary_complete <- stats::complete.cases(panel[,c(deathcols,"rsv_mean01_5pp")])
  panel$low_PCR_volume <- is.na(panel$tests_pcr) | panel$tests_pcr<100
  write_out(panel,"07_regional_matched_monthly.csv")
  by_region <- panel |>
    dplyr::group_by(region) |>
    dplyr::summarise(n_complete=sum(primary_complete),
                     n_missing= sum(!primary_complete),
                     n_PCR_below_100=sum(low_PCR_volume),.groups="drop")
  write_out(by_region,"08_regional_panel_quality.csv")
  obs <- sum(panel$primary_complete)
  panel_ok <- obs >= MIN_PANEL_COVERAGE*nrow(panel) &&
              all(by_region$n_complete >= MIN_MONTHS_PER_REGION) &&
              all(1:10 %in% panel$region[panel$primary_complete])
  if (!panel_ok) {
    panel_status <- paste0("NOT RUN: common sample too sparse (",obs," / ",nrow(panel),
               " region-months); check coverage CSV")
    log_msg(panel_status)
  } else {
    panel_status <- paste0("RAN: ",obs," matched region-months across ten regions")
    log_msg(panel_status)
  }
}

# ========== 6. PANEL Poisson with high-dimensional FEs ========================
# FE: region*year + region*month-of-year + year-month.
# Each region contributes one entire serially correlated unit: only TEN clusters.
# Approximate cluster-t(9) intervals are DESCRIPTIVE. Pairs-bootstrap percentile
# intervals are an extra stability diagnostic, not valid exact few-cluster tests.
prepare_panel <- function(d) {
  d$region_id <- as.character(d$region)
  d$region_year <- paste(d$region_id,d$year,sep="__")
  d$region_moy <- paste(d$region_id,d$moy,sep="__")
  d$yearmonth_fe <- format(d$month,"%Y-%m")
  d$log_days <- log(d$days)
  d
}
panel_model <- function(d, outcome, xvar="rsv_mean01_5pp") {
  if (nrow(d)<150 || length(unique(d$region_id))<5)
    fail("Not enough region-month data for panel model")
  if (anyNA(d[,c(paste0("deaths_",outcome),xvar,"log_days")]))
    fail("Missing required values supplied to fepois")
  dd <- d
  dd$y <- dd[[paste0("deaths_",outcome)]]
  dd$x <- dd[[xvar]]
  m <- fixest::fepois(y ~ x | region_year + region_moy + yearmonth_fe,
                      data=dd, offset=~log_days, nthreads=1)
  if (stats::nobs(m)!=nrow(dd) ||
      !"x" %in% names(stats::coef(m)) ||
      !is.finite(unname(stats::coef(m)[["x"]])))
    fail("Panel fit dropped observations or RSV exposure is not identifiable")
  m
}
panel_effect <- function(m,d,outcome,period,scenario) {
  b <- unname(stats::coef(m)[["x"]])
  # fixest small-sample cluster adjustment: 10 clusters, 9df t-based interval
  # is indicative only; few-cluster inference remains a major limitation.
  v <- stats::vcov(m, vcov=~region_id)
  se <- sqrt(unname(v["x","x"]))
  g <- length(unique(d$region_id))
  critical <- stats::qt(0.975,df=g-1)
  pearson <- stats::residuals(m,type="pearson")
  residual_df <- if (is.numeric(m$nparams) && length(m$nparams)==1L)
    nrow(d)-m$nparams else NA_real_
  disp <- if (is.finite(residual_df) && residual_df>0)
    sum(pearson^2)/residual_df else NA_real_
  if (!is.finite(se) || se<=0) fail("Invalid panel cluster SE")
  data.frame(outcome=outcome,period=period,scenario=scenario,
             n_obs=nrow(d),n_regions=g,coef=b,se_cluster_approx=se,
             IRR=exp(b),CI_low_cluster_t=exp(b-critical*se),
             CI_high_cluster_t=exp(b+critical*se),
             p_cluster_t_approx=2*stats::pt(-abs(b/se),df=g-1),
             dispersion=disp,
             notes="10-cluster approx; refer to bootstrap and sensitivity")
}

if (panel_ok) {
  common <- prepare_panel(panel[panel$primary_complete,,drop=FALSE])
  # Identification check: residual RSV variability after absorbing the SAME
  # fixed effects as the mortality models. Small residual variation means
  # an apparent null/positive association is not an informative test.
  variation <- tryCatch({
    mm <- fixest::feols(rsv_mean01_5pp ~ 1 |
              region_year + region_moy + yearmonth_fe,
              data=common,nthreads=1)
    raw_sd <- stats::sd(common$rsv_mean01_5pp)
    residual_sd <- stats::sd(stats::residuals(mm))
    data.frame(raw_sd=raw_sd,residual_sd=residual_sd,
               residual_fraction=residual_sd/raw_sd)
  },error=function(e) {
    log_msg("Exposure identification diagnostic failed: ",conditionMessage(e))
    NULL
  })
  if (!is.null(variation)) {
    write_out(variation,"08b_regional_exposure_identification.csv")
    if (!is.finite(variation$residual_fraction) ||
        variation$residual_fraction < .02) {
      panel_ok <- FALSE
      panel_status <- "NOT RUN: insufficient regional PCR positivity variation after fixed effects"
      log_msg(panel_status)
    }
  } else {
    panel_ok <- FALSE
    panel_status <- "NOT RUN: cannot validate regional exposure identification"
  }
}
if (panel_ok) {
  panel_results <- list()
  add_model <- function(d,outcome,period,scenario="main",xvar="rsv_mean01_5pp") {
    m <- panel_model(d,outcome,xvar)
    panel_effect(m,d,outcome,period,scenario)
  }
  for (outcome in names(NAMES)) {
    z <- safe_panel(paste(outcome,"main"),function()
      add_model(common,outcome,"2011_2019"))
    if (!is.null(z)) panel_results[[length(panel_results)+1L]] <- z
    later <- common[common$year>=2013,,drop=FALSE]
    z <- safe_panel(paste(outcome,"later"),function()
      add_model(later,outcome,"2013_2019","restricted_2013"))
    if (!is.null(z)) panel_results[[length(panel_results)+1L]] <- z
  }
  write_out(dplyr::bind_rows(panel_results),"09_regional_primary_comparators.csv")
  if (!"mnd" %in% vapply(panel_results,function(z) z$outcome,character(1)))
    panel_status <- "REGIONAL MAIN MND MODEL FAILED: review model_errors.csv"

  # Prespecified MND surveillance specification diagnostics, no best-model
  # selection. Each variant uses the same outcome-specific eligible rows.
  mnd_sens <- list()
  for (sc in c("exclude_CDC_outliers","PCR_tests_at_least_100",
               "adjust_log_PCR_volume", "lag0","lag1","lead2")) {
    z <- safe_panel(paste0("MND_",sc),function() {
      d <- common
      vv <- switch(sc, exclude_CDC_outliers="rsv_clean_mean01_5pp",
             PCR_tests_at_least_100="rsv_mean01_5pp",
             adjust_log_PCR_volume="rsv_mean01_5pp",
             lag0="rsv_lag0_5pp",lag1="rsv_lag1_5pp",
             lead2="rsv_lead2_5pp")
      if (sc=="PCR_tests_at_least_100") d <- d[!d$low_PCR_volume,,drop=FALSE]
      d <- d[is.finite(d[[vv]]),,drop=FALSE]
      d <- prepare_panel(d)
      if (sc=="adjust_log_PCR_volume") {
        d$y <- d$deaths_mnd
        d$x <- d[[vv]]
        d$log_tests <- log(pmax(d$tests_pcr,1))
        m <- fixest::fepois(y ~ x + log_tests | region_year + region_moy + yearmonth_fe,
                  data=d,offset=~log_days,nthreads=1)
        if (stats::nobs(m)!=nrow(d)) fail("Adjusted PCR volume fit dropped rows")
      } else m <- panel_model(d,"mnd",vv)
      panel_effect(m,d,"mnd",paste0(min(d$year),"_",max(d$year)),sc)
    })
    if (!is.null(z)) mnd_sens[[length(mnd_sens)+1L]] <- z
  }
  write_out(dplyr::bind_rows(mnd_sens),"10_regional_MND_sensitivities.csv")

  # Delete one entire HHS region; compare direction and magnitude, NOT p<.05.
  loo <- list()
  for (r in 1:10) {
    d <- prepare_panel(common[common$region!=r,,drop=FALSE])
    z <- safe_panel(paste0("MND_exclude_region_",r),function() {
      m <- panel_model(d,"mnd")
      z <- panel_effect(m,d,"mnd","2011_2019","LOO")
      z$omitted_region <- r
      z
    })
    if (!is.null(z)) loo[[length(loo)+1L]] <- z
  }
  write_out(dplyr::bind_rows(loo),"11_leave_one_region_out.csv")

  # Cluster pairs bootstrap (whole regions resampled, all months and all four
  # outcomes jointly). Fresh region ID per *drawn copy* prevents FE collisions.
  # This is NOT a wild-cluster bootstrap and is NOT a definitive p-value test.
  # With 10 sampling units, percentile bands are approximate diagnostics only.
  if (B_CLUSTER>0L && length(panel_results)>0) {
    set.seed(BOOT_SEED)
    boot_draws <- list()
    for (i in seq_len(B_CLUSTER)) {
      chosen <- sample(1:10,size=10,replace=TRUE)
      if (length(unique(chosen))<5) next
      replicate_data <- dplyr::bind_rows(lapply(seq_along(chosen),function(j) {
        d <- common[common$region==chosen[[j]],,drop=FALSE]
        d$region <- j  # independent replica of sampled regional cluster
        d
      }))
      replicate_data <- prepare_panel(replicate_data)
      b <- setNames(rep(NA_real_,length(NAMES)),names(NAMES))
      for (nm in names(NAMES)) {
        f <- tryCatch(panel_model(replicate_data,nm),error=function(e) NULL)
        if (!is.null(f)) b[[nm]] <- unname(stats::coef(f)[["x"]])
      }
      if (all(is.finite(b))) {
        boot_draws[[length(boot_draws)+1L]] <- data.frame(draw=i,
                 beta_mnd=b[["mnd"]],beta_parkinson=b[["parkinson"]],
                 beta_transport=b[["transport"]],beta_respiratory=b[["respiratory"]])
      }
      if (i%%50L==0L) log_msg("Regional resampling diagnostics: ",i," / ",B_CLUSTER)
    }
    boot <- dplyr::bind_rows(boot_draws)
    write_out(boot,"12_cluster_pairs_bootstrap_draws.csv")
    if (nrow(boot)>=max(80L,ceiling(0.70*B_CLUSTER))) {
      bootsum <- dplyr::bind_rows(lapply(names(NAMES),function(nm) {
        x <- boot[[paste0("beta_",nm)]]
        q <- as.numeric(stats::quantile(x,c(.025,.5,.975),na.rm=TRUE))
        data.frame(outcome=nm,B_success=nrow(boot),
                 IRR_2.5_pct=exp(q[1]),IRR_median=exp(q[2]),IRR_97.5_pct=exp(q[3]))
      }))
      write_out(bootsum,"13_pairs_bootstrap_percentiles.csv")
      ctr <- dplyr::bind_rows(lapply(setdiff(names(NAMES),"mnd"),function(nm) {
        v <- boot$beta_mnd-boot[[paste0("beta_",nm)]]
        q <- as.numeric(stats::quantile(v,c(.025,.5,.975)))
        data.frame(contrast=paste0("MND_vs_",nm),B_success=nrow(boot),
                   log_IRR_difference_median=q[2],
                   ratio_of_IRRs_low=exp(q[1]),ratio_of_IRRs_median=exp(q[2]),
                   ratio_of_IRRs_high=exp(q[3]))
      }))
      write_out(ctr,"14_paired_outcome_contrasts_bootstrap.csv")
    } else log_msg("Too few successful bootstrap draws (",nrow(boot),
                   "); no bootstrap intervals reported")
  }
}

if (length(error_rows)) write_out(dplyr::bind_rows(error_rows),"model_errors.csv")
summary_lines <- c(
  "RSV-MND comparator and regional robustness: analysis summary",
  "===========================================================",
  "Study: US 2011-2019 UCD G12.2 with RSV PCR NREVSS positivity, per 5 pp",
  "Pilot-derived specifications, not prospectively confirmatory.",
  paste("National comparator models:",nrow(national_res)),
  paste("Regional analysis:",panel_status),
  paste("Panel bootstrap draws requested:",B_CLUSTER),
  "Any regional interval uses only 10 HHS clusters; treat as approximate.",
  "Parkinson's is disease comparator; transport accidents are a falsification",
  "outcome, not a validated matched negative control; J09-J18 a benchmark.",
  "National models use Poisson with Newey-West/HAC(3) standard errors.",
  "Regional models use Poisson HDFE and region-clustered SE, plus pairs bootstrap.",
  "Suppressed/absent death cells remain NA; no numerical reconstruction.",
  "RSV surveillance region = laboratory region, not decedent residence.",
  "No individual ALS cases/RSV-attributable deaths identified; no causality.",
  "DO NOT select favourable lags, periods, outcomes or isolated p-values.",
  if (length(error_rows)) paste("Regional model errors:",length(error_rows))
  else "Regional model errors: none recorded",
  "All pre-existing pilot data and scripts are unchanged."
)
writeLines(summary_lines,file.path(OUT,"analysis_summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_session_info.txt"))

# ========== 7. Publication-ready exploratory plots ===========================
plotdf <- national[,c("month","deaths_mnd","pct_pcr","tests_pcr")]
p <- ggplot2::ggplot(plotdf,ggplot2::aes(x=month)) +
  ggplot2::geom_line(ggplot2::aes(y=pct_pcr),linewidth=.55) +
  ggplot2::labs(x="Month",y="RSV PCR positivity (%)",
                title="RSV PCR surveillance, United States, 2011-2019") +
  ggplot2::theme_bw(base_size=11)
ggplot2::ggsave(file.path(OUT,"figure_01_RSV_national.png"),p,width=9,height=4.3,dpi=250)
p2 <- ggplot2::ggplot(national_res,ggplot2::aes(x=IRR,y=outcome))+
  ggplot2::geom_vline(xintercept=1,linetype=2,color="grey60")+
  ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),
                         width=.22,orientation="y")+
  ggplot2::geom_point(size=2)+
  ggplot2::facet_grid(period~season)+
  ggplot2::labs(x="IRR per 5 percentage points RSV PCR positivity (95% HAC CI)",
                y="Outcome",title="National outcome comparator models")+
  ggplot2::theme_bw(base_size=10)
ggplot2::ggsave(file.path(OUT,"figure_02_national_forest.png"),p2,width=11,height=6,dpi=250)
if (panel_ok && exists("panel_results") && length(panel_results)) {
  pr <- dplyr::bind_rows(panel_results)
  p3 <- ggplot2::ggplot(pr,ggplot2::aes(x=IRR,y=outcome))+
    ggplot2::geom_vline(xintercept=1,linetype=2,color="grey60")+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low_cluster_t,
                                       xmax=CI_high_cluster_t),
                         width=.22,orientation="y")+
    ggplot2::geom_point(size=2)+
    ggplot2::facet_wrap(~period)+ggplot2::theme_bw(base_size=11)+
    ggplot2::labs(x="IRR per 5 pp (approximate 10-cluster t interval)",
                  y="Outcome",title="Region-matched robustness analysis")
  ggplot2::ggsave(file.path(OUT,"figure_03_regional_forest.png"),p3,
                  width=9,height=4.8,dpi=250)
}
log_msg("COMPLETE. Read analysis_results/analysis_summary.txt and ",
        "05_national_comparator_results.csv first.")
