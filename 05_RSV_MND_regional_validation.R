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
# 05_RSV_MND_regional_validation.R
# RSV-MND national-to-regional robustness extension (2011-2019)
# Date: 24 October 2026. Exploratory study following national pilot.
# Run from the original project root (containing raw_data/ and analysis_results/):
#   source("05_RSV_MND_regional_validation.R")
# Existing national results remain unchanged. All outputs go into
# regional_validation_results/. Files containing suppressed/absent counts
# are NEVER completed by differencing regional and national totals.

# ---------------------------- configuration -------------------------------
START <- as.Date("2011-01-01")
END <- as.Date("2019-12-01")
B <- 299L                      # region-pairs bootstrap, diagnostic only
SEED <- 20261024L
MIN_PCR_TESTS <- 100            # fractional monthly tests from week allocation
OUT <- "regional_validation_results"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

pkgs <- c("fixest", "dplyr", "sandwich", "ggplot2")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop("Install packages: install.packages(c(",
  paste(shQuote(missing_pkgs), collapse=", "), "))", call.=FALSE)

logfile <- file.path(OUT,"analysis_log.txt")
if(file.exists(logfile)) file.remove(logfile)
msg <- function(...) {s <- paste0(...,collapse=""); cat(s,"\n");
                     cat(s,"\n",file=logfile,append=TRUE)}
fail <- function(...) stop(paste0(...,collapse=""),call.=FALSE)
save_csv <- function(d, name) utils::write.csv(d,file.path(OUT,name),row.names=FALSE,na="")
first_existing <- function(paths, label) {
  ok <- paths[file.exists(paths)]
  if (!length(ok)) fail("Missing ",label,". Looked for: ",paste(paths,collapse="; "))
  ok[[1]]
}
allmonths <- seq(START,END,by="month")
days_in_month <- function(x) {
  y <- as.integer(format(x,"%Y")); m <- as.integer(format(x,"%m"))
  nextm <- as.Date(sprintf("%04d-%02d-01",y+as.integer(m==12), m%%12+1L))
  as.integer(nextm-x)
}
outcomes <- c("mnd","parkinson","transport","respiratory")
codes <- c(mnd="G12.2",parkinson="G20",transport="V01-V99",respiratory="J09-J18")
file_options <- list(
 mnd=c("raw_data/MND_national_monthly.txt",
       "raw_data/MND_G12_2_HHS_monthly_2011_2019.txt"),
 parkinson=c("raw_data/Parkinson_G20_national_monthly.txt",
             "raw_data/Parkinson_G20_HHS_monthly_2011_2019.txt"),
 transport=c("raw_data/Transport_national_monthly.txt",
             "raw_data/Transport_V01_V99_HHS_monthly_2011_2019.txt"),
 respiratory=c("raw_data/Respiratory_national_monthly.txt",
             "raw_data/Respiratory_J09_J18_HHS_monthly_2011_2019.txt")
)
paths <- vapply(outcomes, function(z) first_existing(file_options[[z]],paste(z,"regional WONDER file")), "")
msg("Input files: ", paste(names(paths),paths,sep="=",collapse=" | "))

# --------------------------- A. Read CDC WONDER -----------------------------
# CDC WONDER sometimes hides suppressed regional rows rather than labelling
# them. They must remain unavailable; do NOT impute them as zero or infer
# death counts using independently observed national totals.
read_wonder_regional <- function(path, code) {
  lines <- readLines(path,encoding="UTF-8",warn=FALSE)
  first <- which(grepl("Month Code",lines,fixed=TRUE) &
                 grepl("Deaths",lines,fixed=TRUE) & grepl("\t",lines,fixed=TRUE))
  if (!length(first)) fail("No WONDER TSV header: ",path)
  foot <- lines[grepl("^\\s*\\\"?(Dataset:|ICD-10 Codes:|Group By:|Year/Month:)",lines)]
  code_line <- foot[grepl("ICD-10 Codes:",foot,fixed=TRUE)]
  data_line <- foot[grepl("Dataset:",foot,fixed=TRUE)]
  group_line <- foot[grepl("Group By:",foot,fixed=TRUE)]
  if (!length(code_line) || !length(data_line) || !length(group_line))
    fail("Missing CDC WONDER query metadata in ",path)
  p <- switch(code, "G12.2"="G12\\.2(?![0-9])", "G20"="G20(?![0-9])",
               "V01-V99"="V01\\s*[-:]\\s*V99", "J09-J18"="J09\\s*[-:]\\s*J18")
  if (!grepl(p,toupper(gsub("[–—]","-",code_line[[1]])),perl=TRUE))
    fail("Incorrect ICD code in: ",path," / ",code_line[[1]])
  if (!grepl("Underlying Cause of Death",data_line[[1]],ignore.case=TRUE) ||
      !grepl("HHS Region",group_line[[1]],ignore.case=TRUE) ||
      !grepl("Month",group_line[[1]],ignore.case=TRUE))
    fail("Dataset/grouping incompatible with HHS region-month model: ",path)
  x <- lines[first[[1]]:length(lines)]
  ix <- which(grepl('^\\s*"?---"?\\s*$',x) | grepl('^"?Dataset:',x))
  if (length(ix)) x <- head(x,ix[[1]]-1L)
  d <- utils::read.delim(text=paste(x,collapse="\n"),sep="\t",quote='"',
                         comment.char="",check.names=FALSE,stringsAsFactors=FALSE,
                         colClasses="character",fill=TRUE)
  need <- c("HHS Region Code","Month Code","Deaths")
  if (!all(need %in% names(d))) fail("Missing WONDER fields in: ",path)
  d <- d[grepl("^20[0-9]{2}/(0[1-9]|1[0-2])$",d[["Month Code"]]) &
         !is.na(d[["Month Code"]]),,drop=FALSE]
  d$month <- as.Date(paste0(gsub("/","-",d[["Month Code"]]),"-01"))
  d <- d[d$month>=START & d$month<=END,,drop=FALSE]
  d$region <- suppressWarnings(as.integer(sub("^HHS0?","",trimws(d[["HHS Region Code"]]))))
  if (!nrow(d) || anyNA(d$region) || any(!d$region %in% 1:10))
    fail("Invalid region IDs in ",path)
  if (anyDuplicated(paste(d$region,d$month))) fail("Duplicate region/month in ",path)
  raw_counts <- trimws(d$Deaths)
  d$deaths <- suppressWarnings(as.numeric(gsub(",","",raw_counts,fixed=TRUE)))
  if (any(!is.na(d$deaths) & (d$deaths<0 | d$deaths!=round(d$deaths))))
    fail("Invalid mortality counts in ",path)
  if (any(!is.na(d$deaths) & d$deaths>0 & d$deaths<10))
    fail("Unexpected unsuppressed count below 10 in ",path)
  grid <- expand.grid(region=1:10,month=allmonths)
  grid$month <- as.Date(grid$month,origin="1970-01-01")
  z <- merge(grid,d[,c("region","month","deaths")],by=c("region","month"),all.x=TRUE,
             sort=TRUE)
  if (nrow(z)!=1080L) fail("Regional import did not produce 1080 cells: ",path)
  attr(z,"observed_rows") <- nrow(d)
  attr(z,"suppressed_or_absent") <- sum(is.na(z$deaths))
  z
}

mort <- setNames(lapply(outcomes,function(nm) read_wonder_regional(paths[[nm]],codes[[nm]])),outcomes)
quality <- do.call(rbind,lapply(outcomes,function(nm) {
  x <- mort[[nm]]
  data.frame(outcome=nm,source_path=paths[[nm]],expected_rows=1080L,
             observed_rows=attr(x,"observed_rows"),
             missing_or_suppressed=attr(x,"suppressed_or_absent"),
             reported_death_counts=sum(!is.na(x$deaths)),
             # Report no partial death totals where suppression exists.
             all_cells_complete=all(!is.na(x$deaths)))
}))
save_csv(quality,"00_regional_mortality_audit.csv")
msg("Death-data audit complete. Missing/suppressed: ",
    paste(quality$outcome,quality$missing_or_suppressed,sep="=",collapse=" | "))

# ------------------------- B. RSV and national audit ------------------------
rsvfile <- first_existing(c("analysis_results/01_regional_rsv_monthly.csv",
                            "01_regional_rsv_monthly.csv"),
                            "prior regional PCR surveillance output from script 04 v2")
rsv <- utils::read.csv(rsvfile,stringsAsFactors=FALSE)
required_rsv <- c("region","month","rsv_mean01_5pp","rsv_clean_mean01_5pp",
                  "rsv_lag0_5pp","rsv_lag1_5pp","rsv_lead2_5pp","tests_pcr")
if (!all(required_rsv %in% names(rsv)))
  fail("The RSV output is missing fields: ",
       paste(setdiff(required_rsv,names(rsv)),collapse=", "),". Run script 04 v2 first.")
rsv$month <- as.Date(rsv$month)
if (anyNA(rsv$month) || anyDuplicated(paste(rsv$region,rsv$month)))
  fail("Regional RSV contains invalid dates or duplicate region-months")
rsv <- rsv[rsv$month>=START & rsv$month<=END, ,drop=FALSE]
if (!all(is.finite(rsv$rsv_mean01_5pp)))
  fail("RSV positivity has missing/invalid values; investigate exposure before fitting")

panel <- mort[["mnd"]][,c("region","month")]
for (nm in outcomes) {
  z <- mort[[nm]][,c("region","month","deaths")]
  names(z)[3] <- paste0("deaths_",nm)
  panel <- dplyr::left_join(panel,z,by=c("region","month"))
}
panel <- dplyr::left_join(panel,rsv,by=c("region","month"))
if (nrow(panel)!=1080L || anyNA(panel$region) || anyNA(panel$month))
  fail("Region/month join failed")
panel$year <- as.integer(format(panel$month,"%Y"))
panel$moy <- as.integer(format(panel$month,"%m"))
panel$days <- days_in_month(panel$month)
panel$low_pcr_volume <- is.na(panel$tests_pcr) | panel$tests_pcr < MIN_PCR_TESTS
names_deaths <- paste0("deaths_",outcomes)
panel$complete_for_comparison <- stats::complete.cases(panel[,c(names_deaths,"rsv_mean01_5pp")])
# A common identical sample is used for all primary mortality comparisons.
common <- panel[panel$complete_for_comparison,,drop=FALSE]
coverage <- do.call(rbind,lapply(outcomes,function(nm) {
  data.frame(outcome=nm,region=1:10,
    available=sapply(1:10,function(i) sum(!is.na(panel[[paste0("deaths_",nm)]][panel$region==i]))),
    common_sample=sapply(1:10,function(i) sum(common$region==i)))
}))
save_csv(coverage,"01_regional_coverage_by_outcome.csv")
if (nrow(common)<0.80*1080 || any(table(factor(common$region,levels=1:10))<60))
  fail("Too little common regional panel coverage; see coverage CSV")
msg("Complete common region-month sample: ",nrow(common)," / 1080")
# Reconcile only fully observed months with prior independently sourced national
# MND counts. Do NOT calculate differences for any incomplete month.
nat_path <- first_existing(c("analysis_results/04_national_matched_monthly.csv",
                             "04_national_matched_monthly.csv"),
                           "prior national matched results")
nat <- utils::read.csv(nat_path,stringsAsFactors=FALSE)
nat$month <- as.Date(nat$month)
if (anyDuplicated(nat$month) || !all(c("month",names_deaths) %in% names(nat)))
  fail("Invalid national matched CSV from previous analysis")
recon <- do.call(rbind,lapply(outcomes,function(nm) {
  col <- paste0("deaths_",nm)
  d <- panel[,c("month",col)]
  ss <- dplyr::group_by(d,month) |>
    dplyr::summarise(n_observed=sum(!is.na(.data[[col]])),
                     sum_if_complete=if (all(!is.na(.data[[col]]))) sum(.data[[col]])
                                     else NA_real_,.groups="drop")
  nn <- nat[,c("month",col)]
  names(nn)[2] <- "national_deaths"
  r <- dplyr::left_join(ss,nn,by="month")
  if (anyNA(r$national_deaths) || nrow(r)!=108L) fail("National comparator reconciliation reference missing months")
  # Do not write sums or national totals to reconciliation output; preserve
  # confidentiality for suppressed cells (including totals inferable by subtraction).
  data.frame(outcome=nm,months_complete=sum(r$n_observed==10L),
             months_incomplete=sum(r$n_observed!=10L),
             complete_months_consistent=all(r$sum_if_complete[r$n_observed==10L] ==
                                            r$national_deaths[r$n_observed==10L]))
}))
save_csv(recon,"02_regional_national_reconciliation.csv")
if (any(!recon$complete_months_consistent))
  fail("Fully observed HHS months disagree with national historical files. Check causes/settings")
msg("All fully observed region-month sums reconcile with prior national series.")

# Identifiable within-region deviation after removing national temporal trends.
prepare <- function(d) {
  d$region <- as.integer(d$region)
  d$region_year <- interaction(d$region,d$year,drop=TRUE)
  d$region_moy <- interaction(d$region,d$moy,drop=TRUE)
  d$yearmonth_fe <- format(d$month,"%Y-%m")
  d$log_days <- log(d$days)
  d
}
common <- prepare(common)
write_common_panel <- common[,c("region","month","year","moy",names_deaths,
                               "rsv_mean01_5pp","rsv_clean_mean01_5pp",
                               "tests_pcr","low_pcr_volume")]
# Same observed records already supplied by CDC; no imputed/suppressed values.
save_csv(write_common_panel,"03_common_panel_used.csv")

variation <- function(d,label) {
  m <- fixest::feols(rsv_mean01_5pp ~ 1 | region_year + region_moy + yearmonth_fe,
                    data=d,nthreads=1)
  raw <- stats::sd(d$rsv_mean01_5pp)
  resid <- stats::sd(stats::residuals(m))
  data.frame(period=label,n_obs=nrow(d),raw_sd=raw,residual_sd=resid,
             residual_fraction=resid/raw)
}
variation_table <- dplyr::bind_rows(variation(common,"2011_2019"),
                                   variation(prepare(common[common$year>=2013,,drop=FALSE]),"2013_2019"))
save_csv(variation_table,"04_exposure_identification.csv")
if (any(!is.finite(variation_table$residual_fraction)) ||
    any(variation_table$residual_fraction < 0.02))
  fail("Too little RSV variation after absorbing fixed effects; see 04_exposure_identification.csv")

# -------------------------- C. Regional models ------------------------------
# Panel FE specification: region x year + region x calendar-month + national
# year-month. The last term absorbs all shared shocks and days per month.
# Within-region serial dependence: 10 region-clustered standard errors.
# t(9) CI and associated P values are APPROXIMATE, not confirmatory.
fit_panel <- function(d, outcome, xvar="rsv_mean01_5pp", add_volume=FALSE) {
  dd <- prepare(d)
  if (!all(is.finite(dd[[xvar]])) || anyNA(dd[[paste0("deaths_",outcome)]]))
    fail("Missing inputs in panel model for ",outcome," / ",xvar)
  dd$y <- dd[[paste0("deaths_",outcome)]]
  dd$x <- dd[[xvar]]
  if (add_volume) {
    if (any(!is.finite(dd$tests_pcr)) || any(dd$tests_pcr<=0))
      fail("Invalid PCR tests for volume-adjusted panel")
    dd$log_pcr_tests <- log(dd$tests_pcr)
    form <- y ~ x + log_pcr_tests | region_year + region_moy + yearmonth_fe
  } else form <- y ~ x | region_year + region_moy + yearmonth_fe
  # log(days/month) offset is absorbed exactly by year-month FE; include to
  # keep output comparable to the original national count model.
  fit <- fixest::fepois(form,data=dd,offset=~log_days,nthreads=1,
                       notes=FALSE,warn=TRUE)
  if (stats::nobs(fit)!=nrow(dd) || !"x" %in% names(stats::coef(fit)) ||
      !is.finite(stats::coef(fit)[["x"]]))
    fail("Model failed/dropped rows: ",outcome," / ",xvar)
  fit
}
get_panel_result <- function(f,d,outcome,period,spec) {
  b <- unname(stats::coef(f)[["x"]])
  vv <- stats::vcov(f,vcov=~region)
  se <- sqrt(unname(vv["x","x"]))
  if (!is.finite(se) || se<=0) fail("Invalid cluster SE for ",outcome)
  G <- length(unique(d$region)); tcrit <- stats::qt(.975,G-1)
  data.frame(outcome=outcome,period=period,specification=spec,n_obs=nrow(d),
    n_regions=G,beta=b,SE_cluster=se,IRR=exp(b),
    CI_low_t_approx=exp(b-tcrit*se),CI_high_t_approx=exp(b+tcrit*se),
    p_t_approx=2*stats::pt(-abs(b/se),df=G-1),
    df_clusters=G-1,cluster_warning="ONLY 10 regions; intervals exploratory")
}
models <- list(); errors <- list()
optional <- function(name,fn) {
  tryCatch(fn(),error=function(e) {
    errors[[length(errors)+1L]] <<- data.frame(stage=name,error=conditionMessage(e))
    msg("WARNING: ",name,": ",conditionMessage(e))
    NULL
  })
}
fit_all <- function(d,period,spec="main",xvar="rsv_mean01_5pp",add_volume=FALSE) {
  d <- prepare(d)
  do.call(rbind,lapply(outcomes,function(nm) {
    f <- fit_panel(d,nm,xvar,add_volume)
    get_panel_result(f,d,nm,period,spec)
  }))
}
primary <- fit_all(common,"2011_2019")
later <- fit_all(common[common$year>=2013,,drop=FALSE],"2013_2019")
models[[1]] <- primary
models[[2]] <- later
save_csv(dplyr::bind_rows(models),"05_regional_primary_comparators.csv")
msg("Primary regional models fit for all four outcomes and both time periods.")

# Placebo lead is diagnostic only, NOT evidence of reverse causality.
mnd_sens <- list()
sens <- list(
  exclude_CDC_outliers=list(x="rsv_clean_mean01_5pp",df=common,vol=FALSE),
  PCR_tests_at_least_100=list(x="rsv_mean01_5pp",
     df=common[!common$low_pcr_volume,,drop=FALSE],vol=FALSE),
  adjust_log_PCR_volume=list(x="rsv_mean01_5pp",df=common,vol=TRUE),
  lag0=list(x="rsv_lag0_5pp",df=common,vol=FALSE),
  lag1=list(x="rsv_lag1_5pp",df=common,vol=FALSE),
  lead2_placebo=list(x="rsv_lead2_5pp",df=common,vol=FALSE))
for (sc in names(sens)) {
  ob <- sens[[sc]]
  z <- optional(paste0("MND_",sc),function() {
    dd <- prepare(ob$df[is.finite(ob$df[[ob$x]]),,drop=FALSE])
    ff <- fit_panel(dd,"mnd",ob$x,ob$vol)
    get_panel_result(ff,dd,"mnd","2011_2019",sc)
  })
  if (!is.null(z)) mnd_sens[[length(mnd_sens)+1L]] <- z
}
if(length(mnd_sens)) save_csv(dplyr::bind_rows(mnd_sens),"06_MND_specification_sensitivities.csv")

loo <- list()
for (r in 1:10) {
  z <- optional(paste0("omit_region_",r),function() {
    d <- common[common$region!=r,,drop=FALSE]
    ff <- fit_panel(d,"mnd")
    x <- get_panel_result(ff,d,"mnd","2011_2019","leave_one_region_out")
    x$omitted_region <- r
    x
  })
  if (!is.null(z)) loo[[length(loo)+1L]] <- z
}
if(length(loo)) save_csv(dplyr::bind_rows(loo),"07_MND_leave_one_region_out.csv")

# ------------------- D. Paired regional resampling --------------------------
# Whole-region pairs bootstrap preserves within-region time dependence and
# dependence BETWEEN outcome estimates. Fresh region IDs prevent merging FEs
# for duplicate draws. B=299 is diagnostic; not small-cluster-exact inference.
set.seed(SEED)
boot <- list()
for (i in seq_len(B)) {
  sampled <- sample(1:10,10,replace=TRUE)
  if(length(unique(sampled))<5L) next
  d <- dplyr::bind_rows(lapply(seq_along(sampled),function(j) {
    z <- common[common$region==sampled[[j]],,drop=FALSE]
    z$region <- j
    z
  }))
  d <- prepare(d)
  v <- setNames(rep(NA_real_,4L),outcomes)
  for (nm in outcomes) {
    fit <- tryCatch(suppressWarnings(fit_panel(d,nm)),error=function(e) NULL)
    if (!is.null(fit)) v[[nm]] <- unname(stats::coef(fit)[["x"]])
  }
  if (all(is.finite(v)))
    boot[[length(boot)+1L]] <- data.frame(draw=i,beta_mnd=v[["mnd"]],
      beta_parkinson=v[["parkinson"]],beta_transport=v[["transport"]],
      beta_respiratory=v[["respiratory"]])
  if(i%%50L==0L) msg("Bootstrap draws attempted: ",i," / ",B)
}
boots <- if(length(boot)) dplyr::bind_rows(boot) else
  data.frame(draw=integer(),beta_mnd=numeric(),beta_parkinson=numeric(),
             beta_transport=numeric(),beta_respiratory=numeric())
save_csv(boots,"08_region_pairs_bootstrap_draws.csv")
if(nrow(boots)>=max(80L,ceiling(B*0.7))) {
  pct <- do.call(rbind,lapply(outcomes,function(nm) {
    q <- unname(stats::quantile(boots[[paste0("beta_",nm)]],c(.025,.5,.975)))
    data.frame(outcome=nm,B_success=nrow(boots),
               IRR_low=exp(q[1]),IRR_median=exp(q[2]),IRR_high=exp(q[3]))
  }))
  save_csv(pct,"09_region_pairs_bootstrap_IRR_percentiles.csv")
  contrasts <- do.call(rbind,lapply(setdiff(outcomes,"mnd"),function(nm) {
    d <- boots$beta_mnd - boots[[paste0("beta_",nm)]]
    q <- unname(stats::quantile(d,c(.025,.5,.975)))
    data.frame(contrast=paste0("MND_minus_",nm),B_success=nrow(boots),
               ratio_IRRs_low=exp(q[1]),ratio_IRRs_median=exp(q[2]),
               ratio_IRRs_high=exp(q[3]))
  }))
  save_csv(contrasts,"10_paired_outcome_contrasts.csv")
} else msg("WARNING: fewer than 70% of bootstrap fits succeeded; no percentile intervals.")

# -------------------- E. NATIONAL comparator diagnostics -------------------
# Strong overdispersion in J09-J18 requires caution. HAC(3/6/12) uncertainty
# and df=4/6 trends are exploratory stress tests of the existing national fits.
# They are NOT multiple new independent confirmatory hypothesis tests.
nat$year <- as.integer(format(nat$month,"%Y"))
nat$t <- match(nat$month,allmonths)
nat$month_factor <- factor(as.integer(format(nat$month,"%m")),levels=1:12)
nat$days <- days_in_month(nat$month)
nat$lab_era <- as.integer(nat$month>=as.Date("2015-10-01"))
if(!"rsv_mean01_5pp" %in% names(nat)) fail("Missing national RSV predictor from prior CSV")
national_stress <- list()
for (nm in outcomes) for (period in c("2011_2019","2013_2019")) {
  d <- if (period=="2011_2019") nat else nat[nat$year>=2013,,drop=FALSE]
  for (df in c(4L,6L)) {
    ff <- stats::as.formula(paste0("deaths_",nm,
      " ~ rsv_mean01_5pp + splines::ns(t, df=",df,") + month_factor + lab_era + offset(log(days))"))
    z <- optional(paste0("national_",nm,"_",period,"_trenddf",df),function() {
      fit <- stats::glm(ff,data=d,family=stats::poisson(link="log"),
                        control=stats::glm.control(maxit=100))
      if(!isTRUE(fit$converged)) fail("National GLM failed to converge")
      b <- unname(stats::coef(fit)[["rsv_mean01_5pp"]])
      pear <- stats::residuals(fit,type="pearson")
      if(!is.finite(b)) fail("National RSV coefficient not identifiable")
      do.call(rbind,lapply(c(3L,6L,12L),function(lag) {
        cv <- sandwich::NeweyWest(fit,lag=lag,prewhite=FALSE,adjust=TRUE)
        se <- sqrt(unname(cv["rsv_mean01_5pp","rsv_mean01_5pp"]))
        data.frame(outcome=nm,period=period,trend_df=df,
          seasonal="month_fixed",HAC_lag=lag,beta=b,SE_HAC=se,
          IRR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),
          p_HAC=2*stats::pnorm(-abs(b/se)),
          dispersion=sum(pear^2)/stats::df.residual(fit),
          LjungBox_p12=stats::Box.test(pear,lag=12,type="Ljung-Box")$p.value)
      }))
    })
    if (!is.null(z)) national_stress[[length(national_stress)+1L]] <- z
  }
}
if (length(national_stress)) save_csv(dplyr::bind_rows(national_stress),
 "11_national_comparator_diagnostic_stress_tests.csv")
if(length(errors)) save_csv(dplyr::bind_rows(errors),"model_errors.csv")

# ------------------------- F. Figures and summary --------------------------
primary_plot <- dplyr::bind_rows(models)
p <- ggplot2::ggplot(primary_plot,ggplot2::aes(x=IRR,y=outcome)) +
  ggplot2::geom_vline(xintercept=1,linetype="dashed",color="grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low_t_approx,
                     xmax=CI_high_t_approx),width=.2,orientation="y") +
  ggplot2::geom_point(size=2) + ggplot2::facet_wrap(~period) +
  ggplot2::theme_bw(base_size=11)+
  ggplot2::labs(x="Regional IRR per +5 pp RSV PCR positivity (10-cluster approximate 95% CI)",
                y="Cause of death",title="Regional RSV-mortality associations")
ggplot2::ggsave(file.path(OUT,"figure_01_regional_forest.png"),p,
                width=10,height=5,dpi=250)
if(length(loo)) {
  lp <- dplyr::bind_rows(loo)
  q <- ggplot2::ggplot(lp,ggplot2::aes(x=omitted_region,y=IRR))+
    ggplot2::geom_hline(yintercept=1,linetype=2,color="grey50")+
    ggplot2::geom_point(size=2)+ggplot2::geom_line()+
    ggplot2::scale_x_continuous(breaks=1:10)+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::labs(x="HHS region excluded",y="Estimated MND IRR",
                  title="Leave-one-region-out influence")
  ggplot2::ggsave(file.path(OUT,"figure_02_MND_leave_one_region_out.png"),q,
                  width=7,height=4,dpi=250)
}
summary <- c(
  "REGIONAL RSV-MND VALIDATION: analysis summary",
  "==============================================",
  "Period: 2011-2019; RSV PCR NREVSS; deaths from CDC WONDER UCD.",
  "Exploratory follow-up to earlier national models, NOT preregistered.",
  paste("Common region-months:",nrow(common),"of 1080"),
  paste("Suppressed/absent MND region-months:",sum(is.na(panel$deaths_mnd))),
  paste("Regional MND (2011-2019) IRR:",signif(primary$IRR[primary$outcome=="mnd"],4)),
  paste("Regional MND (2013-2019) IRR:",signif(later$IRR[later$outcome=="mnd"],4)),
  paste("Bootstrap successful region-pairs draws:",nrow(boots),"of",B),
  "Inference caveat: 10 clusters. Cluster-t and pairs-bootstrap intervals are approximate.",
  "A null result does not establish equivalence or absence of an effect.",
  "Regional RSV denotes reporting LABORATORY region; deaths denote residence region.",
  "Fixed effects absorb region-year, region-specific month-of-year and national month-year.",
  "Underlying-cause mortality does not identify RSV-infected people or ALS progression.",
  "Suppressed WONDER cells have NOT been imputed or deduced from national totals.",
  "National stress tests are diagnostics, not post hoc model selection.",
  paste("Non-primary stage errors:",length(errors)),
  "Inspect coverage, identification, comparative models, bootstraps and diagnostics",
  "before making a publication decision."
)
writeLines(summary,file.path(OUT,"analysis_summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_session_info.txt"))
msg("COMPLETE. Read regional_validation_results/analysis_summary.txt first.")
