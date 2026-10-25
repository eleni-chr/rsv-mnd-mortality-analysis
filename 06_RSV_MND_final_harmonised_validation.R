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
# 06_RSV_MND_final_harmonised_validation.R
# Final targeted methodological validation of 2011-2019 RSV-PCR / MND mortality
# Exploratory follow-up to scripts 01-05. Version 1.0, 2026-10-24.
# Run from PROJECT ROOT: source("06_RSV_MND_final_harmonised_validation.R")
# Requires: fixest, sandwich, dplyr, ggplot2. It never modifies existing outputs.
# Inputs: regional_validation_results/03_common_panel_used.csv,
#         analysis_results/01_regional_rsv_monthly.csv,
#         analysis_results/04_national_matched_monthly.csv.
# Do not use fwildclusterboot::boottest on Poisson/fepois: that method handles
# linear/IV models, not these Poisson regressions.

START <- as.Date("2011-01-01")
END   <- as.Date("2019-12-01")
OUT   <- "final_method_validation_results"
# Inversion is computationally more demanding; TRUE is appropriate for final run.
RUN_SCORE_TEST_INVERSION <- TRUE
SCORE_GRID_POINTS <- 41L     # descriptive grid; cannot yield an exact CI endpoint
SCORE_GRID_HALF_WIDTH <- 0.15 # log-IRR, increased if 6*cluster SE is larger
SCORE_LEVEL <- 0.05
REGIONS_EXPECTED <- 1:10
OUTCOMES <- c("mnd", "parkinson", "transport", "respiratory")

pkgs <- c("fixest", "sandwich", "dplyr", "ggplot2")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop("Install packages first: install.packages(c(",
    paste(shQuote(missing_pkgs), collapse=", "), "))", call.=FALSE)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
LOG <- file.path(OUT, "analysis_log.txt")
if (file.exists(LOG)) file.remove(LOG)
logit <- function(...) {
  s <- paste0(..., collapse="")
  cat(s, "\n")
  cat(s, "\n", file=LOG, append=TRUE)
}
fail <- function(...) stop(paste0(..., collapse=""), call.=FALSE)
write_csv <- function(x, filename) utils::write.csv(x, file.path(OUT, filename),
                                             row.names=FALSE, na="")
get_file <- function(p, alt=NULL) {
  candidates <- c(p, alt)
  valid <- candidates[file.exists(candidates)]
  if (!length(valid)) fail("Required input missing: ", p,
                           ". Run scripts 04 v2 and 05 first, in project root.")
  valid[[1L]]
}
get_month <- function(x) {
  z <- as.Date(as.character(x))
  if (anyNA(z)) fail("Unparseable YYYY-MM-DD dates in input")
  z
}
all_months <- seq(START,END,by="month")
month_days <- function(x) {
  y <- as.integer(format(x,"%Y")); m <- as.integer(format(x,"%m"))
  nextmonth <- as.Date(sprintf("%04d-%02d-01", y+(m==12L), m%%12L+1L))
  as.integer(nextmonth - x)
}

# ---------------------- 1. Read previously audited inputs ------------------
common_file <- get_file("regional_validation_results/03_common_panel_used.csv",
                        "03_common_panel_used.csv")
rsv_file <- get_file("analysis_results/01_regional_rsv_monthly.csv",
                      "01_regional_rsv_monthly.csv")
national_file <- get_file("analysis_results/04_national_matched_monthly.csv",
                           "04_national_matched_monthly.csv")
regional <- utils::read.csv(common_file, stringsAsFactors=FALSE)
rsv <- utils::read.csv(rsv_file, stringsAsFactors=FALSE)
national <- utils::read.csv(national_file, stringsAsFactors=FALSE)
regional$month <- get_month(regional$month)
rsv$month <- get_month(rsv$month)
national$month <- get_month(national$month)
for (dname in c("regional","rsv")) {
  d <- get(dname)
  if (!all(c("region","month") %in% names(d)) ||
      anyDuplicated(paste(d$region,d$month)))
    fail(dname," has duplicate or missing region-month keys")
  if (!setequal(unique(as.integer(d$region)),REGIONS_EXPECTED))
    fail(dname," lacks some HHS regions")
}
if (anyDuplicated(national$month)) fail("National file has duplicate month keys")
required_regional <- c("region","month","rsv_mean01_5pp",
                       paste0("deaths_",OUTCOMES))
required_rsv <- c("region","month","rsv_mean01_5pp","tests_pcr")
required_national <- c("month",paste0("deaths_",OUTCOMES))
for (pair in list(list(regional,required_regional,"regional"),
                  list(rsv,required_rsv,"rsv"),
                  list(national,required_national,"national"))) {
  missing <- setdiff(pair[[2]],names(pair[[1]]))
  if(length(missing)) fail(pair[[3]]," missing columns: ",paste(missing,collapse=", "))
}
regional <- regional[regional$month>=START & regional$month<=END,,drop=FALSE]
rsv <- rsv[rsv$month>=START & rsv$month<=END,,drop=FALSE]
national <- national[national$month>=START & national$month<=END,,drop=FALSE]
if (nrow(rsv)!=1080L || nrow(national)!=108L) fail("Wrong study window in prior RSV/national outputs")
if (!all(regional$region %in% REGIONS_EXPECTED) || anyNA(regional$region)) fail("Bad region IDs")
if (anyNA(regional[,required_regional]) ||
    any(!is.finite(as.matrix(regional[,c("rsv_mean01_5pp",paste0("deaths_",OUTCOMES))]))))
  fail("Previously audited complete-case panel contains NA or invalid mortality counts")
if (any(!is.finite(rsv$tests_pcr)) || any(rsv$tests_pcr<=0) ||
    any(!is.finite(rsv$rsv_mean01_5pp)))
  fail("Missing, zero or nonfinite PCR exposure / test volume")
if (nrow(regional)!=1079L) {
  logit("NOTE: prior regional common case sample is ",nrow(regional),
        " rows, not 1079. Check prior analysis summary for changes.")
}
# Compare local predictor to source-file predictor, without silently changing it.
pred_compare <- merge(regional[,c("region","month","rsv_mean01_5pp")],
                      rsv[,c("region","month","rsv_mean01_5pp")],
                      by=c("region","month"), suffixes=c("_prior","_source"))
if (nrow(pred_compare)!=nrow(regional) ||
    max(abs(pred_compare$rsv_mean01_5pp_prior-
              pred_compare$rsv_mean01_5pp_source))>1e-8)
  fail("Regional RSV predictors differ from the source used in script 05")

# ---------------------- 2. Build exactly matched months -------------------
# One CDC-suppressed MND region-month existed in script 05. To compare 10-region
# and national estimates on identical calendar support, drop the WHOLE month
# containing any incomplete panel cell. NEVER infer or reconstruct that death.
count_by_month <- table(factor(format(regional$month,"%Y-%m"),
                            levels=format(all_months,"%Y-%m")))
months_complete <- all_months[as.integer(count_by_month)==10L]
months_dropped <- all_months[as.integer(count_by_month)!=10L]
if (length(months_complete)<96L) fail("Insufficient complete months for matched comparison")
if (length(months_dropped)>5L) fail("More than five months missing: investigate the export")
regional <- regional[regional$month %in% months_complete,,drop=FALSE]
national <- national[national$month %in% months_complete,,drop=FALSE]
rsv <- rsv[rsv$month %in% months_complete,,drop=FALSE]
if (nrow(regional)!=10L*length(months_complete) ||
    nrow(rsv)!=10L*length(months_complete) ||
    nrow(national)!=length(months_complete)) fail("Month harmonisation failed")
# Reconcile only fully observed months. No national-versus-partial differences
# are computed or written for the suppressed month.
for (nm in OUTCOMES) {
  col <- paste0("deaths_",nm)
  s <- stats::aggregate(regional[[col]], list(month=regional$month),sum)
  names(s)[2] <- "sum_complete"
  z <- merge(s,national[,c("month",col)],by="month",all=TRUE)
  if (nrow(z)!=length(months_complete) || anyNA(z) ||
      any(z$sum_complete != z[[col]]))
    fail("National/regional reconciliation failed for ",nm)
}
logit("Fully observed matched months: ",length(months_complete),
      "; region-months: ",nrow(regional),
      "; excluded incomplete calendar months: ",length(months_dropped))

# ---------------------- 3. Harmonise RSV PCR exposure ----------------------
# x_local: exact regional mean-of-current-and-lagged percent positivity / 5,
# from script 04, not a newly chosen lag specification.
# x_nat: contemporaneous PCR-test-weighted average of those SAME regional
# x_local variables. Hence x_local = x_nat + x_dev by construction.
# This new national definition is NOT guaranteed to equal the original pooled
# national positivity series: method-dependent contrast is measured explicitly.
rsv$rsv_month <- rsv$rsv_mean01_5pp
rsv$weighted_value <- rsv$rsv_month * rsv$tests_pcr
nat_rsv <- stats::aggregate(cbind(weighted_value,tests_pcr) ~ month,
                            data=rsv,FUN=sum)
nat_rsv$x_nat <- nat_rsv$weighted_value/nat_rsv$tests_pcr
nat_rsv <- nat_rsv[,c("month","x_nat","tests_pcr")]
regional <- merge(regional, nat_rsv[,c("month","x_nat")],by="month",sort=FALSE)
regional$x_local <- regional$rsv_mean01_5pp
regional$x_dev <- regional$x_local-regional$x_nat
regional$year <- as.integer(format(regional$month,"%Y"))
regional$moy <- as.integer(format(regional$month,"%m"))
regional$ym <- format(regional$month,"%Y-%m")
regional$region_year <- interaction(regional$region,regional$year,drop=TRUE)
regional$region_moy <- interaction(regional$region,regional$moy,drop=TRUE)
regional$log_days <- log(month_days(regional$month))
regional <- regional[order(regional$month,regional$region),,drop=FALSE]
national <- merge(national,nat_rsv[,c("month","x_nat")],by="month",sort=FALSE)
national$year <- factor(format(national$month,"%Y"))
national$moy <- factor(format(national$month,"%m"),levels=sprintf("%02d",1:12))
national$log_days <- log(month_days(national$month))
national <- national[order(national$month),,drop=FALSE]
if (max(abs(regional$x_local-regional$x_nat-regional$x_dev))>1e-10)
  fail("Exposure decomposition identity failed")
# Weighted deviations are zero by construction each month (using RSV tests).
check_d <- merge(regional[,c("month","region","x_dev")],
                 rsv[,c("month","region","tests_pcr")],
                 by=c("month","region"))
zero_dev <- stats::aggregate(I(x_dev*tests_pcr)~month,data=check_d,FUN=sum)
if (max(abs(zero_dev[[2]]))>1e-5) fail("Weighted regional deviations are not zero")
write_csv(data.frame(n_input_regional_months=108L,
                     n_matched_months=nrow(national),
                     n_matched_region_months=nrow(regional),
                     n_excluded_incomplete_months=length(months_dropped),
                     all_reconciliations_pass=TRUE,
                     max_local_exposure_identity_error=max(abs(regional$x_local-regional$x_nat-regional$x_dev)),
                     max_weighted_dev_sum= max(abs(zero_dev[[2]]))),
          "00_matched_sample_audit.csv")
# Withhold month-specific partial totals and any inferred suppressed cell.
write_csv(national[,c("month","x_nat")],"01_national_harmonised_exposure.csv")
write_csv(regional[,c("region","month","x_local","x_nat","x_dev")],
          "02_regional_exposure_decomposition.csv")
# Compare only with separately existing pooled monthly exposure (no outcome).
orig_nat <- utils::read.csv(national_file,stringsAsFactors=FALSE)
orig_nat$month <- get_month(orig_nat$month)
orig_nat <- orig_nat[orig_nat$month %in% months_complete,,drop=FALSE]
if (!"rsv_mean01_5pp" %in% names(orig_nat)) fail("Missing original national PCR positivity")
comparison <- merge(national[,c("month","x_nat")],
                    orig_nat[,c("month","rsv_mean01_5pp")],by="month")
comparison$difference <- comparison$x_nat - comparison$rsv_mean01_5pp
write_csv(comparison,"03_exposure_definition_comparison.csv")
logit("Exposure mean-of-local-versus-pooled national maximum absolute difference (per 5pp): ",
      signif(max(abs(comparison$difference)),4))

# ---------------------- 4. Model and inference helpers ---------------------
# HAC covers dependence across national calendar time (lag 6 primary);
# cluster-t and two-way cluster are approximate, with ONLY 10 regions.
safe_glm <- function(d,nm) {
  d$y <- d[[paste0("deaths_",nm)]]
  fit <- stats::glm(y ~ x_nat + year + moy + offset(log_days),
                    data=d,family=stats::poisson(link="log"),
                    control=stats::glm.control(maxit=100))
  if (!isTRUE(fit$converged) || stats::nobs(fit)!=nrow(d) ||
      !is.finite(stats::coef(fit)[["x_nat"]]))
    fail("National matched Poisson failed for ",nm)
  fit
}
get_nat <- function(f,d,nm,period,lags=c(3L,6L,12L)) {
  b <- unname(stats::coef(f)[["x_nat"]])
  do.call(rbind,lapply(lags,function(k) {
    v <- sandwich::NeweyWest(f,lag=k,prewhite=FALSE,adjust=TRUE)
    se <- sqrt(v["x_nat","x_nat"])
    if (!is.finite(se) || se<=0) fail("National HAC variance invalid")
    data.frame(period=period,outcome=nm,n_months=nrow(d),HAC_lag=k,
               beta=b,se=se,IRR=exp(b),CI_low=exp(b-1.96*se),
               CI_high=exp(b+1.96*se),
               p_normal=2*stats::pnorm(-abs(b/se)),
               poisson_pearson_dispersion=sum(stats::residuals(f,type="pearson")^2)/stats::df.residual(f),
               interpretation="National time-series; year and season FE; HAC approx")
  }))
}
fit_pois <- function(d,nm,type=c("within","joint","pooled_local")) {
  type <- match.arg(type)
  d$y <- d[[paste0("deaths_",nm)]]
  fm <- switch(type,
     within=y ~ x_local | region_year + region_moy + ym,
     joint=y ~ x_nat + x_dev | region_year + region_moy,
     pooled_local=y ~ x_local | region_year + region_moy)
  fit <- fixest::fepois(fm,data=d,offset=~log_days,nthreads=1,
                       fixef.tol=1e-9,notes=FALSE)
  must <- switch(type,within="x_local",joint=c("x_nat","x_dev"),
                 pooled_local="x_local")
  if(stats::nobs(fit)!=nrow(d) || !all(must %in% names(stats::coef(fit))) ||
     any(!is.finite(stats::coef(fit)[must])))
    fail("Regional Poisson fit dropped rows/coefficients: ",nm," / ",type)
  fit
}
get_reg <- function(f,d,nm,period,type) {
  vc_region <- stats::vcov(f,vcov=~region)
  # Two-way VCV accounts (asymptotically) for region and common month errors.
  # This is only diagnostic with G_region=10, not a small-cluster guarantee.
  vc_tw <- tryCatch(stats::vcov(f,vcov=~region+ym),error=function(e) NULL)
  G <- length(unique(d$region))
  ct <- stats::qt(.975,G-1)
  cfs <- stats::coef(f)
  out <- lapply(names(cfs),function(term) {
    b <- unname(cfs[[term]])
    se <- sqrt(unname(vc_region[term,term]))
    se2 <- if (!is.null(vc_tw) && is.finite(vc_tw[term,term]) &&
               vc_tw[term,term]>0) sqrt(unname(vc_tw[term,term])) else NA_real_
    if (!is.finite(se) || se<=0) fail("Invalid cluster variance: ",type," / ",term)
    data.frame(period=period,outcome=nm,model=type,term=term,n_region_months=nrow(d),
       n_regions=G,beta=b,IRR=exp(b),se_region_cluster=se,
       CI_low_region_t=exp(b-ct*se),CI_high_region_t=exp(b+ct*se),
       p_region_t=2*stats::pt(-abs(b/se),df=G-1),
       se_two_way_approx=se2,
       p_two_way_t9_approx=if(is.finite(se2)) 2*stats::pt(-abs(b/se2),G-1) else NA_real_)
  })
  do.call(rbind,out)
}

allnat <- list(); allreg <- list(); joint_tests <- list()
PERIODS <- list(`2011_2019`=2011L, `2013_2019`=2013L)
within_fits <- list()
for (period in names(PERIODS)) {
  lower_year <- PERIODS[[period]]
  dn <- national[as.integer(as.character(national$year))>=lower_year,,drop=FALSE]
  dr <- regional[regional$year>=lower_year,,drop=FALSE]
  dr$region_year <- droplevels(dr$region_year)
  dr$region_moy <- droplevels(dr$region_moy)
  logit("Fitting harmonised outcomes for ",period,"; months=",nrow(dn),
        "; regional records=",nrow(dr))
  for (nm in OUTCOMES) {
    fn <- safe_glm(dn,nm)
    allnat[[length(allnat)+1L]] <- get_nat(fn,dn,nm,period)
    for (model in c("within","joint","pooled_local")) {
      fp <- fit_pois(dr,nm,model)
      allreg[[length(allreg)+1L]] <- get_reg(fp,dr,nm,period,model)
      if (model=="within") within_fits[[paste(period,nm,sep="_")]] <- fp
      if(model=="joint") {
        v <- stats::vcov(fp,vcov=~region)
        b <- stats::coef(fp)
        # This is a test of NATIONAL versus REGIONAL-DEVIATION slopes INSIDE
        # the same panel joint model, NOT the difference of separate models.
        diff <- unname(b["x_nat"]-b["x_dev"])
        sed <- sqrt(unname(v["x_nat","x_nat"]+v["x_dev","x_dev"]-
                            2*v["x_nat","x_dev"]))
        joint_tests[[length(joint_tests)+1L]] <- data.frame(
          period=period,outcome=nm,beta_difference_nat_minus_dev=diff,
          ratio_of_IRRs=exp(diff),se_region_cluster=sed,
          p_t9_approx=2*stats::pt(-abs(diff/sed),9),
          description="Joint panel model within vs national component; 10-region inference approximate")
      }
    }
  }
}
national_results <- do.call(rbind,allnat)
regional_results <- do.call(rbind,allreg)
write_csv(national_results,"04_harmonised_national_models.csv")
write_csv(regional_results,"05_harmonised_regional_models.csv")
write_csv(do.call(rbind,joint_tests),"06_joint_model_coefficient_differences.csv")

# ---------------------- 5. Nonlinear wild CLUSTER SCORE ---------------------
# Restricted Poisson score sign-flip, NOT fwildclusterboot on fepois.
# For each null beta0, fit Poisson WITHOUT x_local, with offset beta0*x_local.
# Efficient score residualises x_local against ALL nuisance fixed effects
# using fitted Poisson means as WLS weights. Cluster influence scores are:
#   S_g = sum_{i in g} (y_i - mu0_i) * (x_i - P_W^mu x_i).
# The approximate null distribution comes from ALL 2^10 Rademacher
# signs on 10 region score contributions; no Monte Carlo noise.
# Sign patterns are exactly enumerated, but the P VALUES ARE NOT EXACT
# finite-sample tests: score symmetry / negligible leverage may fail with G=10.
# Crucially: NO pseudo-Poisson dependent variables, NO linear-model bootstrap.
all_signs <- function(G) {
  Z <- as.matrix(expand.grid(rep(list(c(-1,1)),G)))
  storage.mode(Z) <- "double"
  Z
}
W <- all_signs(10L)
score_test <- function(d,nm,beta0=0) {
  d$y <- d[[paste0("deaths_",nm)]]
  d$restricted_offset <- d$log_days + beta0*d$x_local
  null <- fixest::fepois(y ~ 1 | region_year + region_moy + ym,
                        data=d,offset=~restricted_offset,nthreads=1,
                        fixef.tol=1e-9,notes=FALSE)
  if(stats::nobs(null)!=nrow(d)) fail("Restricted null fit dropped observations")
  mu <- as.numeric(stats::predict(null,type="response"))
  if(length(mu)!=nrow(d) || any(!is.finite(mu)) || any(mu<=0))
    fail("Invalid restricted fitted means")
  d$mu0 <- mu
  # Weighted FE regression partials out the Poisson nuisance tangent space.
  proj <- fixest::feols(x_local ~ 1 | region_year + region_moy + ym,
                        data=d,weights=~mu0,nthreads=1,
                        fixef.tol=1e-9,notes=FALSE)
  if(stats::nobs(proj)!=nrow(d)) fail("Projection dropped observations")
  xr <- as.numeric(stats::residuals(proj))
  if(length(xr)!=nrow(d) || stats::sd(xr)<1e-5) fail("Insufficient residual exposure variation")
  cluster_score <- as.numeric(tapply((d$y-mu)*xr,
                                      factor(d$region,levels=REGIONS_EXPECTED),sum))
  if(length(cluster_score)!=10 || any(!is.finite(cluster_score)))
    fail("Invalid cluster scores")
  denom <- sqrt(sum(cluster_score^2))
  if(!is.finite(denom) || denom<1e-8) fail("Degenerate score test")
  observed_t <- sum(cluster_score)/denom
  null_stats <- as.vector(W%*%cluster_score)/denom
  # Exact enumeration of sign patterns, with floating-point tolerance.
  p <- mean(abs(null_stats)>=abs(observed_t)-1e-12)
  list(p=p,observed_t=observed_t,
       G=10L,n_sign_patterns=nrow(W),
       largest_abs_cluster_score_fraction=max(abs(cluster_score))/denom,
       score_total=sum(cluster_score),
       score_energy=denom,
       scores=cluster_score)
}
score_summary <- list(); score_contributions <- list(); score_grid <- list()
for (period in names(PERIODS)) {
  dr <- regional[regional$year>=PERIODS[[period]],,drop=FALSE]
  dr$region_year <- droplevels(dr$region_year)
  dr$region_moy <- droplevels(dr$region_moy)
  for(nm in OUTCOMES) {
    sc <- score_test(dr,nm,beta0=0)
    score_summary[[length(score_summary)+1L]] <- data.frame(
      period=period,outcome=nm,null_beta=0,
      wild_cluster_score_p=sc$p,observed_studentised_score=sc$observed_t,
      n_sign_patterns=sc$n_sign_patterns,n_clusters=sc$G,
      max_abs_cluster_score_fraction=sc$largest_abs_cluster_score_fraction,
      caveat="Approximate nonlinear cluster-score sign-flip; NOT exact finite-sample inference")
    score_contributions[[length(score_contributions)+1L]] <- data.frame(
      period=period,outcome=nm,region=REGIONS_EXPECTED,
      restricted_efficient_cluster_score=sc$scores)
    logit("Cluster-score ",period," / ",nm,": p=",format(sc$p,digits=5),
          " (1024 sign patterns).")
  }
  if (RUN_SCORE_TEST_INVERSION) {
    nm <- "mnd"
    fk <- within_fits[[paste(period,nm,sep="_")]]
    center <- unname(stats::coef(fk)[["x_local"]])
    vc <- stats::vcov(fk,vcov=~region)
    se <- sqrt(unname(vc["x_local","x_local"]))
    half <- max(SCORE_GRID_HALF_WIDTH,6*se)
    grid <- seq(center-half,center+half,length.out=SCORE_GRID_POINTS)
    for (j in seq_along(grid)) {
      sc <- score_test(dr,nm,beta0=grid[[j]])
      score_grid[[length(score_grid)+1L]] <- data.frame(
        period=period,outcome=nm,index=j,beta0=grid[[j]],IRR0=exp(grid[[j]]),
        p_score=sc$p,accept_at_5pct=sc$p>=SCORE_LEVEL)
      if (j%%10L==0L) logit("Score inversion ",period,": ",j,"/",length(grid))
    }
  }
}
score_results <- do.call(rbind,score_summary)
write_csv(score_results,"07_exactly_enumerated_wild_cluster_score_tests.csv")
write_csv(do.call(rbind,score_contributions),"08_null_cluster_score_contributions.csv")
if(length(score_grid)) {
  grid_df <- do.call(rbind,score_grid)
  write_csv(grid_df,"09_score_test_inversion_grid_MND.csv")
  gridsum <- do.call(rbind,lapply(split(grid_df,grid_df$period),function(x) {
    ok <- x[x$accept_at_5pct,,drop=FALSE]
    # Range of ACCEPTED GRID POINTS is not an exact CI; validity requires
    # the null score approximation and a dense/bracketed grid.
    data.frame(period=unique(x$period),n_grid_points=nrow(x),
      n_accepted=nrow(ok),accepted_IRR_min_grid=if(nrow(ok)) min(ok$IRR0) else NA_real_,
      accepted_IRR_max_grid=if(nrow(ok)) max(ok$IRR0) else NA_real_,
      touches_left_bound=isTRUE(x$accept_at_5pct[1]),
      touches_right_bound=isTRUE(tail(x$accept_at_5pct,1)),
      n_acceptance_segments=if(nrow(ok)) sum(diff(c(FALSE,x$accept_at_5pct))==1) else 0L,
      warning="GRID ACCEPTANCE RANGE ONLY, NOT AN EXACT 95% CI; check bounds/disconnected regions")
  }))
  write_csv(gridsum,"10_score_grid_acceptance_summary.csv")
}

# ---------------------- 6. Matched national/regional figure ---------------
plot_rows <- list()
for(period in names(PERIODS)) {
  nr <- national_results[national_results$period==period &
                         national_results$outcome=="mnd" &
                         national_results$HAC_lag==6L,,drop=FALSE]
  rr <- regional_results[regional_results$period==period &
                         regional_results$outcome=="mnd",,drop=FALSE]
  for(term in c("x_nat","x_dev")) {
    ii <- rr[rr$model=="joint" & rr$term==term,,drop=FALSE]
    plot_rows[[length(plot_rows)+1L]] <- data.frame(period=period,
      analysis=paste0("Joint panel: ",if(term=="x_nat") "national temporal" else "regional deviation"),
      IRR=ii$IRR,low=ii$CI_low_region_t,high=ii$CI_high_region_t)
  }
  w <- rr[rr$model=="within" & rr$term=="x_local",,drop=FALSE]
  plot_rows[[length(plot_rows)+1L]] <- data.frame(period=period,
    analysis="Panel with nationwide year-month FE",IRR=w$IRR,
    low=w$CI_low_region_t,high=w$CI_high_region_t)
  plot_rows[[length(plot_rows)+1L]] <- data.frame(period=period,
    analysis="National monthly series, HAC(6)",IRR=nr$IRR,
    low=nr$CI_low,high=nr$CI_high)
}
plt <- do.call(rbind,plot_rows)
plt$analysis <- factor(plt$analysis,levels=c("National monthly series, HAC(6)",
  "Joint panel: national temporal","Joint panel: regional deviation",
  "Panel with nationwide year-month FE"))
write_csv(plt,"11_figure_data_harmonised_IRRs.csv")
g <- ggplot2::ggplot(plt,ggplot2::aes(x=IRR,y=analysis)) +
    ggplot2::geom_vline(xintercept=1,linetype="dashed",colour="grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin=low,xmax=high),width=0.15,orientation="y") +
    ggplot2::geom_point(size=2.4) + ggplot2::facet_wrap(~period,ncol=2) +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::labs(title="Matched national and regional RSV–MND contrasts",
      subtitle="Common 10-region months; different estimands; intervals approximate",
      x="IRR per +5 percentage points of RSV PCR positivity",y=NULL)
ggplot2::ggsave(file.path(OUT,"figure_01_harmonised_national_vs_regional.png"),
                g,width=11.5,height=5.0,dpi=250)

# ---------------------- 7. Final log and transparency ---------------------
resids <- do.call(rbind,lapply(names(PERIODS),function(period) {
  fk <- within_fits[[paste(period,"mnd",sep="_")]]
  data.frame(period=period,
    n_obs=stats::nobs(fk),
    MND_within_beta=unname(stats::coef(fk)[["x_local"]]),
    MND_within_IRR=exp(unname(stats::coef(fk)[["x_local"]])),
    score_p_zero=score_results$wild_cluster_score_p[
       score_results$period==period & score_results$outcome=="mnd"])
}))
write_csv(resids,"12_MND_final_key_results.csv")
summary <- c(
  "FINAL RSV–MND METHODOLOGICAL VALIDATION (script 06)",
  "=================================================",
  "Study remains exploratory; 2011-2019 and 2013-2019 are previously used periods.",
  paste0("Identically observed calendar months: ", length(months_complete),"/108"),
  paste0("Region-month observations per outcome: ",nrow(regional)),
  paste0("Excluded complete calendar months due to at least one suppressed MND cell: ",
         length(months_dropped)),
  "No suppressed cell was imputed, recovered or reported by subtraction.",
  "Exposure is test-volume-weighted mean of SAME regional lag-01 PCR series.",
  "National and region-level predictors differ from original pooled national metric; see 03 CSV.",
  "National Poisson with year and month-of-year effects; HAC 3/6/12, focus HAC 6.",
  "Joint panel: common national PCR series + local-minus-national deviation; region-year and region-calendar-month FE.",
  "Within panel: region-year, region-calendar-month, and full national year-month FE.",
  "Joint model national and within coefficients not causally interpretable; common national shocks unmeasured.",
  "Nonlinear restricted Poisson efficient-score wild sign flips: all 1024 patterns, G=10 regions.",
  "EXACT enumeration of multipliers does NOT make this an exact finite-sample hypothesis test.",
  "The score test uses Poisson-null fitted means and weighted nuisance-FE orthogonalisation.",
  "Score-grid acceptance ranges, if generated, are approximate and may be bounded/disconnected.",
  "Ten clusters remains an inference limitation even if coefficients and p-values agree.",
  "Reporting-laboratory region differs from mortality residence region (ecological exposure misclassification).",
  "Do not infer causality or biological protection; non-significance does not establish no association.",
  "Do not use the numerical results from exploratory sensitivity choices to make post hoc confirmations.",
  "IMPORTANT: compare model estimands, not merely opposite signs or p-values.",
  "Primary outputs: 04, 05, 06, 07, 09, 10, 12 plus figure_01."
)
writeLines(summary,file.path(OUT,"analysis_summary.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_session_info.txt"))
logit("COMPLETE: outputs saved to ",OUT)
