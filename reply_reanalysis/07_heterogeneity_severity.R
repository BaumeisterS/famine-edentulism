# =====================================================================
#  07_heterogeneity_severity.R — pooling famines of very
#  different severity (about 400-1200 kcal/day) and survivor selection
#  (famine-to-survey interval up to ~89 years, Mexico).
#
#  Study-specific estimates for four designs (RDiT, reduced form, local
#  ITS, dose-response) -> random-effects meta-analysis (REML): tau^2, I^2,
#  prediction interval; univariable meta-regression on (i) documented
#  minimum daily energy intake as printed in the paper's Table 1, (ii)
#  famine duration, (iii) famine-to-survey interval (survivor-selection
#  proxy) and (iv) outcome type (tooth count vs binary self-report);
#  duration strata for the design-based estimators.
#  Between-famine meta-regression is ecological and has 8-10 data points:
#  exploratory. The within-country severity test is 06 (CHARLS).
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages(library(metafor))
d  <- make_analytic(load_pooled(4))

# study characteristics (paper Table 1 as printed; survey years from the Table 1 note)
char <- data.frame(study = 1:10, label = study_labels,
  kcal_floor = c(1000, NA, 1000, 1000, NA, 1000, 400, 1000, 1200, 1000),
  kcal_note  = c("<1000 in severely hit regions (1959-61)", "no Henan-specific figure",
                 "<1000 (labouring classes)", "<1000 (labouring classes)", "acute scarcity; no figure",
                 "~1000-1200 urban; rural lower", "400-800 (winter 1944)", "1000-1250", "<1200 urban", "1000-1500 rations"),
  survey_year = c(2013, 2008, 2017, 2003, 2003, 2003, 2013, 2013, 2013, 2013),
  tooth_count = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0))
fs <- tapply(d$famine_start, d$study, `[`, 1); dur <- tapply(d$famine_dur, d$study, `[`, 1)
char$famine_start <- as.numeric(fs[as.character(char$study)])
char$duration <- as.numeric(dur[as.character(char$study)])
char$interval <- char$survey_year - char$famine_start

# study-specific estimates
pub <- read.csv(file.path(GH, "tables", "rdd_by_study.csv"))
pub$se <- (pub$ub - pub$lb) / (2 * zc)
one <- function(k) {
  dd <- d[d$study == k, ]; md <- mean(dd$dose15[dd$exposed15 == 1])
  g <- function(f, sc = 1) { r <- tryCatch(f(dd$edent), error = function(e) c(est = NA, se = NA)); sc * r }
  rf <- g(make_rf(dd)); it <- g(make_its(dd, weighted = TRUE, win = 15, centre = "cutoff"))
  ds <- g(make_dose(dd, win = 15, df = 4), md)
  data.frame(study = k,
    RDiT_est = 100 * pub$est[pub$study == k], RDiT_se = 100 * pub$se[pub$study == k],
    RF_est = 100 * rf["est"], RF_se = 100 * rf["se"],
    ITS_est = 100 * it["est"], ITS_se = 100 * it["se"],
    Dose_est = 100 * ds["est"], Dose_se = 100 * ds["se"])
}
ss <- merge(char, do.call(rbind, lapply(1:10, one)), by = "study")
print(ss[, c("label","kcal_floor","duration","interval","RDiT_est","RDiT_se","RF_est","RF_se","ITS_est","ITS_se","Dose_est","Dose_se")],
      digits = 3, row.names = FALSE)
wcsv(ss, "07_study_specific_estimates.csv")

# random-effects meta-analysis + meta-regression
des <- c("RDiT", "RF", "ITS", "Dose")
ma <- do.call(rbind, lapply(des, function(z) {
  yi <- ss[[paste0(z, "_est")]]; sei <- ss[[paste0(z, "_se")]]; ok <- is.finite(yi) & is.finite(sei)
  m <- rma(yi = yi[ok], sei = sei[ok], method = "REML"); pr <- predict(m)
  data.frame(design = z, k = sum(ok), pooled_pp = m$b[1], lb = m$ci.lb, ub = m$ci.ub, tau2 = m$tau2,
             I2 = m$I2, Q = m$QE, Q_p = m$QEp, pi_lb = pr$pi.lb, pi_ub = pr$pi.ub)
}))
print(ma, digits = 3, row.names = FALSE); wcsv(ma, "07_random_effects.csv")

mods <- c("kcal_floor", "duration", "interval", "tooth_count")
mr <- do.call(rbind, lapply(des, function(z) do.call(rbind, lapply(mods, function(mo) {
  yi <- ss[[paste0(z, "_est")]]; sei <- ss[[paste0(z, "_se")]]; x <- ss[[mo]]
  ok <- is.finite(yi) & is.finite(sei) & is.finite(x)
  m <- tryCatch(rma(yi = yi[ok], sei = sei[ok], mods = ~ x[ok], method = "REML"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.frame(design = z, moderator = mo, k = sum(ok), slope = m$b[2], se = m$se[2],
             lb = m$ci.lb[2], ub = m$ci.ub[2], p = m$pval[2], R2 = m$R2)
}))))
mr$slope_unit <- ifelse(mr$moderator == "kcal_floor", "pp per 100 kcal/day", ifelse(mr$moderator == "interval", "pp per 10 years", ifelse(mr$moderator == "duration", "pp per famine year", "pp tooth-count vs binary")))
sc <- ifelse(mr$moderator == "kcal_floor", 100, ifelse(mr$moderator == "interval", 10, 1))
mr$slope <- mr$slope * sc; mr$se <- mr$se * sc; mr$lb <- mr$lb * sc; mr$ub <- mr$ub * sc
print(mr, digits = 3, row.names = FALSE); wcsv(mr, "07_meta_regression.csv")

# duration strata (short <= 2 y vs long) for the design-based estimators (pooled within stratum)
strat <- do.call(rbind, lapply(list(Short = which(char$duration <= 2), Long = which(char$duration > 2)), function(ix) {
  dd <- d[d$study %in% ix, ]; md <- mean(dd$dose15[dd$exposed15 == 1])
  r <- list(RDiT = make_rd(dd)(dd$edent), RF = make_rf(dd)(dd$edent),
            ITS_local = make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(dd$edent),
            Dose = make_dose(dd)(dd$edent) * md)
  data.frame(estimator = names(r), est_pp = 100 * sapply(r, `[`, 1), se_pp = 100 * sapply(r, `[`, 2),
             studies = paste(study_labels[ix], collapse = ", "))
}))
strat$stratum <- rep(c("Short (<= 2 y)", "Long (> 2 y)"), each = 4)
strat$lb <- strat$est_pp - zc * strat$se_pp; strat$ub <- strat$est_pp + zc * strat$se_pp
print(strat[, c("stratum","estimator","est_pp","se_pp","lb","ub")], digits = 3, row.names = FALSE)
wcsv(strat, "07_duration_strata_design_based.csv")

# leave-one-out for the design-based pooled estimators (incl. Mexico)
loo <- do.call(rbind, lapply(1:10, function(k) {
  dd <- d[d$study != k, ]; md <- mean(dd$dose15[dd$exposed15 == 1])
  data.frame(excluded = study_labels[k],
    RDiT = 100 * make_rd(dd)(dd$edent)["est"],
    ITS_local = 100 * make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(dd$edent)["est"],
    Dose = 100 * md * make_dose(dd)(dd$edent)["est"])
}))
print(loo, digits = 3, row.names = FALSE); wcsv(loo, "07_leave_one_out_design_based.csv")
cat("DONE 07\n")
