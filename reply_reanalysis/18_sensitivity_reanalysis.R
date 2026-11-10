# =====================================================================
#  18_sensitivity_reanalysis.R — further sensitivity analyses on the real
#  data (run alongside the overnight simulation; single core)
#  (a) placebo cut-offs for the RDiT (cut-off moved to cohorts with no
#      famine exposure in the window): design-based calibration of the
#      RDiT standard error against real cohort features
#  (b) placebo famine dates for the dose-response analysis (famine shifted
#      by 9-20 years, real exposure adjusted for): calibration of the
#      dose-response standard error
#  (c) lack-of-fit of the local ITS (+/-15 cohorts, study-specific linear
#      trends) against study-specific cubic trends
#  (d) dose-response with alternative susceptibility windows (ages 1-8,
#      1-12, 3-8)
#  (e) graded severity DiD (CHARLS), leave one birth province out
#  (f) honest RDiT on the study- and covariate-adjusted outcome over a grid
#      of curvature bounds
#  (g) RDiT at the age-11/12 boundary (where susceptibility would end under
#      the flat-to-11 scenario) and at the age-15/16 boundary (placebo)
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages({ library(fixest); library(RDHonest) })
sink(file.path(OUT, "18_sensitivity_reanalysis.log"), split = TRUE)
d <- make_analytic(load_pooled(4))

# (a) ---------------------------------------------------------------
pc <- do.call(rbind, lapply(c(-35:-20, 19:35), function(cc) {
  dd <- d; dd$run <- d$run - cc; dd$D <- as.integer(dd$run >= 0)
  w <- abs(dd$run) <= H_BW; nl <- sum(w & dd$run < 0); nr <- sum(w & dd$run >= 0)
  if (min(nl, nr) < 1000) return(NULL)
  r <- tryCatch(make_rd(dd)(dd$edent), error = function(e) c(est = NA, se = NA))
  data.frame(cutoff_shift = cc, n_left = nl, n_right = nr, est_pp = 100 * r["est"], se_pp = 100 * r["se"])
}))
pc$t <- pc$est_pp / pc$se_pp
cat("\n(a) Placebo cut-offs (RDiT, h = 9.52):\n"); print(pc, digits = 3, row.names = FALSE)
pa <- data.frame(n_placebos = nrow(pc), sd_est_pp = sd(pc$est_pp), mean_se_pp = mean(pc$se_pp), sd_t = sd(pc$t),
                 share_abs_t_gt_1.96 = mean(abs(pc$t) > zc), actual_est_pp = 100 * make_rd(d)(d$edent)["est"])
print(pa, digits = 3); wcsv(pc, "18a_placebo_cutoffs.csv"); wcsv(pa, "18a_placebo_cutoffs_summary.csv")

# (b) ---------------------------------------------------------------
pd <- do.call(rbind, lapply(c(-20:-9, 9:20), function(k) {
  dd <- d; dd$famine_start <- d$famine_start + k; dd$famine_end <- d$famine_end + k
  dd$dose_pk <- overlap_years(dd$byear, dd$famine_start, dd$famine_end, 1, 5)
  dd$dose_real <- d$dose15; dd$dose615_real <- d$dose615
  r <- tryCatch(make_dose(dd, dose_var = "dose_pk", extra = c("dose_real", "dose615_real"))(dd$edent), error = function(e) c(est = NA, se = NA))
  data.frame(shift = k, est_pp_per_year = 100 * r["est"], se_pp = 100 * r["se"])
}))
pd$t <- pd$est_pp_per_year / pd$se_pp
cat("\n(b) Placebo famine dates (dose-response, per famine year at ages 1-5):\n"); print(pd, digits = 3, row.names = FALSE)
pb <- data.frame(n_placebos = sum(is.finite(pd$t)), sd_est = sd(pd$est_pp_per_year, na.rm = TRUE), mean_se = mean(pd$se_pp, na.rm = TRUE),
                 sd_t = sd(pd$t, na.rm = TRUE), share_abs_t_gt_1.96 = mean(abs(pd$t) > zc, na.rm = TRUE), actual = 100 * make_dose(d)(d$edent)["est"])
print(pb, digits = 3); wcsv(pd, "18b_placebo_famine_dates.csv"); wcsv(pb, "18b_placebo_famine_dates_summary.csv")

# (c) ---------------------------------------------------------------
cd <- d |> filter(abs(run) <= 15) |> group_by(study, run) |> summarise(ybar = mean(edent), n = n(), .groups = "drop") |>
  mutate(post = as.integer(run >= 0), tp = run * post)
m1 <- lm(ybar ~ post + factor(study) + factor(study):run + factor(study):tp, weights = n, data = cd)
m2 <- update(m1, . ~ . + factor(study):I(run^2) + factor(study):I(run^3))
m3 <- update(m1, . ~ . + factor(study):I(run^2) + factor(study):I(run^3) + factor(study):I(tp * run) + factor(study):I(tp * run^2))
lof <- data.frame(model = c("linear (as used)", "+ cubic trend", "+ cubic on each side"),
                  post_pp = 100 * c(coef(m1)["post"], coef(m2)["post"], coef(m3)["post"]),
                  se_pp = 100 * c(sqrt(sandwich::vcovHC(m1, type = "HC3")["post", "post"]), sqrt(sandwich::vcovHC(m2, type = "HC3")["post", "post"]),
                                  sqrt(sandwich::vcovHC(m3, type = "HC3")["post", "post"])))
a12 <- anova(m1, m2); a13 <- anova(m1, m3)
lof$F_vs_linear <- c(NA, a12$F[2], a13$F[2]); lof$df1 <- c(NA, a12$Df[2], a13$Df[2]); lof$df2 <- c(NA, a12$Res.Df[2], a13$Res.Df[2])
lof$p_vs_linear <- c(NA, a12$`Pr(>F)`[2], a13$`Pr(>F)`[2])
# reference: binomial residual scale (weighted residual variance x mean n should be about p(1-p) if only sampling noise)
lof$resid_disp <- c(sum(m1$weights * resid(m1)^2) / m1$df.residual, sum(m2$weights * resid(m2)^2) / m2$df.residual, sum(m3$weights * resid(m3)^2) / m3$df.residual) /
  weighted.mean(cd$ybar * (1 - cd$ybar), cd$n)
cat("\n(c) ITS lack of fit (+/-15 cohorts):\n"); print(lof, digits = 3, row.names = FALSE); wcsv(lof, "18c_its_lack_of_fit.csv")

# (d) ---------------------------------------------------------------
wins <- list(`1-5 (main)` = c(1, 5), `0-5` = c(0, 5), `1-8` = c(1, 8), `3-8` = c(3, 8), `1-12` = c(1, 12))
dw <- do.call(rbind, lapply(names(wins), function(nm) {
  a <- wins[[nm]]; dd <- d; dd$dosew <- overlap_years(d$byear, d$famine_start, d$famine_end, a[1], a[2])
  r <- make_dose(dd, dose_var = "dosew")(dd$edent); md <- mean(dd$dosew[dd$dosew >= 1])
  data.frame(window = nm, per_year_pp = 100 * r["est"], lb = 100 * (r["est"] - zc * r["se"]), ub = 100 * (r["est"] + zc * r["se"]),
             mean_dose_exposed = md, at_mean_pp = 100 * r["est"] * md, at_mean_ub = 100 * (r["est"] + zc * r["se"]) * md)
}))
cat("\n(d) Dose-response by susceptibility window:\n"); print(dw, digits = 3, row.names = FALSE); wcsv(dw, "18d_dose_windows.csv")

# (e) ---------------------------------------------------------------
ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning", Heilongjing = "Heilongjiang")
ly$prov <- ifelse(lab %in% names(std), std[lab], lab); map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
stopifnot(abs(ly$edr1960_abs[ly$prov == "Shanxi"] - 2.20) < 0.01, abs(ly$edr1960_abs[ly$prov == "Shaanxi"] - 1.90) < 0.01)
ch <- readRDS(file.path(LOC, "charls_did.rds")) |> filter(!is.na(edent), !is.na(age), !is.na(male), !is.na(birthprov)) |>
  mutate(prov = ifelse(birthprov %in% names(map), map[birthprov], birthprov)) |>
  left_join(ly[, c("prov", "edr1960_abs")], by = "prov") |> filter(!is.na(edr1960_abs)) |> mutate(I = edr1960_abs / 10)
Ih <- read.csv(file.path(OUT, "11_severity_implied_effects.csv"))$intensity[1]
s <- ch |> filter(byear %in% c(1954:1960, 1961:1967)) |> mutate(ex = as.integer(byear %in% 1954:1960))
loo <- do.call(rbind, lapply(c("(none)", sort(unique(s$prov))), function(p) {
  ss <- if (p == "(none)") s else s[s$prov != p, ]
  f <- feols(edent ~ I:ex | prov + byear, data = ss, cluster = ~prov)
  data.frame(left_out = p, n = nobs(f), per10_pp = 100 * coef(f)[1], se_pp = 100 * se(f)[1], implied_high_pp = 100 * coef(f)[1] * Ih / 10,
             implied_high_ub = 100 * (coef(f)[1] + zc * se(f)[1]) * Ih / 10)
}))
cat("\n(e) Graded severity DiD vs post-famine controls, leave one province out:\n"); print(loo, digits = 3, row.names = FALSE)
wcsv(loo, "18e_severity_did_leave_one_province_out.csv")

# (f) ---------------------------------------------------------------
yres <- resid(lm(edent ~ male + factor(father_ses) + factor(study), data = d))
dd <- data.frame(y = yres, run = d$run)
h0 <- RDHonest(y ~ run, data = dd, cutoff = 0, kern = "triangular", opt.criterion = "MSE", sclass = "H"); M0 <- as.data.frame(h0$coefficients)$M[1]
hg <- do.call(rbind, lapply(c(0.5, 1, 2, 4), function(mm) {
  f <- if (mm == 1) h0 else RDHonest(y ~ run, data = dd, cutoff = 0, kern = "triangular", opt.criterion = "MSE", sclass = "H", M = M0 * mm)
  c1 <- as.data.frame(f$coefficients)
  data.frame(M_multiple = mm, M = M0 * mm, est_pp = 100 * c1$estimate[1], lb = 100 * c1$conf.low[1], ub = 100 * c1$conf.high[1],
             ub_one_sided = 100 * c1$conf.high.onesided[1], bandwidth = c1$bandwidth[1])
}))
cat("\n(f) Honest RDiT, adjusted outcome, M grid:\n"); print(hg, digits = 3, row.names = FALSE); wcsv(hg, "18f_rdhonest_adjusted_M_grid.csv")

# (g) ---------------------------------------------------------------
bd <- do.call(rbind, lapply(list(list("age 11/12 boundary", 6), list("age 15/16 boundary (placebo)", 10)), function(L) {
  do.call(rbind, lapply(c(4, 5), function(h) {
    dd <- d; dd$run <- d$run + L[[2]]; dd$D <- as.integer(dd$run >= 0)
    r <- make_rd(dd, h = h)(dd$edent)
    # right side (D = 1) = cohorts aged <= 11 (or <= 15) at onset = the younger side; harm there predicts a positive jump
    data.frame(boundary = L[[1]], h = h, effect_pp = 100 * r["est"], lb = 100 * (r["est"] - zc * r["se"]), ub = 100 * (r["est"] + zc * r["se"]))
  }))
}))
cat("\n(g) RDiT at older-age boundaries (younger side minus older side; harm on the younger side predicts a positive value):\n")
print(bd, digits = 3, row.names = FALSE); wcsv(bd, "18g_older_age_boundaries.csv")
cat("DONE 18\n"); sink()
