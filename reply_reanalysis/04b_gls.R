# =====================================================================
#  04b_gls.R — combine the four valid designs (RDiT, reduced form, local
#  ITS, dose-response) into ONE estimate of a famine effect delta, for each
#  assumed effect shape, using the saved plasmode replications:
#    theta_obs - b0 = a * delta + e,   e ~ N(0, Sigma)
#  b0 = mean estimate under the null DGP (small curvature bias), Sigma =
#  covariance of the four estimators under the null (captures their
#  correlation: they share data), a = response of each estimator per unit
#  of true delta (slope over delta = 0, 1, 2, 3, 5 pp). GLS:
#    delta_hat = (a' S^-1 a)^-1 a' S^-1 (theta_obs - b0).
# =====================================================================
source("01_common.R")
z <- readRDS(file.path(OUT, "04_plasmode_raw.rds"))
valid <- c("RDiT", "RF", "ITS_local", "Dose")
deltas <- c(1, 2, 3, 5) / 100
b0  <- colMeans(z$null[, valid, "est"])
Sig <- cov(z$null[, valid, "est"]); Si <- solve(Sig)
th  <- z$obs["est", valid] - b0

gls <- do.call(rbind, lapply(c("U", "G", "R"), function(shp) {
  m <- sapply(deltas, function(dl) colMeans(z$cells[[paste(shp, dl * 100)]][, valid, "est"]))  # 4 x 4
  a <- sapply(seq_along(valid), function(i) coef(lm(c(b0[i], m[i, ]) ~ c(0, deltas)))[2])
  v <- 1 / as.numeric(t(a) %*% Si %*% a); est <- v * as.numeric(t(a) %*% Si %*% th); se <- sqrt(v)
  ind <- th / a; ind_se <- sqrt(diag(Sig)) / a
  data.frame(shape = c(U = "uniform", G = "graded (dose-proportional)", R = "graded + recovery spill-over")[shp],
             combined_pp = 100 * est, se_pp = 100 * se, lb = 100 * (est - zc * se), ub = 100 * (est + zc * se),
             p_two_sided_vs_0 = 2 * pnorm(-abs(est / se)),
             p_one_sided_delta_ge_1pp = pnorm((est - 0.01) / se),
             p_one_sided_delta_ge_2pp = pnorm((est - 0.02) / se),
             tost_p_2pp = tost_p(est, se, 0.02),
             a_RDiT = a[1], a_RF = a[2], a_ITS_local = a[3], a_Dose = a[4],
             implied_RDiT_pp = 100 * ind[1], implied_RF_pp = 100 * ind[2],
             implied_ITS_local_pp = 100 * ind[3], implied_Dose_pp = 100 * ind[4],
             implied_RDiT_se = 100 * ind_se[1], implied_RF_se = 100 * ind_se[2],
             implied_ITS_local_se = 100 * ind_se[3], implied_Dose_se = 100 * ind_se[4])
}))
print(gls, digits = 3, row.names = FALSE); wcsv(gls, "04_plasmode_gls.csv")
cat("\nNull-DGP correlation of the four valid designs:\n"); print(round(cor(z$null[, valid, "est"]), 2))
wcsv(as.data.frame(round(cor(z$null[, valid, "est"]), 3)), "04_plasmode_null_correlation.csv")
# Monte-Carlo upper bound for the joint 0/R1 result (exact binomial, one-sided 95%)
cat(sprintf("\nJoint check: 0 of %d replications per scenario -> one-sided 95%% upper bound %.2f%%\n",
            z$R1, 100 * (1 - 0.05^(1 / z$R1))))
cat("DONE 04b\n")
