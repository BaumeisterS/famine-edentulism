# =====================================================================
#  17b_plasmode_rd_summary.R — summarise the overnight plasmode (17) from
#  whatever blocks exist in output/17_chunks/ (can be run on partial
#  results). Per cell x estimator: mean, SD, mean SE, Monte Carlo SEs,
#  response per unit effect, rejection of 0, coverage of the true mean
#  effect, and the share of replications whose upper 95% limit lies below
#  +2 pp ("excludes 2 pp"). GLS combination with the conventional, robust
#  bias-corrected or honest RDiT, under the base and the cohort-shock DGP.
# =====================================================================
source("01_common.R")
CH <- file.path(OUT, "17_chunks")
fs <- list.files(CH, pattern = "\\.rds$", full.names = TRUE)
cellname <- sub("_[0-9]{4}-[0-9]{4}\\.rds$", "", basename(fs))
dat <- lapply(split(fs, cellname), function(f) do.call(rbind, lapply(f, readRDS)))
ob <- readRDS(file.path(OUT, "17_observed.rds")); obs <- ob$obs; TM <- ob$true_mean
SC <- data.frame(cell = c("null", "U2", "G2", "R2", "S1_2", "S2_2", "shock_null", "shock_G2", "G_RR2", "G_nonCN2", "G1", "S2_3"),
                 dgp = c("base", "base", "base", "base", "base", "base", "shock", "shock", "base", "base", "base", "base"),
                 shape = c("U", "U", "G", "R", "S1", "S2", "U", "G", "G_RR", "G_nonCN", "G", "S2"),
                 delta = c(0, .02, .02, .02, .02, .02, 0, .02, .02, .02, .01, .03), stringsAsFactors = FALSE)
SC <- SC[SC$cell %in% names(dat), ]
cat("Replications per cell:\n"); print(sapply(dat[SC$cell], nrow))

EST <- c(paste0(c("RDiT", "RDiT_h7", "RDiT_h5", "RDiT_h384", "RDiT_h3", "Donut1", "RF", "ITS_local", "Dose")),
         "rdr_conv", "rdr_robust", "hon_raw", "hon_dm")
get_e <- function(m, e) {
  if (e == "rdr_conv") { est <- m[, "rdr_conv_est"]; se <- m[, "rdr_conv_se"]; lb <- est - zc * se; ub <- est + zc * se; ub1 <- est + qnorm(.95) * se }
  else if (e == "rdr_robust") { est <- m[, "rdr_bc_est"]; se <- m[, "rdr_rob_se"]; lb <- est - zc * se; ub <- est + zc * se; ub1 <- est + qnorm(.95) * se }
  else if (grepl("^hon_", e)) { est <- m[, paste0(e, "_est")]; se <- m[, paste0(e, "_se")]; lb <- m[, paste0(e, "_lb")]; ub <- m[, paste0(e, "_ub")]; ub1 <- m[, paste0(e, "_ub1")] }
  else { est <- m[, paste0(e, "_est")]; se <- m[, paste0(e, "_se")]; lb <- est - zc * se; ub <- est + zc * se; ub1 <- est + qnorm(.95) * se }
  list(est = est, se = se, lb = lb, ub = ub, ub1 = ub1)
}
nullof <- function(dgp) if (dgp == "shock" && "shock_null" %in% names(dat)) "shock_null" else "null"
tab <- do.call(rbind, lapply(seq_len(nrow(SC)), function(i) {
  s <- SC[i, ]; m <- dat[[s$cell]]; m0 <- dat[[nullof(s$dgp)]]
  do.call(rbind, lapply(EST, function(e) {
    x <- get_e(m, e); x0 <- get_e(m0, e); ok <- is.finite(x$est); R <- sum(ok)
    tru <- if (s$delta > 0) s$delta * TM[[s$shape]] else 0
    pr <- function(v) c(mean(v[ok]), sqrt(mean(v[ok]) * (1 - mean(v[ok])) / R))
    rej <- pr(x$lb > 0 | x$ub < 0); exc2 <- pr(x$ub < 0.02); exc2_1 <- pr(x$ub1 < 0.02); cov <- pr(x$lb <= tru & x$ub >= tru)
    data.frame(cell = s$cell, dgp = s$dgp, shape = s$shape, delta_pp = 100 * s$delta, true_pooled_pp = 100 * tru, estimator = e, R = R,
               mean_pp = 100 * mean(x$est[ok]), mcse_mean_pp = 100 * sd(x$est[ok]) / sqrt(R), sd_pp = 100 * sd(x$est[ok]),
               mean_se_pp = 100 * mean(x$se[ok]), mean_ci_width_pp = 100 * mean((x$ub - x$lb)[ok]),
               response = if (s$delta > 0) (mean(x$est[ok]) - mean(x0$est, na.rm = TRUE)) / s$delta else NA,
               reject0 = rej[1], reject0_mcse = rej[2], excl2_twosided = exc2[1], excl2_onesided = exc2_1[1], excl2_mcse = exc2[2],
               coverage_true = cov[1], coverage_mcse = cov[2],
               mean_h = if (e == "rdr_conv" || e == "rdr_robust") mean(m[ok, "rdr_h"]) else if (grepl("^hon_", e)) mean(m[ok, paste0(e, "_h")]) else NA)
  }))
}))
print(tab[, c("cell", "estimator", "R", "mean_pp", "sd_pp", "mean_se_pp", "response", "reject0", "excl2_twosided", "excl2_onesided", "coverage_true", "mean_h")],
      digits = 3, row.names = FALSE)
wcsv(tab, "17_plasmode_rd_summary.csv")

# observed values on the real data (pp)
ob_tab <- data.frame(estimator = EST, t(sapply(EST, function(e) { x <- get_e(t(as.matrix(obs)), e); 100 * c(est = x$est, lb = x$lb, ub = x$ub, ub1 = x$ub1) })))
print(ob_tab, digits = 3, row.names = FALSE); wcsv(ob_tab, "17_observed_estimators.csv")

# GLS combinations
gls_fit <- function(set, shp_cells, null_cell, infl = NULL, label = "") {
  m0 <- dat[[null_cell]]; E0 <- sapply(set, function(e) get_e(m0, e)$est); ok0 <- stats::complete.cases(E0); E0 <- E0[ok0, , drop = FALSE]
  b0 <- colMeans(E0); Sig <- cov(E0)
  if (!is.null(infl)) { f <- ifelse(set %in% names(infl), infl[set], 1); Sig <- Sig * outer(f, f) }   # SD inflation from real-data placebos
  Si <- solve(Sig)
  th <- sapply(set, function(e) get_e(t(as.matrix(obs)), e)$est) - b0
  do.call(rbind, lapply(shp_cells, function(cc) {
    s <- SC[SC$cell == cc, ]; E1 <- sapply(set, function(e) get_e(dat[[cc]], e)$est)
    a <- (colMeans(E1, na.rm = TRUE) - b0) / s$delta
    v <- 1 / as.numeric(t(a) %*% Si %*% a); est <- v * as.numeric(t(a) %*% Si %*% th); se <- sqrt(v)
    r <- th - a * est; J <- as.numeric(t(r) %*% Si %*% r)
    data.frame(set = paste0(paste(set, collapse = "+"), label), null = null_cell, shape_cell = cc, R_null = nrow(E0), R_cell = nrow(E1),
               combined_pp = 100 * est, lb = 100 * (est - zc * se), ub = 100 * (est + zc * se),
               p_ge1 = pnorm((est - 0.01) / se), p_ge2 = pnorm((est - 0.02) / se), J = J, J_p = 1 - pchisq(J, length(set) - 1),
               weights = paste(sprintf("%.2f", v * as.numeric(t(a) %*% Si)), collapse = " / "),
               responses = paste(sprintf("%.2f", a), collapse = " / "))
  }))
}
base_cells <- intersect(c("U2", "G2", "R2", "S1_2", "S2_2", "G_RR2", "G_nonCN2"), names(dat))
sets <- list(c("RDiT", "RF", "ITS_local", "Dose"), c("rdr_robust", "RF", "ITS_local", "Dose"), c("hon_dm", "RF", "ITS_local", "Dose"),
             c("RF", "ITS_local", "Dose"), c("RDiT", "RF", "ITS_local"))
gls <- do.call(rbind, lapply(sets, function(st) tryCatch(gls_fit(st, base_cells, "null"), error = function(e) NULL)))
# variances recalibrated to the real-data placebo dispersion (SD of placebo t: RDiT 1.13 [18a], ITS 1.05 [18h],
# dose-response 1.34 [18b]; reduced form has no placebo analogue and gets the RDiT factor)
pf <- function(f, col, default) { x <- tryCatch(read.csv(file.path(OUT, f))[[col]][1], error = function(e) NA); if (is.finite(x)) max(1, x) else default }
INFL <- c(RDiT = pf("18a_placebo_cutoffs_summary.csv", "sd_t", 1.13), RF = pf("18a_placebo_cutoffs_summary.csv", "sd_t", 1.13),
          ITS_local = pf("18h_its_placebo_cutoffs_summary.csv", "sd_t", 1.05), Dose = pf("18b_placebo_famine_dates_summary.csv", "sd_t", 1.34))
INFL["rdr_robust"] <- INFL["RDiT"]; INFL["hon_dm"] <- 1
gls <- rbind(gls, do.call(rbind, lapply(sets[1:3], function(st) tryCatch(gls_fit(st, base_cells, "null", infl = INFL, label = " [placebo-calibrated variances]"),
                                                                        error = function(e) NULL))))
if (all(c("shock_null", "shock_G2") %in% names(dat)))
  gls <- rbind(gls, do.call(rbind, lapply(sets[1:3], function(st) tryCatch(gls_fit(st, "shock_G2", "shock_null"), error = function(e) NULL))))
print(gls, digits = 3, row.names = FALSE); wcsv(gls, "17_plasmode_rd_gls.csv")
cat("DONE 17b\n")
