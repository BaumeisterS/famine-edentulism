# =====================================================================
#  04c_plasmode_susceptibility.R — susceptibility beyond age 5:
#  susceptibility to famine does not stop at age 5 (permanent-tooth crowns
#  mineralise until about age 8; canines, premolars and second molars erupt
#  into adolescence), so children just LEFT of the RDiT cut-off (aged 6+ at
#  onset) are also harmed and the sharp contrast is diluted.
#    S1 declining: weight 1 at ages 0-5, falling linearly to 0 at age 12
#    S2 flat     : weight 1 at every age 0-11
#  effect_i = beta x sum over famine years of weight(age in that year),
#  scaled so that the mean effect among the exposed (aged 1-5) is delta.
#  Same null DGP and estimators as 04_plasmode.R (identical p0).
# =====================================================================
source("01_common.R")
set.seed(20261110)
t0 <- Sys.time()
d <- make_analytic(load_pooled(4)); n <- nrow(d)
EST <- list(RDiT = make_rd(d), RF = make_rf(d),
            ITS_local = make_its(d, weighted = TRUE, win = 15, trend = "study", centre = "cutoff"),
            Dose = make_dose(d, win = 15, df = 4))
mean_dose <- mean(d$dose15[d$exposed15 == 1]); scale <- c(RDiT = 1, RF = 1, ITS_local = 1, Dose = mean_dose)
run_all <- function(y) { r <- sapply(EST, function(f) f(y)); r["est", ] <- r["est", ] * scale; r["se", ] <- r["se", ] * scale; r }
S <- do.call(cbind, lapply(sort(unique(d$study)), function(k) { M <- matrix(0, n, 6); r <- d$study == k; M[r, ] <- ns(d$byear[r], df = 6); M }))
Xg <- drop_collinear(cbind(1, S, cov_block(d))); p0 <- as.numeric(glm.fit(Xg, d$edent, family = binomial())$fitted.values)

wdose <- function(wfun) {   # weighted famine years over ages 0-15
  out <- numeric(n)
  for (a in 0:15) out <- out + wfun(a) * overlap_years(d$byear, d$famine_start, d$famine_end, a, a)
  out
}
ex <- d$exposed15 == 1
shapes <- list(S1_declining = wdose(function(a) ifelse(a <= 5, 1, pmax(0, (12 - a) / 7))),
               S2_flat_0_11 = wdose(function(a) as.numeric(a <= 11)))
shapes <- lapply(shapes, function(z) z / mean(z[ex]))
cat(sprintf("Share of the pre-famine cohort (left of cut-off, within h) receiving an effect: S1 %.2f, S2 %.2f\n",
            mean(shapes$S1_declining[d$run < 0 & d$run >= -H_BW] > 0), mean(shapes$S2_flat_0_11[d$run < 0 & d$run >= -H_BW] > 0)))
obs <- run_all(d$edent)
R1 <- as.integer(Sys.getenv("PLASMODE_R1", 200)); deltas <- c(1, 2, 3) / 100
cells <- list()
for (sh in names(shapes)) for (dl in deltas) {
  a <- array(NA_real_, c(R1, length(EST), 2), dimnames = list(NULL, names(EST), c("est", "se")))
  for (r in seq_len(R1)) a[r, , ] <- t(run_all(rbinom(n, 1, pmin(pmax(p0 + dl * shapes[[sh]], 0), 1))))
  cells[[paste(sh, dl * 100)]] <- a
  cat(sprintf("  %s delta=%g done (%.1f min)\n", sh, dl * 100, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
z <- readRDS(file.path(OUT, "04_plasmode_raw.rds"))      # null replications (same DGP) for bias and covariance
b0 <- colMeans(z$null[, names(EST), "est"]); Sig <- cov(z$null[, names(EST), "est"]); Si <- solve(Sig)
summ <- do.call(rbind, lapply(names(cells), function(k) { s <- strsplit(k, " ")[[1]]
  do.call(rbind, lapply(names(EST), function(e) { est <- cells[[k]][, e, "est"]; se <- cells[[k]][, e, "se"]
    data.frame(shape = s[1], delta_pp = as.numeric(s[2]), estimator = e, mean_est_pp = 100 * mean(est),
               response = (mean(est) - b0[e]) / (as.numeric(s[2]) / 100), power = mean(abs(est / se) > zc),
               p_le_obs = mean(est <= obs["est", e])) })) }))
print(summ, digits = 3, row.names = FALSE); wcsv(summ, "04c_susceptibility_summary.csv")
joint <- do.call(rbind, lapply(names(cells), function(k) data.frame(cell = k,
  p_all_le_obs = mean(apply(cells[[k]][, , "est"], 1, function(v) all(v <= obs["est", ]))))))
print(joint, row.names = FALSE); wcsv(joint, "04c_susceptibility_joint.csv")
th <- obs["est", ] - b0
gls <- do.call(rbind, lapply(names(shapes), function(sh) {
  m <- sapply(deltas, function(dl) colMeans(cells[[paste(sh, dl * 100)]][, , "est"]))
  a <- sapply(seq_along(EST), function(i) coef(lm(c(b0[i], m[i, ]) ~ c(0, deltas)))[2])
  v <- 1 / as.numeric(t(a) %*% Si %*% a); est <- v * as.numeric(t(a) %*% Si %*% th); se <- sqrt(v)
  data.frame(shape = sh, combined_pp = 100 * est, lb = 100 * (est - zc * se), ub = 100 * (est + zc * se),
             p_ge_1pp = pnorm((est - 0.01) / se), p_ge_2pp = pnorm((est - 0.02) / se),
             a_RDiT = a[1], a_RF = a[2], a_ITS_local = a[3], a_Dose = a[4]) }))
print(gls, digits = 3, row.names = FALSE); wcsv(gls, "04c_susceptibility_gls.csv")
cat("Total:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\nDONE 04c\n")
