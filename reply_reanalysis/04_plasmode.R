# =====================================================================
#  04_plasmode.R — plasmode simulation: graded exposure and the effect sizes
#  the designs can detect or exclude
#
#  Null data-generating process (DGP): logistic model fitted to the real
#  analytic data with study-specific natural splines of birth year (df = 6;
#  keeps each study's observed smooth age/cohort gradient), sex, SEP and
#  study FE — and NO famine term. On top of the null probabilities p0 we
#  inject a known famine effect delta (absolute risk difference, pp) of
#  three shapes:
#    U  uniform    : +delta for everyone aged 1-5 in a famine year
#    G  graded     : proportional to famine years at ages 1-5,
#                    scaled so the mean among the exposed is delta
#    R  graded + recovery spill-over: as G plus half-strength exposure for
#                    the two post-famine 'recovery' years at ages 1-5
#                    (cohorts born just after the famine); contaminates the
#                    post-famine control
#  Every estimator is the identical closure used on the observed data.
# =====================================================================
source("01_common.R")
set.seed(20261109)
t0 <- Sys.time()
d  <- make_analytic(load_pooled(4)); n <- nrow(d)
dpub <- d   # ITS on the analytic sample (differs negligibly from the published sample)

# ---- estimators (7) -------------------------------------------------
EST <- list(
  RDiT       = make_rd(d),
  Donut1     = make_rd(d, donut = 1),
  RF         = make_rf(d),
  ITS_pub    = make_its(d, weighted = FALSE, trend = "common", centre = "famine_start"),
  ITS_wtd    = make_its(d, weighted = TRUE,  trend = "common", centre = "famine_start"),
  ITS_local  = make_its(d, weighted = TRUE,  win = 15, trend = "study", centre = "cutoff"),
  Dose       = make_dose(d, win = 15, df = 4))
mean_dose <- mean(d$dose15[d$exposed15 == 1])
scale <- c(RDiT = 1, Donut1 = 1, RF = 1, ITS_pub = 1, ITS_wtd = 1, ITS_local = 1, Dose = mean_dose)
run_all <- function(y) {
  r <- sapply(EST, function(f) f(y)); r["est", ] <- r["est", ] * scale; r["se", ] <- r["se", ] * scale; r
}
obs <- run_all(d$edent)
cat("Observed (pp):\n"); print(round(100 * obs, 3))

# ---- null DGP -------------------------------------------------------
S <- do.call(cbind, lapply(sort(unique(d$study)), function(k) {
  M <- matrix(0, n, 6); r <- d$study == k; M[r, ] <- ns(d$byear[r], df = 6); M }))
Xg <- cbind(1, S, cov_block(d)); Xg <- drop_collinear(Xg)
g0 <- glm.fit(Xg, d$edent, family = binomial())
p0 <- as.numeric(g0$fitted.values)
cat(sprintf("Null DGP fitted: mean p0 = %.4f (observed %.4f)\n", mean(p0), mean(d$edent)))

# ---- effect shapes --------------------------------------------------
ex <- d$exposed15 == 1
shape <- list(
  U = as.numeric(ex),
  G = d$dose15 / mean(d$dose15[ex]),
  R = { z <- d$dose15 + 0.5 * d$rec15; z / mean(z[ex]) })
cat(sprintf("Recovery spill-over: share of post-famine RF controls receiving effect under R = %.3f\n",
            mean(shape$R[d$cohort3 %in% 2 & d$treatment %in% 0] > 0)))

sim_cell <- function(delta, shp, R) {
  dp <- delta * shape[[shp]]
  out <- array(NA_real_, c(R, length(EST), 2), dimnames = list(NULL, names(EST), c("est", "se")))
  for (r in seq_len(R)) {
    y <- rbinom(n, 1, pmin(pmax(p0 + dp, 0), 1))
    out[r, , ] <- t(run_all(y))
  }
  out
}

R0 <- as.integer(Sys.getenv("PLASMODE_R0", 600)); R1 <- as.integer(Sys.getenv("PLASMODE_R1", 300))
deltas <- c(1, 2, 3, 5) / 100
null <- sim_cell(0, "U", R0)
cat("Null done:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cells <- list()
for (shp in c("U", "G", "R")) for (dl in deltas) {
  cells[[paste(shp, dl * 100)]] <- sim_cell(dl, shp, R1)
  cat(sprintf("  %s delta=%g done (%.1f min)\n", shp, dl * 100,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(obs = obs, null = null, cells = cells, mean_dose = mean_dose, R0 = R0, R1 = R1),
        file.path(OUT, "04_plasmode_raw.rds"))

# ---- summaries ------------------------------------------------------
summ <- function(a, delta, shp) {
  do.call(rbind, lapply(names(EST), function(e) {
    est <- a[, e, "est"]; se <- a[, e, "se"]; z <- est / se
    data.frame(shape = shp, delta_pp = 100 * delta, estimator = e,
               mean_est_pp = 100 * mean(est), sd_est_pp = 100 * sd(est), mean_se_pp = 100 * mean(se),
               power_two_sided = mean(abs(z) > zc), power_positive = mean(z > zc),
               p_le_observed = mean(est <= obs["est", e]))
  }))
}
tab <- rbind(summ(null, 0, "null"),
             do.call(rbind, lapply(names(cells), function(k) {
               s <- strsplit(k, " ")[[1]]; summ(cells[[k]], as.numeric(s[2]) / 100, s[1]) })))
b0 <- setNames(tab$mean_est_pp[tab$shape == "null"], tab$estimator[tab$shape == "null"])
tab$attenuation <- ifelse(tab$delta_pp > 0, (tab$mean_est_pp - b0[tab$estimator]) / tab$delta_pp, NA)
print(tab, digits = 3, row.names = FALSE)
wcsv(tab, "04_plasmode_summary.csv")

# joint: probability that ALL valid designs are at or below their observed values
valid <- c("RDiT", "RF", "ITS_local", "Dose")
joint <- do.call(rbind, c(list(data.frame(shape = "null", delta_pp = 0,
  p_all_le_obs = mean(apply(null[, valid, "est"], 1, function(v) all(v <= obs["est", valid])))) ),
  lapply(names(cells), function(k) { s <- strsplit(k, " ")[[1]]
    data.frame(shape = s[1], delta_pp = as.numeric(s[2]),
               p_all_le_obs = mean(apply(cells[[k]][, valid, "est"], 1, function(v) all(v <= obs["est", valid])))) })))
print(joint, row.names = FALSE); wcsv(joint, "04_plasmode_joint.csv")

# ---- GLS combination of the valid designs: see 04b_gls.R ------------
#  (separate script so it runs from the saved replications; the first
#   version of this block mis-aligned the null means and stopped here)
source(file.path(RS, "04b_gls.R"))
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\nDONE 04\n")
