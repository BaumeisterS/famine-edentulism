# =====================================================================
#  16_gls_sensitivity.R — how fragile is the GLS combination (04b/04c)?
#  For each effect shape: (i) all four designs, (ii) leaving out each
#  design in turn, (iii) all variances and covariances inflated by 1.5
#  (allowing for cohort-level shocks the plasmode DGP does not contain).
#  Reports the combined estimate, 95% CI and one-sided p for >= 1 and
#  >= 2 pp.
# =====================================================================
source("01_common.R")
z <- readRDS(file.path(OUT, "04_plasmode_raw.rds"))
valid <- c("RDiT", "RF", "ITS_local", "Dose")
b0 <- colMeans(z$null[, valid, "est"]); Sig0 <- cov(z$null[, valid, "est"]); th0 <- z$obs["est", valid] - b0
g1 <- read.csv(file.path(OUT, "04_plasmode_gls.csv")); g2 <- read.csv(file.path(OUT, "04c_susceptibility_gls.csv"))
A <- rbind(as.matrix(g1[, c("a_RDiT", "a_RF", "a_ITS_local", "a_Dose")]), as.matrix(g2[, c("a_RDiT", "a_RF", "a_ITS_local", "a_Dose")]))
shapes <- c(g1$shape, g2$shape)
comb <- function(a, th, S) { Si <- solve(S); v <- 1 / as.numeric(t(a) %*% Si %*% a); est <- v * as.numeric(t(a) %*% Si %*% th)
  c(est = 100 * est, lb = 100 * (est - zc * sqrt(v)), ub = 100 * (est + zc * sqrt(v)),
    p_ge1 = pnorm((est - 0.01) / sqrt(v)), p_ge2 = pnorm((est - 0.02) / sqrt(v))) }
res <- do.call(rbind, lapply(seq_along(shapes), function(i) {
  a <- A[i, ]
  vars <- list(`all four` = 1:4, `without RDiT` = 2:4, `without reduced form` = c(1, 3, 4),
               `without ITS` = c(1, 2, 4), `without dose-response` = 1:3)
  out <- do.call(rbind, lapply(names(vars), function(nm) { k <- vars[[nm]]
    data.frame(shape = shapes[i], variant = nm, t(comb(a[k], th0[k], Sig0[k, k, drop = FALSE]))) }))
  rbind(out, data.frame(shape = shapes[i], variant = "all four, variances x 1.5", t(comb(a, th0, 1.5 * Sig0))))
}))
print(res, digits = 3, row.names = FALSE); wcsv(res, "16_gls_sensitivity.csv")
cat("DONE 16\n")
