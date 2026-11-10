# =====================================================================
#  10_rd_inference.R — inference for the primary RDiT
#  The published headline interval (+0.13; -1.50 to +1.77) equals the
#  conventional rdrobust interval (0.134 +/- 1.96 x 0.833). Here: (i) rdrobust
#  conventional / bias-corrected / robust estimates and intervals for the
#  published specification; (ii) RDHonest honest (bias-aware) intervals at
#  the estimated curvature bound M and over a grid of M, with ONE-SIDED
#  95% upper limits (the relevant test for 'can an effect >= 2 pp be
#  excluded').
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages({ library(rdrobust); library(RDHonest) })
d <- make_analytic(load_pooled(4))
cm <- model.matrix(~ male + factor(father_ses) + factor(study), data = d)[, -1]
t0 <- Sys.time()
r <- rdrobust(d$edent, d$run, c = 0, p = 1, kernel = "triangular", bwselect = "mserd", covs = cm, masspoints = "adjust")
rr <- data.frame(inference = rownames(r$coef), est_pp = 100 * r$coef[, 1], se_pp = 100 * r$se[, 1],
                 lb = 100 * r$ci[, 1], ub = 100 * r$ci[, 2], p = r$pv[, 1], h_left = r$bws[1, 1], h_right = r$bws[1, 2])
print(rr, digits = 3, row.names = FALSE); wcsv(rr, "10_rdrobust_inference.csv")
cat("rdrobust time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

dd <- d[, c("edent", "run")]
h0 <- RDHonest(edent ~ run, data = dd, cutoff = 0, kern = "triangular", opt.criterion = "MSE", sclass = "H")
cat("\nRDHonest default fit:\n"); print(h0)
co <- as.data.frame(h0$coefficients); cat("columns:", paste(names(co), collapse = ", "), "\n")
Mhat <- if (!is.null(co$M)) co$M[1] else if (!is.null(h0$M)) h0$M else NA
cat("Estimated M:", Mhat, "\n")
grid <- if (is.finite(Mhat)) Mhat * c(0.5, 1, 2, 4) else numeric(0)
hon <- do.call(rbind, lapply(c(NA, grid), function(M) {
  f <- if (is.na(M)) h0 else RDHonest(edent ~ run, data = dd, cutoff = 0, kern = "triangular",
                                      opt.criterion = "MSE", sclass = "H", M = M)
  c1 <- as.data.frame(f$coefficients)
  data.frame(M = if (is.na(M)) Mhat else M, M_multiple = if (is.na(M)) 1 else M / Mhat,
             est_pp = 100 * c1$estimate[1], se_pp = 100 * c1$std.error[1],
             max_bias_pp = if (!is.null(c1$maximum.bias)) 100 * c1$maximum.bias[1] else NA,
             lb = 100 * c1$conf.low[1], ub = 100 * c1$conf.high[1],
             ub_one_sided_95 = if (!is.null(c1$conf.high.onesided)) 100 * c1$conf.high.onesided[1] else NA,
             bandwidth = if (!is.null(c1$bandwidth)) c1$bandwidth[1] else NA)
}))
print(hon, digits = 3, row.names = FALSE); wcsv(hon, "10_rdhonest_M_grid.csv")
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\nDONE 10\n")
