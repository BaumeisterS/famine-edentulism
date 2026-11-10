# =====================================================================
#  18h_its_placebo_cutoffs.R — placebo cut-offs for the local ITS (as 18a
#  for the RDiT): famine onset shifted so that the +/-15-cohort window
#  around the placebo cut-off excludes the famine-exposed cohorts.
# =====================================================================
source("01_common.R")
d <- make_analytic(load_pooled(4))
pi <- do.call(rbind, lapply(c(-40:-25, 24:40), function(cc) {
  dd <- d; dd$famine_start <- d$famine_start + cc
  r <- tryCatch(make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(dd$edent), error = function(e) c(est = NA, se = NA))
  data.frame(shift = cc, est_pp = 100 * r["est"], se_pp = 100 * r["se"])
}))
pi$t <- pi$est_pp / pi$se_pp; pi <- pi[is.finite(pi$t), ]
print(pi, digits = 3, row.names = FALSE)
s <- data.frame(n_placebos = nrow(pi), sd_est_pp = sd(pi$est_pp), mean_se_pp = mean(pi$se_pp), sd_t = sd(pi$t), share_abs_t_gt_1.96 = mean(abs(pi$t) > zc))
print(s, digits = 3); wcsv(pi, "18h_its_placebo_cutoffs.csv"); wcsv(s, "18h_its_placebo_cutoffs_summary.csv"); cat("DONE 18h\n")
