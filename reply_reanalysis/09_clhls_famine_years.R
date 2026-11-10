# =====================================================================
#  09_clhls_famine_years.R — the pooled data date the CLHLS (Henan) famine
#  1942-1945, the paper's Table 1 gives 1942-1943. Re-estimate the pooled
#  designs with the documented end year (affects famine-exposed / post-
#  famine cohort membership and dose; not the RDiT cut-off, which depends
#  on the famine start only).
# =====================================================================
source("01_common.R")
p  <- load_pooled(4)
p2 <- p; i <- p2$study == 2
p2$famine_end[i] <- 1943
p2$famine_dur[i] <- 2
p2$cohort3[i] <- with(p2[i, ], ifelse(byear < famine_start - 5, 0, ifelse(byear <= famine_end - 1, 1, ifelse(byear > famine_end, 2, NA))))
res <- do.call(rbind, lapply(list(`CLHLS famine 1942-1945 (pooled data)` = p, `CLHLS famine 1942-1943 (Table 1)` = p2), function(pp) {
  d <- make_analytic(pp); md <- mean(d$dose15[d$exposed15 == 1])
  r <- list(RDiT = make_rd(d)(d$edent), RF = make_rf(d)(d$edent),
            ITS_local = make_its(d, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(d$edent),
            Dose = make_dose(d)(d$edent) * md)
  data.frame(estimator = names(r), est_pp = 100 * sapply(r, `[`, 1), se_pp = 100 * sapply(r, `[`, 2))
}))
res$version <- rep(c("CLHLS famine 1942-1945 (pooled data)", "CLHLS famine 1942-1943 (Table 1)"), each = 4)
res$lb <- res$est_pp - zc * res$se_pp; res$ub <- res$est_pp + zc * res$se_pp
print(res[, c("version", "estimator", "est_pp", "lb", "ub")], digits = 3, row.names = FALSE)
wcsv(res, "09_clhls_famine_years.csv"); cat("DONE 09\n")
