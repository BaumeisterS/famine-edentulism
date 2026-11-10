# =====================================================================
#  05_outcome_threshold.R — outcome definition (CLHLS tooth-count threshold)
#
#  CLHLS records the number of remaining natural teeth (0-32). The pooled
#  data code CLHLS edentulism as <= 4 remaining teeth (>= 28 lost) and the
#  region DiD as 0 remaining; the paper's text describes >= 30 teeth lost
#  (<= 2 remaining). Here every design is re-run with the
#  CLHLS outcome re-coded as k = 0 (complete edentulism), k = 2 (the
#  paper's stated definition) and k = 4 (published coding), plus a
#  HARMONISED-OUTCOME analysis that drops CLHLS so that all nine remaining
#  studies share the identical binary self-report question ("lost all
#  teeth"). Optional rdrobust runs (published method) with RUN_RDROBUST=1.
# =====================================================================
source("01_common.R")
RUN_RDROBUST <- Sys.getenv("RUN_RDROBUST", "0") == "1"
if (RUN_RDROBUST) suppressPackageStartupMessages(library(rdrobust))

est_set <- function(d, studies = NULL, lab) {
  dd <- if (is.null(studies)) d else d[d$study %in% studies, ]
  one <- length(unique(dd$study)) == 1
  f <- list(
    RDiT       = make_rd(dd),
    RF         = make_rf(dd),
    ITS_local  = make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff"),
    Dose       = make_dose(dd, win = 15, df = 4))
  md <- mean(dd$dose15[dd$exposed15 == 1])
  do.call(rbind, lapply(names(f), function(nm) {
    r <- tryCatch(f[[nm]](dd$edent), error = function(e) c(est = NA, se = NA))
    sc <- if (nm == "Dose") md else 1
    data.frame(analysis = lab, estimator = nm, est_pp = 100 * sc * r["est"], se_pp = 100 * sc * r["se"])
  }))
}

out <- list(); prev <- list()
for (k in c(0, 2, 4)) {
  d <- make_analytic(load_pooled(k))
  lab <- sprintf("CLHLS edentulous = <=%d teeth", k)
  out[[length(out) + 1]] <- est_set(d, NULL, paste("Pooled 10 studies;", lab))
  out[[length(out) + 1]] <- est_set(d, 2, paste("CLHLS only;", lab))
  prev[[length(prev) + 1]] <- d[d$study == 2, ] |> group_by(cohort3) |>
    summarise(n = n(), prevalence_pct = 100 * mean(edent), .groups = "drop") |>
    mutate(threshold = k)
  if (RUN_RDROBUST) {
    cm <- model.matrix(~ male + factor(father_ses) + factor(study), data = d)[, -1]
    r <- rdrobust(d$edent, d$run, c = 0, p = 1, kernel = "triangular", bwselect = "mserd",
                  covs = cm, masspoints = "adjust")
    out[[length(out) + 1]] <- data.frame(analysis = paste("Pooled 10 studies;", lab),
      estimator = "rdrobust (published method)", est_pp = 100 * r$coef[1], se_pp = 100 * r$se[1])
  }
}
# harmonised outcome: the nine studies with the identical binary question
d4 <- make_analytic(load_pooled(4))
out[[length(out) + 1]] <- est_set(d4, setdiff(1:10, 2), "Harmonised outcome: 9 studies, binary 'lost all teeth' (CLHLS excluded)")
if (RUN_RDROBUST) {
  dh <- d4[d4$study != 2, ]
  cm <- model.matrix(~ male + factor(father_ses) + factor(study), data = dh)[, -1]
  r <- rdrobust(dh$edent, dh$run, c = 0, p = 1, kernel = "triangular", bwselect = "mserd",
                covs = cm, masspoints = "adjust")
  out[[length(out) + 1]] <- data.frame(analysis = "Harmonised outcome: 9 studies, binary 'lost all teeth' (CLHLS excluded)",
    estimator = "rdrobust (published method)", est_pp = 100 * r$coef[1], se_pp = 100 * r$se[1])
}
res <- do.call(rbind, out)
res$lb <- res$est_pp - zc * res$se_pp; res$ub <- res$est_pp + zc * res$se_pp
res$p  <- 2 * pnorm(-abs(res$est_pp / res$se_pp))
res$contains_plus2 <- res$lb <= 2 & res$ub >= 2
print(res, digits = 3, row.names = FALSE)
wcsv(res, if (RUN_RDROBUST) "05_outcome_threshold_with_rdrobust.csv" else "05_outcome_threshold.csv")
pv <- do.call(rbind, prev); print(pv); wcsv(pv, "05_clhls_prevalence_by_threshold.csv")
cat("DONE 05\n")
