# =====================================================================
#  02_reproduce.R — validate the re-analysis estimators against the
#  published results before any new analysis relies on them, and add the
#  observed dose-response (graded-exposure) estimates.
# =====================================================================
source("01_common.R")
set.seed(1)
p  <- load_pooled(4)
d  <- make_analytic(p)
cat("Analytic N =", nrow(d), "(published 165,586)\n")

# ---- (a) rebuilt CLHLS rows (k = 4) must equal the published pooled rows ----
key <- function(x) with(x, paste(byear, round(age), male, edent, treatment, cohort3, father_ses))
pub <- p[p$study == 2, ]; reb <- build_clhls_rows(4)
same <- identical(sort(key(pub)), sort(key(reb)))
cat("Rebuilt CLHLS rows identical to published pooled CLHLS rows:", same,
    "| n =", nrow(pub), "vs", nrow(reb), "\n")

# ---- (b) estimators vs published ----
rd   <- make_rd(d)(d$edent)
rd0  <- make_rd(d, covs = FALSE)(d$edent)
rf   <- make_rf(d)(d$edent)
its  <- make_its(d)(d$edent)
d_its_pub <- p |> filter(!is.na(edent), !is.na(byear)) |> mutate(run = 0) |> as.data.frame()
its_pub <- make_its(d_its_pub)(d_its_pub$edent)
n_rf <- sum(d$cohort3 %in% 1:2 & !is.na(d$treatment) & !is.na(d$age))

rep <- data.frame(
  quantity = c("RDiT local-linear h=9.52, covariates (fast)", "RDiT, no covariates (fast)",
               "Reduced form, age cubic (LPM)", "ITS level shift (analytic sample)",
               "ITS level shift (published sample: non-missing edent & birth year)"),
  est_pp   = 100 * c(rd["est"], rd0["est"], rf["est"], its["est"], its_pub["est"]),
  se_pp    = 100 * c(rd["se"],  rd0["se"],  rf["se"],  its["se"],  its_pub["se"]),
  published = c("+0.134 (SE 0.833) rdrobust conventional", "—",
                "+0.82 (SE 0.82), n = 23,697", "—", "-4.01 (SE 2.78)"),
  n = c(NA, NA, n_rf, NA, NA))
print(rep, digits = 3)
wcsv(rep, "02_reproduction.csv")

# ---- (c) dose-response: famine years at ages 1-5 (graded exposure) ----
cat("\nDose distribution among famine-exposed (aged 1-5 in a famine year):\n")
print(table(study = d$study[d$exposed15 == 1], dose = d$dose15[d$exposed15 == 1]))
mean_dose <- mean(d$dose15[d$exposed15 == 1])
cat(sprintf("Mean famine years at ages 1-5 among the exposed: %.2f\n", mean_dose))

specs <- expand.grid(win = c(10, 15, 20), df = c(3, 4, 5, 6), stringsAsFactors = FALSE)
dr <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  f <- make_dose(d, win = specs$win[i], df = specs$df[i])(d$edent)
  data.frame(model = "dose ages 1-5", win = specs$win[i], df = specs$df[i],
             per_year_pp = 100 * f["est"], se_pp = 100 * f["se"])
}))
alt <- rbind(
  data.frame(model = "dose ages 0-5 (incl. infancy)", win = 15, df = 4,
             per_year_pp = 100 * make_dose(d, dose_var = "dose05")(d$edent)["est"],
             se_pp = 100 * make_dose(d, dose_var = "dose05")(d$edent)["se"]),
  data.frame(model = "dose ages 1-5, adjusted for dose at ages 6-15", win = 15, df = 4,
             per_year_pp = 100 * make_dose(d, extra = "dose615")(d$edent)["est"],
             se_pp = 100 * make_dose(d, extra = "dose615")(d$edent)["se"]),
  data.frame(model = "dose ages 6-15 (older-age exposure; negative-control window)", win = 15, df = 4,
             per_year_pp = 100 * make_dose(d, dose_var = "dose615")(d$edent)["est"],
             se_pp = 100 * make_dose(d, dose_var = "dose615")(d$edent)["se"]))
dr <- rbind(dr, alt)
dr$lb <- dr$per_year_pp - zc * dr$se_pp; dr$ub <- dr$per_year_pp + zc * dr$se_pp
dr$p  <- 2 * pnorm(-abs(dr$per_year_pp / dr$se_pp))
dr$implied_mean_exposed_pp <- dr$per_year_pp * mean_dose     # effect at the exposed cohort's mean dose
dr$implied_ub_pp <- dr$ub * mean_dose
cat("\nDose-response (pp per famine year at ages 1-5; implied effect at the exposed mean dose):\n")
print(dr, digits = 3, row.names = FALSE)
wcsv(dr, "02_dose_response_pooled.csv")

# by study (main spec)
ds <- do.call(rbind, lapply(1:10, function(k) {
  f <- tryCatch(make_dose(d, studies = k)(d$edent), error = function(e) c(est = NA, se = NA))
  data.frame(study = k, label = study_labels[k], per_year_pp = 100 * f["est"], se_pp = 100 * f["se"],
             mean_dose_exposed = mean(d$dose15[d$exposed15 == 1 & d$study == k]))
}))
ds$lb <- ds$per_year_pp - zc * ds$se_pp; ds$ub <- ds$per_year_pp + zc * ds$se_pp
print(ds, digits = 3, row.names = FALSE)
wcsv(ds, "02_dose_response_by_study.csv")
saveRDS(list(mean_dose = mean_dose, rd = rd, rf = rf, its = its), file.path(OUT, "02_observed.rds"))
cat("DONE 02\n")
