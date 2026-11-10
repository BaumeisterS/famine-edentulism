# =====================================================================
#  03_its_specifications.R — the interrupted time series (ITS) re-examined
#
#  The published ITS (-4.01 pp, SE 2.78) is reproduced by an unweighted
#  fit of study x birth-year cell means (the published code evaluates n
#  after edent has been summarised, so every cell has n = 1). This script
#  (i) reproduces the published estimate, (ii) gives the
#  cell-size-weighted estimate, and (iii) shows how the level shift depends
#  on the specification of the secular birth-cohort trend (common linear
#  trend over 60-90 cohorts vs study-specific local trends).
# =====================================================================
source("01_common.R")
p <- load_pooled(4)
dpub <- p |> filter(!is.na(edent), !is.na(byear)) |> as.data.frame()   # published ITS sample

specs <- list(
  list(lab = "A. Published as coded: unweighted, common linear trend, all cohorts (gap at famine onset)",
       w = FALSE, win = NULL, tr = "common", deg = 1, ce = "famine_start"),
  list(lab = "B. Cell-size weighted, common linear trend, all cohorts (gap at famine onset)",
       w = TRUE,  win = NULL, tr = "common", deg = 1, ce = "famine_start"),
  list(lab = "C. Unweighted, common linear trend, all cohorts (jump at cut-off)",
       w = FALSE, win = NULL, tr = "common", deg = 1, ce = "cutoff"),
  list(lab = "D. Weighted, common linear trend, all cohorts (jump at cut-off)",
       w = TRUE,  win = NULL, tr = "common", deg = 1, ce = "cutoff"),
  list(lab = "E. Weighted, common QUADRATIC trend, all cohorts (jump at cut-off)",
       w = TRUE,  win = NULL, tr = "common", deg = 2, ce = "cutoff"),
  list(lab = "F. Weighted, study-specific linear trends, all cohorts (jump at cut-off)",
       w = TRUE,  win = NULL, tr = "study", deg = 1, ce = "cutoff"),
  list(lab = "G. Weighted, study-specific linear trends, +/-15 cohorts (jump at cut-off)",
       w = TRUE,  win = 15, tr = "study", deg = 1, ce = "cutoff"),
  list(lab = "H. Weighted, study-specific linear trends, +/-10 cohorts (jump at cut-off)",
       w = TRUE,  win = 10, tr = "study", deg = 1, ce = "cutoff"),
  list(lab = "I. Weighted, study-specific quadratic trends, +/-15 cohorts (jump at cut-off)",
       w = TRUE,  win = 15, tr = "study", deg = 2, ce = "cutoff"),
  list(lab = "J. Unweighted, study-specific linear trends, +/-15 cohorts (jump at cut-off)",
       w = FALSE, win = 15, tr = "study", deg = 1, ce = "cutoff"))

res <- do.call(rbind, lapply(specs, function(s) {
  f <- make_its(dpub, weighted = s$w, win = s$win, trend = s$tr, degree = s$deg, centre = s$ce)(dpub$edent)
  data.frame(spec = s$lab, est_pp = 100 * f["est"], se_pp = 100 * f["se"],
             lb = 100 * (f["est"] - zc * f["se"]), ub = 100 * (f["est"] + zc * f["se"]),
             p = 2 * pnorm(-abs(f["est"] / f["se"])))
}))
print(res, digits = 3, row.names = FALSE)
wcsv(res, "03_its_specifications.csv")

# study-specific local ITS (weighted, +/-15, linear) - heterogeneity of the jump
bys <- do.call(rbind, lapply(1:10, function(k) {
  dd <- dpub[dpub$study == k, ]
  f <- tryCatch(make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(dd$edent),
                error = function(e) c(est = NA, se = NA))
  data.frame(study = study_labels[k], est_pp = 100 * f["est"], se_pp = 100 * f["se"])
}))
print(bys, digits = 3, row.names = FALSE)
wcsv(bys, "03_its_local_by_study.csv")

# lack-of-fit of the common linear trend: cell-level residual curvature
key <- paste(dpub$study, dpub$byear); cm <- tapply(dpub$edent, key, mean); cn <- tapply(dpub$edent, key, length)
cd <- data.frame(key = names(cm), y = as.numeric(cm), n = as.numeric(cn))
cd$study <- as.numeric(sub(" .*", "", cd$key)); cd$byear <- as.numeric(sub(".* ", "", cd$key))
cd <- merge(cd, unique(dpub[, c("study", "famine_start")]), by = "study")
cd$time <- cd$byear - cd$famine_start; cd$post <- as.integer(cd$byear >= cd$famine_start - 5)
m1 <- lm(y ~ time + post + I(time * post) + factor(study), data = cd, weights = n)
m2 <- update(m1, . ~ . + factor(study):poly(time, 3))
cat("\nLack of fit of the common linear trend (weighted): F-test vs study-specific cubic trends\n")
print(anova(m1, m2))
cat("DONE 03\n")
