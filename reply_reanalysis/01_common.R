# =====================================================================
#  01_common.R — shared helpers for the re-analysis (sourced by 02-07)
#  Re-analysis of Baumeister et al. 2026, J Clin Periodontol 53:1673-80
#  (doi:10.1111/jcpe.70190) for the authors' reply to correspondence.
#  Run every script from this folder (Rscript 02_reproduce.R, ...).
#
#  All estimators are "precomputed closures": the design matrix, weights
#  and bread are built once from the data, and the returned function maps
#  an outcome vector y (aligned with the rows of the data) to c(est, se)
#  in milliseconds. This makes plasmode simulation over thousands of
#  replications feasible and guarantees that the observed estimate and the
#  simulated estimates come from the identical estimator.
# =====================================================================
suppressPackageStartupMessages({
  library(dplyr); library(splines); library(sandwich); library(lmtest)
})
RS  <- normalizePath(Sys.getenv("REANALYSIS_DIR", "."), winslash = "/")   # this folder
# GH = root of the published analysis repository (famine_edentulism_analysis.R, its tables/ and the
# pooled analysis file it writes); set FAMINE_REPO_DIR, or keep this folder inside that repository
GH  <- Sys.getenv("FAMINE_REPO_DIR", "")
if (!nzchar(GH)) for (cand in c(file.path(RS, ".."), file.path(RS, "..", "..", "github")))
  if (file.exists(file.path(cand, "famine_edentulism_analysis.R"))) { GH <- normalizePath(cand, winslash = "/"); break }
if (!nzchar(GH)) stop("Set FAMINE_REPO_DIR to the folder that holds famine_edentulism_analysis.R")
OUT <- file.path(RS, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOC <- file.path(RS, "data_local")
H_BW <- 9.52                       # MSE-optimal bandwidth of the published primary RDiT (years)
zc   <- qnorm(0.975)
study_labels <- c("CHARLS","CLHLS","LASI","SAGE-India","SAGE-Mexico","SAGE-Russia",
                  "SHARE-Netherlands","SHARE-Germany","SHARE-Austria","SHARE-Spain")
num <- function(x) as.numeric(haven::zap_labels(x))

# ---------------------------------------------------------------------
# Data: the published harmonised pooled file, optionally with the CLHLS
# outcome re-coded from the raw natural-tooth count (k = max remaining
# teeth counted as edentulous; k = 4 reproduces the published coding).
# ---------------------------------------------------------------------
build_clhls_rows <- function(k) {
  cl <- readRDS(file.path(LOC, "clhls_slim.rds"))
  cl |>
    transmute(study = 2, country = 1, continent = 1,
              edent = ifelse(is.na(teeth), NA_real_, as.numeric(teeth <= k)),
              treatment = case_when(byear %in% 1937:1942 ~ 1, byear %in% 1943:1948 ~ 0),
              byear, famine_start = 1942, famine_end = 1945, age, male = as.numeric(male),
              father_occupation, father_education) |>
    mutate(study_f = factor("CLHLS", levels = study_labels),
           famine_dur = famine_end - famine_start + 1,
           short_famine = as.integer(famine_dur <= 2),
           event_time = byear - famine_start,
           cohort3 = case_when(byear <  famine_start - 5 ~ 0,
                               byear >= famine_start - 5 & byear <= famine_end - 1 ~ 1,
                               byear >  famine_end ~ 2),
           father_ses = case_when(
             !is.na(father_education) & father_education == 1 ~ 1,
             !is.na(father_education) & father_education == 2 ~ 2,
             !is.na(father_education) & father_education == 3 ~ 3,
             !is.na(father_occupation) & father_occupation == 1 ~ 1,
             !is.na(father_occupation) & father_occupation == 0 ~ 2,
             TRUE ~ 1),
           age_c = age - 70, age_c2 = age_c^2, age_c3 = age_c^3)
}

load_pooled <- function(clhls_k = 4) {
  p <- readRDS(file.path(GH, "famine_edent_pooled.rds"))
  p <- p |> mutate(across(c(study, edent, treatment, byear, famine_start, famine_end, age, male,
                            cohort3, father_ses, age_c, age_c2, age_c3, famine_dur), num))
  if (!identical(clhls_k, 4)) {
    keep <- intersect(names(p), names(build_clhls_rows(clhls_k)))
    p <- bind_rows(p[p$study != 2, keep], build_clhls_rows(clhls_k)[, keep])
  }
  p
}

# famine calendar years [fs, fe] that fall in ages [a_lo, a_hi] of a child born in byear
overlap_years <- function(byear, fs, fe, a_lo, a_hi) {
  lo <- pmax(fs, byear + a_lo); hi <- pmin(fe, byear + a_hi); pmax(0, hi - lo + 1)
}

# analytic sample exactly as in the published script + design variables
make_analytic <- function(p) {
  p |> filter(!is.na(edent), !is.na(byear), !is.na(male)) |>
    mutate(run      = byear - (famine_start - 5),          # published RDiT running variable
           D        = as.integer(run >= 0),
           dose15   = overlap_years(byear, famine_start, famine_end, 1, 5),   # famine years at ages 1-5
           dose05   = overlap_years(byear, famine_start, famine_end, 0, 5),   # incl. infancy
           dose615  = overlap_years(byear, famine_start, famine_end, 6, 15),  # older-age exposure
           rec15    = overlap_years(byear, famine_end + 1, famine_end + 2, 1, 5), # 2 recovery years
           exposed15 = as.integer(dose15 >= 1),            # identical to cohort3 == 1
           cell     = paste(study, byear)) |>
    as.data.frame()
}

drop_collinear <- function(X, w = NULL) {
  Xq <- if (is.null(w)) X else X * sqrt(w)
  q <- qr(Xq, tol = 1e-9); keep <- sort(q$pivot[seq_len(q$rank)]); X[, keep, drop = FALSE]
}

# ---------------------------------------------------------------------
# 1. Sharp RDiT: local linear, triangular kernel, bandwidth h, covariates
#    (sex, early-life SEP, study FE) entered linearly as in rdrobust's
#    covariate-adjusted estimator; HC1 variance.
# ---------------------------------------------------------------------
make_rd <- function(d, h = H_BW, covs = TRUE, donut = 0) {
  idx <- if (donut > 0) which(abs(d$run) <= h & abs(d$run) > donut) else which(abs(d$run) <= h)
  s <- d[idx, ]; w <- pmax(0, 1 - abs(s$run) / h)
  X <- cbind(int = 1, D = s$D, run = s$run, Drun = s$D * s$run)
  if (covs) X <- cbind(X, cov_block(s))
  X <- drop_collinear(X, w)
  B <- solve(crossprod(X, X * w)); P <- B %*% t(X * w)
  n <- nrow(X); k <- ncol(X); j <- which(colnames(X) == "D")
  function(y) {
    yy <- y[idx]; b <- as.numeric(P %*% yy); e <- as.numeric(yy - X %*% b)
    V <- B %*% crossprod(X * (w * e)) %*% B * n / (n - k)
    c(est = b[j], se = sqrt(V[j, j]))
  }
}

# ---------------------------------------------------------------------
# 2. Reduced form: famine-exposed vs immediately post-famine cohort,
#    adjusted for a cubic in current age, sex, SEP, study FE (published
#    Table S2 '+ Age cubic' model). HC1 variance.
# ---------------------------------------------------------------------
make_rf <- function(d, studies = NULL) {
  in_s <- if (is.null(studies)) rep(TRUE, nrow(d)) else d$study %in% studies
  idx <- which(d$cohort3 %in% 1:2 & !is.na(d$treatment) & !is.na(d$age) & in_s)
  s <- d[idx, ]
  X <- drop_collinear(cbind(model.matrix(~ treatment + age_c + age_c2 + age_c3, data = s), cov_block(s)))
  B <- solve(crossprod(X)); P <- B %*% t(X)
  n <- nrow(X); k <- ncol(X); j <- which(colnames(X) == "treatment")
  function(y) {
    yy <- y[idx]; b <- as.numeric(P %*% yy); e <- as.numeric(yy - X %*% b)
    V <- B %*% crossprod(X * e) %*% B * n / (n - k)
    c(est = b[j], se = sqrt(V[j, j]))
  }
}

# ---------------------------------------------------------------------
# 3. Interrupted time series: birth-year cells by study, weighted by n,
#    level shift + slope change at the cut-off (famine start - 5), study FE,
#    HC3 variance (published specification).
# ---------------------------------------------------------------------
#  weighted = FALSE reproduces the published estimate (-4.01, SE 2.78):
#  the published code evaluates the cell size inside the same summarise()
#  call that replaces edent by its mean (summarise(edent = mean(edent),
#  n = sum(!is.na(edent)))), so every cell enters with n = 1 (unweighted
#  cell means). weighted = TRUE weights each cell by its number of
#  respondents. trend = "study" gives
#  study-specific trends and slope changes; win restricts to cells within
#  +/- win years of the cut-off (local ITS); degree = 2 adds curvature.
#  centre = "famine_start" reproduces the published parameterisation
#  (time = birth year - famine start, so 'post' is the gap between the
#  two segments at famine ONSET, 5 years after the cut-off); centre =
#  "cutoff" makes 'post' the jump at the cut-off itself.
make_its <- function(d, weighted = TRUE, win = NULL, trend = c("common", "study"), degree = 1,
                     centre = c("famine_start", "cutoff")) {
  trend <- match.arg(trend); centre <- match.arg(centre)
  key <- paste(d$study, d$byear); u <- unique(key); g <- match(key, u)
  first <- match(u, key)
  cd <- data.frame(study = d$study[first], byear = d$byear[first], fs = d$famine_start[first])
  nc <- tabulate(g, nbins = length(u))
  rel <- cd$byear - (cd$fs - 5)                  # years relative to the cut-off
  cd$time <- if (centre == "cutoff") rel else cd$byear - cd$fs
  cd$post <- as.integer(rel >= 0); cd$tp <- cd$time * cd$post
  use <- if (is.null(win)) rep(TRUE, nrow(cd)) else abs(rel) <= win
  if (length(unique(cd$study)) == 1) trend <- "single"
  f <- switch(trend,
    single = if (degree == 1) ~ time + post + tp else ~ time + I(time^2) + post + tp + I(tp * time),
    common = if (degree == 1) ~ time + post + tp + factor(study) else
             ~ time + I(time^2) + post + tp + I(tp * time) + factor(study),
    study  = if (degree == 1) ~ post + factor(study) + factor(study):time + factor(study):tp else
             ~ post + factor(study) + factor(study):time + factor(study):I(time^2) +
               factor(study):tp + factor(study):I(tp * time))
  wts <- if (weighted) nc else rep(1, length(nc))
  X <- drop_collinear(model.matrix(f, data = cd[use, ]), wts[use])
  w <- wts[use]; B <- solve(crossprod(X, X * w)); P <- B %*% t(X * w)
  hii <- w * rowSums((X %*% B) * X); j <- which(colnames(X) == "post")
  function(y) {
    ybar <- (as.numeric(rowsum(y, g, reorder = TRUE)) / nc)[use]
    b <- as.numeric(P %*% ybar); e <- ybar - as.numeric(X %*% b)
    V <- B %*% crossprod(X * (w * e / pmax(1 - hii, 1e-8))) %*% B
    c(est = b[j], se = sqrt(V[j, j]))
  }
}

# covariate block that silently drops factors with a single level (e.g. SEP in SAGE)
cov_block <- function(s, with_study = TRUE) {
  terms <- c("male", if (length(unique(s$father_ses)) > 1) "factor(father_ses)",
             if (with_study && length(unique(s$study)) > 1) "factor(study)")
  model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = s)[, -1, drop = FALSE]
}

# ---------------------------------------------------------------------
# 4. Dose-response (NEW): edentulism on the number of famine years that
#    fell in ages 1-5 (graded exposure), with a study-specific natural
#    spline of birth year (absorbs the smooth age/cohort gradient), sex,
#    SEP and study FE. Identification comes from the tent-shaped dose
#    profile across birth years, which a smooth monotone age gradient
#    cannot mimic. Window: famine start - win to famine end + win.
#    Cluster-robust (CR1) by study x birth-year cell (exposure is
#    assigned at the cohort level).
# ---------------------------------------------------------------------
make_dose <- function(d, win = 15, df = 4, dose_var = "dose15", extra = character(0),
                      studies = NULL) {
  in_s <- if (is.null(studies)) rep(TRUE, nrow(d)) else d$study %in% studies
  idx <- which(d$byear >= d$famine_start - win & d$byear <= d$famine_end + win & in_s)
  s <- d[idx, ]; us <- sort(unique(s$study))
  S <- do.call(cbind, lapply(us, function(k) {
    M <- matrix(0, nrow(s), df); r <- s$study == k
    M[r, ] <- ns(s$byear[r], df = df); colnames(M) <- paste0("ns", k, "_", seq_len(df)); M }))
  base <- cbind(int = 1, dose = s[[dose_var]])
  if (length(extra)) base <- cbind(base, as.matrix(s[, extra, drop = FALSE]))
  X <- drop_collinear(cbind(base, S, cov_block(s)))
  B <- solve(crossprod(X)); P <- B %*% t(X)
  n <- nrow(X); k <- ncol(X); j <- which(colnames(X) == "dose")
  cl <- match(s$cell, unique(s$cell)); G <- max(cl)
  adj <- G / (G - 1) * (n - 1) / (n - k)
  function(y) {
    yy <- y[idx]; b <- as.numeric(P %*% yy); e <- as.numeric(yy - X %*% b)
    V <- B %*% crossprod(rowsum(X * e, cl, reorder = FALSE)) %*% B * adj
    c(est = b[j], se = sqrt(V[j, j]))
  }
}

# helpers ---------------------------------------------------------------
tost_p <- function(est, se, m) max(pnorm((est - m) / se), 1 - pnorm((est + m) / se))
fmt_ci <- function(est, se, digits = 2) sprintf(paste0("%+.", digits, "f (%+.", digits, "f, %+.", digits, "f)"),
                                               est, est - zc * se, est + zc * se)
wcsv <- function(x, f) write.csv(x, file.path(OUT, f), row.names = FALSE)
