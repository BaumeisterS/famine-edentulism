# =====================================================================
#  15_final_checks.R — final checks for the reply
#  (a) honest RDiT on a study- and covariate-adjusted outcome (the
#      published RDHonest fit has no study FE or covariates)
#  (b) graded severity DiD (CHARLS) with inference allowing cohort-level
#      shocks: two-way clustering (province, birth year) and a
#      cohort-aggregated estimator (one slope per birth year)
#  (c) agreement test (J statistic) among the four GLS inputs, per shape
#  (d) CLHLS DiD coding check saved to file (published coding vs missing
#      tooth counts excluded); exposed vs mean of pre- and post-famine
#      cohorts (removes a Henan-specific linear cohort trend); age at
#      interview of the Henan-born post-famine respondents
#  (e) duration strata with the CLHLS famine dated 1942-43 (short famine)
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages({ library(fixest); library(RDHonest); library(sandwich) })
set.seed(20261113)
sink(file.path(OUT, "15_final_checks.log"), split = TRUE)

# (a) ---------------------------------------------------------------
d <- make_analytic(load_pooled(4))
d$y_dm  <- d$edent - ave(d$edent, d$study)
d$y_res <- resid(lm(edent ~ male + factor(father_ses) + factor(study), data = d, na.action = na.exclude))
hon <- do.call(rbind, lapply(c("y_dm", "y_res"), function(v) {
  dd <- data.frame(y = d[[v]], run = d$run); dd <- dd[is.finite(dd$y), ]
  f <- RDHonest(y ~ run, data = dd, cutoff = 0, kern = "triangular", opt.criterion = "MSE", sclass = "H")
  c1 <- as.data.frame(f$coefficients)
  data.frame(outcome = c(y_dm = "study-demeaned", y_res = "residualised on sex, SES, study")[v],
             est_pp = 100 * c1$estimate[1], lb = 100 * c1$conf.low[1], ub = 100 * c1$conf.high[1],
             ub_one_sided_95 = 100 * c1$conf.high.onesided[1], bandwidth = c1$bandwidth[1], M = c1$M[1])
}))
print(hon, digits = 3, row.names = FALSE); wcsv(hon, "15_rdhonest_adjusted.csv")

# (b) ---------------------------------------------------------------
ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning", Heilongjing = "Heilongjiang")
ly$prov <- ifelse(lab %in% names(std), std[lab], lab); map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
stopifnot(abs(ly$edr1960_abs[ly$prov == "Shanxi"] - 2.20) < 0.01, abs(ly$edr1960_abs[ly$prov == "Shaanxi"] - 1.90) < 0.01)
cd <- readRDS(file.path(LOC, "charls_did.rds")) |> filter(!is.na(edent), !is.na(age), !is.na(male), !is.na(birthprov)) |>
  mutate(prov = ifelse(birthprov %in% names(map), map[birthprov], birthprov)) |>
  left_join(ly[, c("prov", "edr1960_abs")], by = "prov") |> filter(!is.na(edr1960_abs)) |> mutate(I = edr1960_abs / 10)
Ih <- read.csv(file.path(OUT, "11_severity_implied_effects.csv"))$intensity[1]
cohort_slopes <- function(s) {     # one intensity slope per birth year (province-level variation within cohort)
  t(sapply(sort(unique(s$byear)), function(y) { m <- lm(edent ~ I, data = s[s$byear == y, ]); c(byear = y, b = coef(m)[2]) }))
}
gd <- do.call(rbind, lapply(list(
    list("exposed 1954-60 vs pre-famine 1947-53", 1954:1960, 1947:1953),
    list("exposed 1954-60 vs post-famine 1961-67", 1954:1960, 1961:1967),
    list("exposed 1954-60 vs post-famine 1962-67", 1954:1960, 1962:1967),
    list("school-age 1947-53 vs post-famine 1962-67", 1947:1953, 1962:1967)), function(L) {
  s <- cd |> filter(byear %in% c(L[[2]], L[[3]])) |> mutate(ex = as.integer(byear %in% L[[2]]))
  f <- feols(edent ~ I:ex | prov + byear, data = s)
  b <- coef(f)[1]; se1 <- se(summary(f, cluster = ~prov))[1]; se2 <- se(summary(f, cluster = ~prov + byear))[1]
  cs <- as.data.frame(cohort_slopes(s)); names(cs) <- c("byear", "b"); cs$ex <- as.integer(cs$byear %in% L[[2]])
  ca <- lm(b ~ ex, data = cs); bca <- coef(ca)[2]; seca <- sqrt(vcov(ca)[2, 2]); dfca <- ca$df.residual
  data.frame(contrast = L[[1]], per10_pp = 100 * b, se_prov_pp = 100 * se1, se_twoway_pp = 100 * se2,
             lb_twoway = 100 * (b - zc * se2), ub_twoway = 100 * (b + zc * se2),
             cohort_agg_pp = 100 * bca, cohort_agg_se_pp = 100 * seca, cohort_df = dfca,
             lb_cohort = 100 * (bca - qt(0.975, dfca) * seca), ub_cohort = 100 * (bca + qt(0.975, dfca) * seca),
             implied_high_twoway = sprintf("%+.2f (%+.2f, %+.2f)", 100 * b * Ih / 10, 100 * (b - zc * se2) * Ih / 10, 100 * (b + zc * se2) * Ih / 10),
             implied_high_cohort = sprintf("%+.2f (%+.2f, %+.2f)", 100 * bca * Ih / 10, 100 * (bca - qt(0.975, dfca) * seca) * Ih / 10,
                                           100 * (bca + qt(0.975, dfca) * seca) * Ih / 10))
}))
print(gd, digits = 3, row.names = FALSE); wcsv(gd, "15_graded_did_cohort_inference.csv")

# (c) ---------------------------------------------------------------
z <- readRDS(file.path(OUT, "04_plasmode_raw.rds"))
valid <- c("RDiT", "RF", "ITS_local", "Dose")
b0 <- colMeans(z$null[, valid, "est"]); Sig <- cov(z$null[, valid, "est"]); Si <- solve(Sig)
th <- z$obs["est", valid] - b0
g1 <- read.csv(file.path(OUT, "04_plasmode_gls.csv")); g2 <- read.csv(file.path(OUT, "04c_susceptibility_gls.csv"))
A <- rbind(as.matrix(g1[, c("a_RDiT", "a_RF", "a_ITS_local", "a_Dose")]), as.matrix(g2[, c("a_RDiT", "a_RF", "a_ITS_local", "a_Dose")]))
jt <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
  a <- A[i, ]; v <- 1 / as.numeric(t(a) %*% Si %*% a); est <- v * as.numeric(t(a) %*% Si %*% th)
  r <- th - a * est; J <- as.numeric(t(r) %*% Si %*% r); w <- v * as.numeric(t(a) %*% Si)
  data.frame(shape = c(g1$shape, g2$shape)[i], combined_pp = 100 * est, J = J, df = 3, p = 1 - pchisq(J, 3),
             w_RDiT = w[1], w_RF = w[2], w_ITS_local = w[3], w_Dose = w[4])
}))
cat("\nGLS weights (sum to 1 after scaling by the response) and J test of agreement:\n")
print(jt, digits = 3, row.names = FALSE); wcsv(jt, "15_gls_agreement_weights.csv")

# (d) ---------------------------------------------------------------
cl <- readRDS(file.path(LOC, "clhls_slim.rds"))
cod <- do.call(rbind, lapply(c("published coding (missing count = dentate)", "missing tooth count excluded"), function(lab) {
  s <- cl |> mutate(edent = if (grepl("published", lab)) as.integer(!is.na(g21) & g21 == 0) else ifelse(is.na(teeth), NA, as.integer(teeth == 0)),
                    pre = byear %in% 1931:1936, exposed = byear %in% 1937:1942) |>
    filter(pre | exposed, !is.na(edent), !is.na(age), !is.na(male), !is.na(region)) |>
    mutate(treat = as.integer(region == "henan"), exposed = as.integer(exposed))
  f <- lm(edent ~ treat * exposed, data = s)
  data.frame(coding = lab, did_pp = 100 * coef(f)["treat:exposed"],
             se_cr_pp = 100 * sqrt(vcovCL(f, cluster = s$region)["treat:exposed", "treat:exposed"]), n = nrow(s),
             missing_tooth_count_in_window = sum(is.na(cl$teeth[cl$byear %in% 1931:1942 & !is.na(cl$age) & !is.na(cl$male) & !is.na(cl$region)])))
}))
print(cod, digits = 3, row.names = FALSE); wcsv(cod, "15_clhls_did_coding.csv")

avg <- do.call(rbind, lapply(c(0, 2, 4), function(k) {
  s <- cl |> mutate(edent = ifelse(is.na(teeth), NA, as.numeric(teeth <= k)),
                    grp = ifelse(byear %in% 1931:1936, "pre", ifelse(byear %in% 1937:1942, "exposed", ifelse(byear %in% 1943:1948, "post", NA)))) |>
    filter(!is.na(grp), !is.na(edent), !is.na(age), !is.na(male), !is.na(region)) |>
    mutate(treat = as.integer(region == "henan"), grp = factor(grp, levels = c("exposed", "pre", "post")))
  f <- lm(edent ~ treat * grp, data = s); V <- vcovHC(f, type = "HC1")
  L <- setNames(rep(0, length(coef(f))), names(coef(f))); L[c("treat:grppre", "treat:grppost")] <- -0.5
  est <- sum(L * coef(f)); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  data.frame(threshold = sprintf("<= %d teeth", k), exposed_minus_mean_pre_post_pp = 100 * est, se_hc1_pp = 100 * se,
             lb = 100 * (est - zc * se), ub = 100 * (est + zc * se), n_henan_post = sum(s$treat == 1 & s$grp == "post"))
}))
print(avg, digits = 3, row.names = FALSE); wcsv(avg, "15_clhls_exposed_vs_mean_pre_post.csv")

hp <- cl |> filter(region == "henan", byear %in% 1943:1948, !is.na(teeth), !is.na(age), !is.na(male))
ages <- data.frame(n = nrow(hp), age_min = min(hp$age), age_median = median(hp$age), age_max = max(hp$age),
                   share_below_65 = mean(hp$age < 65), born_1943 = sum(hp$byear == 1943))
print(ages); wcsv(ages, "15_clhls_henan_post_ages.csv")

# (e) ---------------------------------------------------------------
p2 <- load_pooled(4); i <- p2$study == 2
p2$famine_end[i] <- 1943; p2$famine_dur[i] <- 2
p2$cohort3[i] <- with(p2[i, ], ifelse(byear < famine_start - 5, 0, ifelse(byear <= famine_end - 1, 1, ifelse(byear > famine_end, 2, NA))))
d2 <- make_analytic(p2)
dur <- tapply(d2$famine_dur, d2$study, `[`, 1)
strat <- do.call(rbind, lapply(list(Short = as.integer(names(dur)[dur <= 2]), Long = as.integer(names(dur)[dur > 2])), function(ix) {
  dd <- d2[d2$study %in% ix, ]; md <- mean(dd$dose15[dd$exposed15 == 1])
  r <- list(RDiT = make_rd(dd)(dd$edent), RF = make_rf(dd)(dd$edent),
            ITS_local = make_its(dd, weighted = TRUE, win = 15, trend = "study", centre = "cutoff")(dd$edent),
            Dose = make_dose(dd)(dd$edent) * md)
  data.frame(estimator = names(r), est_pp = 100 * sapply(r, `[`, 1), se_pp = 100 * sapply(r, `[`, 2),
             studies = paste(study_labels[ix], collapse = ", "))
}))
strat$stratum <- rep(c("Short (<= 2 y), CLHLS 1942-43", "Long (> 2 y)"), each = 4)
strat$lb <- strat$est_pp - zc * strat$se_pp; strat$ub <- strat$est_pp + zc * strat$se_pp
print(strat[, c("stratum", "estimator", "est_pp", "se_pp", "lb", "ub", "studies")], digits = 3, row.names = FALSE)
wcsv(strat, "15_duration_strata_clhls_short.csv")
cat("DONE 15\n"); sink()
