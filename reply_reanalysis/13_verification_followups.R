# =====================================================================
#  13_verification_followups.R — follow-up analyses for the reply
#  (1) equivalence (TOST, +/-2 pp) restated with robust and honest inference
#  (2) CHARLS binary severity DiD with permutation-SD and HC1 intervals;
#      graded severity DiD with permutation-SD intervals
#  (3) RD at the famine-END boundary of the ages 1-5 exposure (the edge
#      at the end of the famine): exposed (born <= famine end - 1) vs not
#      exposed at ages 1-5 (born >= famine end)
#  (4) relative scale in CHARLS (1.2% baseline)
#  (5) survivor-selection bias grid: masked effect = m x (r_d - p_s)
#  (6) CLHLS checks: one row per respondent; coding of missing tooth
#      counts in the DiD
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages({ library(fixest); library(haven) })
set.seed(20261111)
out <- list()

# (1) ---------------------------------------------------------------
rr <- read.csv(file.path(OUT, "10_rdrobust_inference.csv")); hn <- read.csv(file.path(OUT, "10_rdhonest_M_grid.csv"))[1, ]
tost <- function(est, se, m = 2) c(p_upper = pnorm((est - m) / se), p_lower = 1 - pnorm((est + m) / se))
t1 <- rbind(conventional = tost(rr$est_pp[1], rr$se_pp[1]), robust = tost(rr$est_pp[3], rr$se_pp[3]))
t1 <- data.frame(inference = rownames(t1), t1, tost_p = apply(t1, 1, max))
# honest: one-sided limits from RDHonest (bias-aware); equivalence needs both one-sided 95% limits inside +/-2
t1 <- rbind(t1, data.frame(inference = "honest (one-sided 95% limits)", p_upper = NA, p_lower = NA, tost_p = NA))
t1$note <- c("", "", sprintf("one-sided 95%% upper limit %+.2f; lower limit %+.2f (two-sided 90%% interval)", hn$ub_one_sided_95,
                             hn$est_pp - (hn$ub_one_sided_95 - hn$est_pp)))
print(t1, digits = 3); out$tost <- t1; wcsv(t1, "13_tost_restated.csv")

# (2) ---------------------------------------------------------------
did <- read.csv(file.path(OUT, "06_severity_did_inference.csv")); ch <- did[did$study == "CHARLS", ]
ci2 <- data.frame(control = ch$control, est = ch$est_unadj_pp,
  lb_cr = ch$est_unadj_pp - zc * ch$se_cr_pp, ub_cr = ch$est_unadj_pp + zc * ch$se_cr_pp,
  lb_hc1 = ch$est_unadj_pp - zc * ch$se_hc1_pp, ub_hc1 = ch$est_unadj_pp + zc * ch$se_hc1_pp,
  lb_perm = ch$est_unadj_pp - zc * ch$perm_sd_pp, ub_perm = ch$est_unadj_pp + zc * ch$perm_sd_pp)
print(ci2, digits = 3); wcsv(ci2, "13_charls_binary_did_intervals.csv")

ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning", Heilongjing = "Heilongjiang")
ly$prov <- ifelse(lab %in% names(std), std[lab], lab); map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
cd <- readRDS(file.path(LOC, "charls_did.rds")) |> filter(!is.na(edent), !is.na(age), !is.na(male), !is.na(birthprov)) |>
  mutate(prov = ifelse(birthprov %in% names(map), map[birthprov], birthprov)) |>
  left_join(ly[, c("prov", "edr1960_abs")], by = "prov") |> filter(!is.na(edr1960_abs))
Ih <- read.csv(file.path(OUT, "11_severity_implied_effects.csv"))$intensity[1]
g2 <- do.call(rbind, lapply(list(pre = 1947:1953, post = 1961:1967), function(cw) {
  s <- cd |> filter(byear %in% c(cw, 1954:1960)) |> mutate(ex = as.integer(byear %in% 1954:1960), I = edr1960_abs / 10)
  b <- coef(feols(edent ~ I:ex | prov + byear, data = s))[1]
  pv <- unique(s[, c("prov", "I")])
  perm <- replicate(2000, { Ip <- setNames(sample(pv$I), pv$prov); s$Ip <- unname(Ip[s$prov])
    coef(feols(edent ~ Ip:ex | prov + byear, data = s, notes = FALSE))[1] })
  data.frame(est_per10 = 100 * b, perm_sd = 100 * sd(perm), lb = 100 * (b - zc * sd(perm)), ub = 100 * (b + zc * sd(perm)),
             implied_high_mean = 100 * b * Ih / 10, implied_lb = 100 * (b - zc * sd(perm)) * Ih / 10,
             implied_ub = 100 * (b + zc * sd(perm)) * Ih / 10)
}))
g2$control <- c("pre-famine (1947-53)", "post-famine (1961-67)"); print(g2, digits = 3); wcsv(g2, "13_graded_did_perm_intervals.csv")

# (3) ---------------------------------------------------------------
d <- make_analytic(load_pooled(4))
de <- d; de$run <- de$byear - de$famine_end; de$D <- as.integer(de$run >= 0)   # right = not exposed at ages 1-5
fe <- make_rd(de)(de$edent)
fe2 <- make_rd(de[de$study != 2, ])(de$edent[de$study != 2])
end_rd <- data.frame(analysis = c("pooled 10 studies", "9 studies without CLHLS"),
                     jump_pp = 100 * c(fe["est"], fe2["est"]), se_pp = 100 * c(fe["se"], fe2["se"]))
end_rd$famine_effect_pp <- -end_rd$jump_pp      # harm predicts a DROP after the exposure window
end_rd$lb <- end_rd$famine_effect_pp - zc * end_rd$se_pp; end_rd$ub <- end_rd$famine_effect_pp + zc * end_rd$se_pp
print(end_rd, digits = 3); wcsv(end_rd, "13_famine_end_boundary_rd.csv")

# (4) ---------------------------------------------------------------
p0 <- mean(d$edent[d$study == 1 & d$cohort3 %in% 2])
cb <- read.csv(file.path(OUT, "12_charls_specific_bounds.csv"))
rel <- data.frame(bound = c(cb$design, "graded severity DiD (pre), high-mortality mean", "graded severity DiD (post), high-mortality mean"),
                  upper_pp = c(cb$ub, read.csv(file.path(OUT, "11_severity_implied_effects.csv"))$ub[1], g2$implied_ub[2]))
rel$baseline_pct <- 100 * p0; rel$rr_upper <- (p0 + rel$upper_pp / 100) / p0
print(rel, digits = 3); wcsv(rel, "13_charls_relative_scale.csv")

# (5) ---------------------------------------------------------------
grid <- expand.grid(extra_mortality_pct = c(2, 5, 10, 20), risk_gap_pp = c(5, 10, 20, 40))
grid$masked_effect_pp <- grid$extra_mortality_pct / 100 * grid$risk_gap_pp
print(reshape(grid, idvar = "extra_mortality_pct", timevar = "risk_gap_pp", direction = "wide"))
wcsv(grid, "13_survivor_bias_grid.csv")

# (6) ---------------------------------------------------------------
DAT <- Sys.getenv("FAMINE_DATA_DIR", unset = "c:/arbeit/Greifswald/Daten/")
cl <- read_dta(paste0(DAT, "CLHLS/clhls_2008_2018_longitudinal_dataset_released_version1.dta"),
               col_select = any_of(c("id", "v_bthyr", "g21", "a41", "trueage", "a1")))
s6 <- cl |> mutate(byear = as.numeric(v_bthyr)) |> filter(byear %in% 1931:1942, !is.na(trueage), a1 %in% 1:2, !is.na(a41))
chk <- data.frame(rows = nrow(cl), unique_ids = length(unique(cl$id)),
                  did_window_rows = nrow(s6),
                  did_window_missing_or_invalid_teeth = sum(is.na(s6$g21) | s6$g21 > 32))
print(chk); wcsv(chk, "13_clhls_checks.csv")
cat("DONE 13\n")
