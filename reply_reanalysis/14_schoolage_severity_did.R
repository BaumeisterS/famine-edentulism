# =====================================================================
#  14_schoolage_severity_did.R — SCHOOL-AGE famine exposure (CHARLS).
#  The RDiT compares exposure up to age 5 with exposure from age 6, so it
#  cannot detect harm from school-age exposure (Gu et al.'s largest OR;
#  the extended-susceptibility scenario). The graded severity DiD can:
#  cohorts born 1947-53 (aged 6-12 at the 1959 onset) vs the unexposed
#  post-famine cohorts born 1962-67 (1961 dropped: in utero/infant in the
#  final famine year), by 1960 provincial excess death rate (Lin & Yang
#  2000, Table 3). Birth-year and province fixed effects absorb the age
#  gap common to all provinces. Same contrast for the early-childhood
#  cohorts (1954-60) for comparison. Permutation-SD intervals.
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages(library(fixest))
set.seed(20261112)
ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning", Heilongjing = "Heilongjiang")
ly$prov <- ifelse(lab %in% names(std), std[lab], lab); map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
stopifnot(abs(ly$edr1960_abs[ly$prov == "Shanxi"] - 2.20) < 0.01, abs(ly$edr1960_abs[ly$prov == "Shaanxi"] - 1.90) < 0.01)
cd <- readRDS(file.path(LOC, "charls_did.rds")) |> filter(!is.na(edent), !is.na(age), !is.na(male), !is.na(birthprov)) |>
  mutate(prov = ifelse(birthprov %in% names(map), map[birthprov], birthprov)) |>
  left_join(ly[, c("prov", "edr1960_abs")], by = "prov") |> filter(!is.na(edr1960_abs))
Ih <- read.csv(file.path(OUT, "11_severity_implied_effects.csv"))$intensity[1]
run <- function(treated, control, lab) {
  s <- cd |> filter(byear %in% c(treated, control)) |> mutate(ex = as.integer(byear %in% treated), I = edr1960_abs / 10)
  f <- feols(edent ~ I:ex | prov + byear, data = s, cluster = ~prov); b <- coef(f)[1]
  pv <- unique(s[, c("prov", "I")])
  perm <- replicate(2000, { Ip <- setNames(sample(pv$I), pv$prov); s$Ip <- unname(Ip[s$prov])
    coef(feols(edent ~ Ip:ex | prov + byear, data = s, notes = FALSE))[1] })
  sp <- sd(perm)
  data.frame(contrast = lab, n = nobs(f), per10_pp = 100 * b, cr_se_pp = 100 * se(f)[1], perm_sd_pp = 100 * sp,
             lb = 100 * (b - zc * sp), ub = 100 * (b + zc * sp), p_perm = (1 + sum(abs(perm) >= abs(b))) / 2001,
             implied_high_mean_pp = 100 * b * Ih / 10, implied_lb = 100 * (b - zc * sp) * Ih / 10, implied_ub = 100 * (b + zc * sp) * Ih / 10)
}
res <- rbind(run(1947:1953, 1962:1967, "school-age exposure (born 1947-53, aged 6-12 at onset) vs post-famine (1962-67)"),
             run(1954:1960, 1962:1967, "early-childhood exposure (born 1954-60) vs post-famine (1962-67)"))
print(res, digits = 3, row.names = FALSE); wcsv(res, "14_schoolage_severity_did.csv")
cat("DONE 14\n")
