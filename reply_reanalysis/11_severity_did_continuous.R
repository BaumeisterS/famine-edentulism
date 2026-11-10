# =====================================================================
#  11_severity_did_continuous.R — GRADED severity DiD (CHARLS, Great
#  Famine 1959-61). Replaces the published binary high/low province split
#  by a continuous provincial famine intensity: the 1960 excess death rate
#  (EDR1960 = 1960 crude death rate - mean 1956-58, per 1,000), computed
#  from the State Statistical Bureau provincial death rates printed in
#  Lin & Yang (2000, Econ J 110:136-158, Table 3, p. 147), the definition of
#  Chen & Zhou (2007, J Health Econ 26:659-681). Chongqing-born are given
#  Sichuan's value and Hainan-born Guangdong's (both were part of those
#  provinces in 1960); Tibet is not in the table.
#  Robustness: mean excess 1959-61; 1956-57 baseline (famine mortality
#  had begun in 1958 in Sichuan, Gansu, Yunnan, Guizhou, Ningxia); rank.
#  Inference: cluster-robust by province (28 clusters, intensity varies in
#  every cluster) and permutation of intensity across provinces.
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages(library(fixest))
set.seed(20261109)
ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
# The printed labels are kept in the parsed table and are easy to mis-merge:
# '(4) Shansi' = Shanxi (North region); '(24) Shanxi' = Shaanxi (North-West); 'Neimonggu' =
# Inner Mongolia; 'Lianing' = Liaoning; 'Heilongjing' = Heilongjiang; a
# control character follows the asterisked municipalities.
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning",
         Heilongjing = "Heilongjiang")
ly$province_LY <- ifelse(lab %in% names(std), std[lab], lab)
stopifnot(abs(ly$edr1960_abs[ly$province_LY == "Shanxi"]  - 2.20) < 0.01,   # values documented by the source check
          abs(ly$edr1960_abs[ly$province_LY == "Shaanxi"] - 1.90) < 0.01,
          abs(ly$edr1960_abs[ly$province_LY == "Anhui"]   - 56.70) < 0.01,
          length(unique(ly$province_LY)) == 28)
ly$edr_base5657 <- ly$`1960` - (ly$`1956` + ly$`1957`) / 2
ly$edr_mean_base5657 <- (ly$`1959` + ly$`1960` + ly$`1961`) / 3 - (ly$`1956` + ly$`1957`) / 2
cat("Lin & Yang provinces:", nrow(ly), "\n"); print(ly[order(-ly$edr1960_abs), c("province_LY", "edr1960_abs", "edr_mean5961_abs", "edr_base5657")][1:10, ])
map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
HIGH <- c("Anhui","Sichuan","Chongqing","Guizhou","Gansu","Qinghai","Henan","Hunan","Guangxi")

ch <- readRDS(file.path(LOC, "charls_did.rds")) |>
  filter(!is.na(edent), !is.na(age), !is.na(male), !is.na(birthprov)) |>
  mutate(prov_ly = ifelse(birthprov %in% names(map), map[birthprov], birthprov))
unm <- setdiff(unique(ch$prov_ly), ly$province_LY); cat("Birth provinces without a Lin & Yang value:", paste(unm, collapse = ", "), "\n")
ch <- ch |> left_join(ly |> transmute(prov_ly = province_LY, edr1960 = edr1960_abs, edr5961 = edr_mean5961_abs,
                                      edr5657 = edr_base5657, edrm5657 = edr_mean_base5657), by = "prov_ly") |>
  filter(!is.na(edr1960)) |>
  mutate(rank_pct = rank(edr1960) / n(), high = as.integer(birthprov %in% HIGH))
# province-level rank (not person-level)
pr <- ch |> distinct(prov_ly, edr1960) |> mutate(rank_pct = rank(edr1960) / n()); ch$rank_pct <- pr$rank_pct[match(ch$prov_ly, pr$prov_ly)]
cat(sprintf("CHARLS with intensity: n = %d in %d provinces; EDR1960 range %.1f to %.1f per 1,000\n",
            nrow(ch), length(unique(ch$prov_ly)), min(ch$edr1960), max(ch$edr1960)))
cat("Mean EDR1960, documented high-mortality vs other birth provinces (person-weighted):",
    round(tapply(ch$edr1960, ch$high, mean), 2), "\n")

win <- list(pre = 1947:1953, exposed = 1954:1960, post = 1961:1967)
fit_cont <- function(dat, ivar, ctrl = "pre", adj = FALSE, n_perm = 2000) {
  s <- dat |> filter(byear %in% c(win[[ctrl]], win$exposed)) |>
    mutate(ex = as.integer(byear %in% win$exposed), I = .data[[ivar]] / 10)   # per 10 excess deaths / 1,000
  f <- if (adj) feols(edent ~ I:ex + male + age + I(age^2) + I(age^3) | prov_ly + byear, data = s, cluster = ~prov_ly) else
                feols(edent ~ I:ex | prov_ly + byear, data = s, cluster = ~prov_ly)
  b <- coef(f)["I:ex"]; se <- se(f)["I:ex"]
  # permutation: shuffle province intensities across provinces
  pv <- unique(s[, c("prov_ly", "I")])
  perm <- replicate(n_perm, { Ip <- setNames(sample(pv$I), pv$prov_ly); s$Ip <- unname(Ip[s$prov_ly])
    coef(feols(edent ~ Ip:ex | prov_ly + byear, data = s, notes = FALSE))[1] })
  data.frame(intensity = ivar, control = ctrl, adjusted = adj, n = nobs(f),
             per10_pp = 100 * b, se_pp = 100 * se, lb = 100 * (b - zc * se), ub = 100 * (b + zc * se),
             p_cr = 2 * pnorm(-abs(b / se)), p_perm = (1 + sum(abs(perm) >= abs(b))) / (1 + n_perm))
}
res <- rbind(
  fit_cont(ch, "edr1960", "pre"), fit_cont(ch, "edr1960", "pre", adj = TRUE), fit_cont(ch, "edr1960", "post"),
  fit_cont(ch, "edr5961", "pre"), fit_cont(ch, "edr5657", "pre"), fit_cont(ch, "edrm5657", "pre"),
  fit_cont(ch, "rank_pct", "pre", n_perm = 2000))
res$per10_pp[res$intensity == "rank_pct"] <- res$per10_pp[res$intensity == "rank_pct"] / 10     # rank: per full range (0->1) of rank
res$note <- ifelse(res$intensity == "rank_pct", "effect for lowest -> highest-ranked province (rank 0 -> 1)", "pp per 10 excess deaths per 1,000")
for (v in c("se_pp", "lb", "ub")) res[[v]][res$intensity == "rank_pct"] <- res[[v]][res$intensity == "rank_pct"] / 10
print(res, digits = 3, row.names = FALSE); wcsv(res, "11_severity_did_continuous.csv")

# implied effect at the person-weighted mean intensity of the 9 documented high-mortality provinces and at Anhui
b <- res$per10_pp[1] / 10; s <- res$se_pp[1] / 10; Ih <- mean(ch$edr1960[ch$high == 1]); Il <- mean(ch$edr1960[ch$high == 0])
imp <- data.frame(contrast = c("high-mortality provinces (mean) vs zero excess", "high vs low (difference of means)", "Anhui (56.7) vs zero excess"),
                  intensity = c(Ih, Ih - Il, 56.7))
imp$effect_pp <- b * imp$intensity; imp$lb <- (b - zc * s) * imp$intensity; imp$ub <- (b + zc * s) * imp$intensity
print(imp, digits = 3, row.names = FALSE); wcsv(imp, "11_severity_implied_effects.csv")

# intensity terciles (non-parametric dose)
s3 <- ch |> filter(byear %in% c(win$pre, win$exposed)) |> mutate(ex = as.integer(byear %in% win$exposed),
  terc = cut(edr1960, quantile(unique(edr1960), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "middle", "high")))
ft <- feols(edent ~ i(terc, ex, ref = "low") | prov_ly + byear, data = s3, cluster = ~prov_ly)
tt <- as.data.frame(coeftable(ft)); tt$term <- rownames(tt); print(tt); wcsv(tt, "11_severity_terciles.csv")

# event study with continuous intensity (ref = 1953) and pre-trend Wald test
es <- ch |> filter(byear >= 1947, byear <= 1967) |> mutate(I = edr1960 / 10)
fe <- feols(edent ~ i(byear, I, ref = 1953) + male | prov_ly + byear, data = es, cluster = ~prov_ly)
ct <- as.data.frame(coeftable(fe)); ct$byear <- as.numeric(sub(":I", "", sub("byear::", "", rownames(ct))))
pre <- rownames(ct)[!is.na(ct$byear) & ct$byear < 1953]
W <- as.numeric(t(coef(fe)[pre]) %*% solve(vcov(fe)[pre, pre]) %*% coef(fe)[pre])
cat(sprintf("Continuous-intensity event study: pre-trend Wald = %.2f on %d df, p = %.3f\n", W, length(pre), pchisq(W, length(pre), lower.tail = FALSE)))
ct$period <- ifelse(ct$byear < 1954, "pre", ifelse(ct$byear <= 1960, "exposed", "post")); wcsv(ct, "11_severity_eventstudy.csv")
print(round(ct[, c("byear", "Estimate", "Std. Error")] * c(1, 100, 100), 2))
cat("DONE 11\n")
