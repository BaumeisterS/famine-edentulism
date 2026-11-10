# =====================================================================
#  12_extra_checks.R — additional checks
#  (1) Gu et al. ORs re-expressed at CHARLS's own post-famine prevalence
#      (Table S6 used the pooled 8.6% baseline) vs CHARLS-specific estimates
#  (2) information share of each study in the pooled RDiT; pooled RDiT
#      without CHARLS
#  (3) continuous-intensity event study: reference-independent pre-trend
#      test (equality of the pre-period coefficients; linear trend) and the
#      exposed-minus-pre contrast of the event-study coefficients
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages(library(fixest))
d <- make_analytic(load_pooled(4))

# (1) ---------------------------------------------------------------
p_ch <- mean(d$edent[d$study == 1 & d$cohort3 %in% 2]); p_pool <- mean(d$edent[d$cohort3 %in% 2 & d$treatment %in% 0])
or2rd <- function(or, p0) { o <- or * p0 / (1 - p0); o / (1 + o) - p0 }
gu <- data.frame(exposure_group = c("fetal", "preschool", "school-age"), adj_OR = c(2.097, 2.625, 3.907))
gu$rd_at_pooled_baseline_pp <- 100 * or2rd(gu$adj_OR, 0.086)
gu$rd_at_charls_baseline_pp <- 100 * or2rd(gu$adj_OR, p_ch)
cat(sprintf("CHARLS post-famine prevalence %.2f%% (pooled RF reference %.2f%%)\n", 100 * p_ch, 100 * p_pool))
print(gu, digits = 3, row.names = FALSE); wcsv(gu, "12_gu_benchmark_charls_baseline.csv")
ss <- read.csv(file.path(OUT, "07_study_specific_estimates.csv")); c1 <- ss[ss$study == 1, ]
chs <- data.frame(design = c("RDiT", "ITS (correctly specified)", "Dose-response (mean dose)"),
                  est = c(c1$RDiT_est, c1$ITS_est, c1$Dose_est), se = c(c1$RDiT_se, c1$ITS_se, c1$Dose_se))
chs$ub <- chs$est + zc * chs$se; print(chs, digits = 3, row.names = FALSE); wcsv(chs, "12_charls_specific_bounds.csv")

# (2) ---------------------------------------------------------------
w <- 1 / ss$RDiT_se^2; share <- data.frame(study = ss$label, info_share_pct = 100 * w / sum(w))
print(share[order(-share$info_share_pct), ], digits = 3, row.names = FALSE); wcsv(share, "12_rdit_information_share.csv")
dx <- d[d$study != 1, ]; r9 <- make_rd(dx)(dx$edent)
cat(sprintf("Pooled RDiT without CHARLS: %+.2f pp (95%% CI %+.2f to %+.2f)\n", 100 * r9["est"],
            100 * (r9["est"] - zc * r9["se"]), 100 * (r9["est"] + zc * r9["se"])))
wcsv(data.frame(analysis = "pooled RDiT, 9 studies without CHARLS", est_pp = 100 * r9["est"], se_pp = 100 * r9["se"]),
     "12_rdit_without_charls.csv")

# (3) ---------------------------------------------------------------
ly <- read.csv(file.path(RS, "inputs", "linyang2000_table3_province_death_rates.csv"), check.names = FALSE)
lab <- gsub("[^A-Za-z ]", "", ly$province_LY)
std <- c(Shansi = "Shanxi", Shanxi = "Shaanxi", Neimonggu = "Inner Mongolia", Lianing = "Liaoning", Heilongjing = "Heilongjiang")
ly$prov <- ifelse(lab %in% names(std), std[lab], lab)
map <- c(Chongqing = "Sichuan", Hainan = "Guangdong")
ch <- readRDS(file.path(LOC, "charls_did.rds")) |> filter(!is.na(edent), !is.na(male), !is.na(birthprov)) |>
  mutate(prov = ifelse(birthprov %in% names(map), map[birthprov], birthprov)) |>
  left_join(ly[, c("prov", "edr1960_abs")], by = "prov") |> filter(!is.na(edr1960_abs), byear >= 1947, byear <= 1967) |>
  mutate(I = edr1960_abs / 10)
fe <- feols(edent ~ i(byear, I, ref = 1953) + male | prov + byear, data = ch, cluster = ~prov)
b <- coef(fe); V <- vcov(fe); nm <- names(b)[grepl("byear::", names(b))]
yr <- as.numeric(sub(":I", "", sub("byear::", "", nm))); pre <- nm[yr < 1953]; exl <- nm[yr >= 1954 & yr <= 1960]
C <- t(sapply(pre[-1], function(k) { v <- setNames(rep(0, length(b)), names(b)); v[k] <- 1; v[pre[1]] <- -1; v }))
W_eq <- as.numeric(t(C %*% b) %*% solve(C %*% V %*% t(C)) %*% (C %*% b))
x <- yr[yr < 1953]; cl <- setNames(rep(0, length(b)), names(b)); cl[pre] <- (x - mean(x)) / sum((x - mean(x))^2)
tr <- sum(cl * b); tr_se <- sqrt(as.numeric(t(cl) %*% V %*% cl))
ce <- setNames(rep(0, length(b)), names(b)); ce[exl] <- 1 / length(exl); ce[pre] <- -1 / length(pre)
dd <- sum(ce * b); dd_se <- sqrt(as.numeric(t(ce) %*% V %*% ce))
out <- data.frame(test = c("pre-period coefficients equal (reference-free)", "linear pre-trend slope (pp per 10 deaths/1,000 per birth year)",
                           "exposed minus pre-period mean (pp per 10 deaths/1,000)"),
                  stat = c(W_eq, 100 * tr, 100 * dd), se = c(NA, 100 * tr_se, 100 * dd_se),
                  p = c(pchisq(W_eq, nrow(C), lower.tail = FALSE), 2 * pnorm(-abs(tr / tr_se)), 2 * pnorm(-abs(dd / dd_se))))
print(out, digits = 3, row.names = FALSE); wcsv(out, "12_severity_eventstudy_tests.csv")
cat("DONE 12\n")
