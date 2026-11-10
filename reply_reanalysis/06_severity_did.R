# =====================================================================
#  06_severity_did.R — the region-of-birth x famine-INTENSITY (severity)
#  difference-in-differences, re-examined (published CLHLS estimate
#  +2.10 pp, p = 0.07) and results by famine severity.
#
#  CHARLS: 9 documented high-mortality birth provinces vs the rest
#          (Great Famine 1959-61); pre 1947-53, exposed 1954-60, post 1961-67.
#  CLHLS : Henan-born vs born elsewhere (Henan famine 1942-43);
#          pre 1931-36, exposed 1937-42, post 1943-48. ONE treated province.
#  Inference: (a) cluster-robust by province (as published); (b) individual-
#  level HC1 (captures the treated cells' own sampling noise, ignores
#  province shocks: a lower bound on uncertainty); (c) placebo-province
#  permutation inference (Conley & Taber 2011; MacKinnon & Webb 2020) —
#  raw and studentised; (d) event-study pre-trend Wald tests.
# =====================================================================
source("01_common.R")
suppressPackageStartupMessages(library(fixest))
set.seed(20261109)
HIGH <- c("Anhui","Sichuan","Chongqing","Guizhou","Gansu","Qinghai","Henan","Hunan","Guangxi")

did_fit <- function(s, adj) {
  s$exposed <- as.integer(s$exposed)
  f <- if (adj) lm(edent ~ treat * exposed + age + I(age^2) + I(age^3) + male, data = s) else
               lm(edent ~ treat * exposed, data = s)
  b  <- coef(f)["treat:exposed"]
  se_cr <- sqrt(vcovCL(f, cluster = s$region)["treat:exposed", "treat:exposed"])
  se_hc <- sqrt(vcovHC(f, type = "HC1")["treat:exposed", "treat:exposed"])
  c(est = unname(b), se_cr = se_cr, se_hc = se_hc)
}
did_raw <- function(y, tr, ex) {        # fast 2x2 DiD of cell means
  (mean(y[tr == 1 & ex == 1]) - mean(y[tr == 1 & ex == 0])) - (mean(y[tr == 0 & ex == 1]) - mean(y[tr == 0 & ex == 0]))
}
did_t <- function(y, tr, ex) {          # studentised by the binomial-cell SE
  v <- function(z) var(z) / length(z)
  dd <- did_raw(y, tr, ex)
  se <- sqrt(v(y[tr == 1 & ex == 1]) + v(y[tr == 1 & ex == 0]) + v(y[tr == 0 & ex == 1]) + v(y[tr == 0 & ex == 0]))
  c(did = dd, t = dd / se)
}

analyse <- function(dat, study, treated_set, ctrl = "pre", k_lab = "", n_perm = 5000, min_prov = 40) {
  s <- dat |> filter(.data[[ctrl]] | exposed, !is.na(edent), !is.na(age), !is.na(male), !is.na(region))
  s$treat <- as.integer(s$region %in% treated_set)
  rU <- did_fit(s, FALSE); rA <- did_fit(s, TRUE)
  obs <- did_t(s$edent, s$treat, s$exposed)
  # placebo provinces
  provs <- s |> count(region, exposed) |> group_by(region) |> filter(n() == 2, min(n) >= 10) |>
    summarise(n = sum(n), .groups = "drop") |> filter(n >= min_prov)
  if (length(treated_set) == 1) {                        # CLHLS: each other province as placebo
    pl <- setdiff(provs$region, treated_set)
    perm <- t(sapply(pl, function(r) did_t(s$edent, as.integer(s$region == r), s$exposed)))
  } else {                                              # CHARLS: random sets of the same size
    pool <- provs$region
    perm <- t(replicate(n_perm, { r <- sample(pool, length(treated_set))
      did_t(s$edent, as.integer(s$region %in% r), s$exposed) }))
  }
  pr <- function(stat, o) (1 + sum(abs(perm[, stat]) >= abs(o))) / (1 + nrow(perm))
  data.frame(study = study, outcome = k_lab, control = ctrl,
    n = nrow(s), n_treated_pre = sum(s$treat == 1 & !s$exposed), n_treated_exposed = sum(s$treat == 1 & s$exposed),
    est_unadj_pp = 100 * rU["est"], se_cr_pp = 100 * rU["se_cr"], se_hc1_pp = 100 * rU["se_hc"],
    p_cr = 2 * pnorm(-abs(rU["est"] / rU["se_cr"])), p_hc1 = 2 * pnorm(-abs(rU["est"] / rU["se_hc"])),
    est_adj_pp = 100 * rA["est"], se_adj_cr_pp = 100 * rA["se_cr"], se_adj_hc1_pp = 100 * rA["se_hc"],
    n_placebo = nrow(perm), perm_sd_pp = 100 * sd(perm[, "did"]),
    p_perm_raw = pr("did", obs["did"]), p_perm_studentised = pr("t", obs["t"]))
}

# ---------------- CHARLS ----------------
ch <- readRDS(file.path(LOC, "charls_did.rds")) |>
  transmute(edent, byear, age, male, region = birthprov,
            pre = byear %in% 1947:1953, exposed = byear %in% 1954:1960, post = byear %in% 1961:1967)
res <- list()
res[[1]] <- analyse(ch, "CHARLS", HIGH, "pre",  "self-report: lost all teeth")
res[[2]] <- analyse(ch, "CHARLS", HIGH, "post", "self-report: lost all teeth")

# ---------------- CLHLS ----------------
cl0 <- readRDS(file.path(LOC, "clhls_slim.rds"))
cat("CLHLS birth-province labels (n):\n"); print(sort(table(cl0$region), decreasing = TRUE)[1:12])
for (k in c(0, 2, 4)) {
  cl <- cl0 |> transmute(edent = ifelse(is.na(teeth), NA, as.numeric(teeth <= k)), byear, age, male, region,
                         pre = byear %in% 1931:1936, exposed = byear %in% 1937:1942, post = byear %in% 1943:1948)
  for (ctrl in c("pre", "post"))
    res[[length(res) + 1]] <- analyse(cl, "CLHLS", "henan", ctrl, sprintf("<= %d natural teeth", k))
}
tab <- do.call(rbind, res)
print(tab, digits = 3, row.names = FALSE)
wcsv(tab, "06_severity_did_inference.csv")

# ---------------- event-study pre-trend tests ----------------
pretrend <- function(dat, treated_set, gs, study) {
  e <- dat |> filter(byear >= gs - 7, byear <= gs + 13, !is.na(edent), !is.na(male), !is.na(region)) |>
    mutate(treat = as.integer(region %in% treated_set))
  m <- feols(edent ~ i(byear, treat, ref = gs - 1) + male | region + byear, data = e)
  nm <- names(coef(m)); pre <- nm[grepl("byear::", nm) & as.numeric(sub(":treat", "", sub("byear::", "", nm))) < gs - 1]
  out <- lapply(list(CR_province = vcov(m, cluster = ~region), HC1 = vcov(m, vcov = "hetero")), function(V) {
    b <- coef(m)[pre]; W <- as.numeric(t(b) %*% solve(V[pre, pre]) %*% b)
    c(wald = W, df = length(pre), p = pchisq(W, length(pre), lower.tail = FALSE)) })
  data.frame(study = study, vcov = names(out), wald = sapply(out, `[`, 1), df = sapply(out, `[`, 2),
             p = sapply(out, `[`, 3), n_pre_coefs = length(pre))
}
cl_k0 <- cl0 |> transmute(edent = ifelse(is.na(teeth), NA, as.numeric(teeth == 0)), byear, age, male, region)
pt <- rbind(pretrend(ch, HIGH, 1954, "CHARLS"), pretrend(cl_k0, "henan", 1937, "CLHLS (0 teeth)"))
print(pt, digits = 3, row.names = FALSE); wcsv(pt, "06_pretrend_tests.csv")

# ---------------- CLHLS Henan cell sizes ----------------
cs <- cl0 |> filter(byear >= 1930, byear <= 1950, !is.na(teeth)) |>
  mutate(henan = region == "henan") |> count(byear, henan) |>
  tidyr::pivot_wider(names_from = henan, values_from = n, names_prefix = "henan_")
wcsv(cs, "06_clhls_henan_cell_sizes.csv")
cat("DONE 06\n")
